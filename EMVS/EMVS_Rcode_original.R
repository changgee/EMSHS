##########################################################################
# **************** Auxiliary functions********************************** #
##########################################################################

# ******************** E-step ****************************************** #

# Mean field approximation to obtain expectation of MRF(mu,Sigma)

MFA<-function(start,mu,Sigma){

	E<-start
	Eps<-1
	E_new<-E
	p=nrow(Sigma)
	while (Eps>0.001){
		for ( i in 1:p )
		{
			E[i]<-mu[i]+sum(Sigma[,i]*E)
			E[i]<-1/(1+exp(-E[i]))
		}
		Eps<-sum(abs(E_new-E))
		E_new<-E
			    }
	E_new
}

# Find a transition point in the MRF prior within a sequence of sparsity parameters "mus" and 
# smoothness matrix "Sigma"

find_transition<-function(mus,Sigma,start){
	
	expect<-numeric(length(mus))
	for (i in (1:length(mus))){
		expect[i]<-sum(MFA(rep(start,ncol(Sigma)),rep(mus[i],ncol(Sigma)),Sigma))
			}
	index<-which.max(expect[-1]-expect[-length(expect)])
	plot(expect~mus,pch=19,ylab="Sum of E(gamma)",xlab=expression(theta[1]))
	abline(v=mus[index+1],lwd=2,col=2,lty=2)
	list(expect=expect,transition=mus[index])

}

# E-step for MRF

# t     ... temperature for deterministic annealing
# mu    ... sparsity parameter theta_1 in MRF prior
# Sigma ... smoothness matrix theta_2 in MRF prior


E_v0v1_MRF<-function(beta_k,sigma_k,v0,v1,mu,Sigma,t){

	mu_star<-as.numeric(t*1/2*log(v0/v1)-t*(v0-v1)/(2*sigma_k^2*v1*v0)*beta_k^2+t*mu)
      p_star<-MFA(rep(0.5,length(beta_k)),as.numeric(mu_star),t*Sigma)
	list(inv_var=p_star/v1+(1-p_star)/v0,prob=p_star)		

}


# E-step for beta-binomial and fixed prior

E_v0v1<-function(beta_k,sigma_k,v0,v1,p,t){
	
	p_star<-(p*dnorm(beta_k,0,sd=sigma_k*sqrt(v1)))^t/((p*dnorm(beta_k,0,sd=sigma_k*sqrt(v1)))^t+((1-p)*dnorm(beta_k,0,sd=sigma_k*sqrt(v0)))^t)
	list(inv_var=p_star/v1+(1-p_star)/v0,prob=p_star)		

}

# ******************** M-step ****************************************** #


# (1) UPDATE BETA

# inv_var ... conditional posterior averaged spike and slab precisions
# p>n     ... Woodburry Sherman
# p<=n    ... classical ridge computation

M_beta<-function(tXY,X,tXX,inv_var){

    n = nrow(X) 
      if (ncol(X)>n){
     
         X_Psinv<-t(t(X)/(tausq<-as.numeric(inv_var)))           
         (1/tausq*diag(length(tausq))-t(X_Psinv)%*%solve(diag(n)+X_Psinv%*%t(X),tol=1e-400)%*%X_Psinv)%*%tXY
        
       
                      } else
              {
            Psi<-as.numeric(inv_var)*diag(length(inv_var))
            solve(tXX+Psi,tol=1e-400)%*%tXY  
              } 
}


# (2) UPDATE SIGMA

# eta,lambda... parameters of the IG prior

M_sigma<-function(Y,X,beta_k,inv_var,eta,lambda){
   
        e<-Y-X%*%beta_k
        n<-length(Y)
        p<-dim(X)[2]
     
        as.numeric(sqrt((t(e)%*%e+t(beta_k*as.numeric(inv_var))%*%beta_k+lambda)/(n+p+eta)) )
}

# (3) UPDATE THETA

# beta-binomial prior

M_p<-function(post,a,b){
	(sum(post)+a-1)/(b+a+length(post)-2)
}

# logistic regression product prior

# thetas assigned the inverse logistic beta distribution with parameters a,b
# Z is the group identification matrix

logpi<-function(x,a,b){
	sum(a*x-(a+b)*log(1+exp(x)))
}

Q_logistic<-function(x,Z,probs,a,b){
	-(sum(probs*Z%*%x)-sum(log(1+exp(Z%*%x)))+sum(logpi(x,a,b)))
}


M_theta<-function(b_k,Z,post,a,b){
	o<-optim(b_k,Q_logistic,probs=post,a=a,b=b,Z=Z)
	o$par
}


# (4) UPDATE v1


M_v1<-function(v1,beta,sigma,post,a,b){
	o<-optimize(f_3,c(0,10e4),beta=beta,sigma=sigma,post=post,a=a,b=b)
	o$minimum
	
}



# *************** Posterior model evaluation *************************** #

# logarithm of the prior on the model space

log_prior<-function(gamma,type,probs,a,b,Sigma){
  
  return<-numeric(0)
  if (type=="fixed"){
    if (sum(probs==0.5)==length(probs)){
      return<-sum(gamma)*log(0.5)}else{
        return<-sum(gamma*log(probs)+(1-gamma)*log(1-probs))  }
  }
  if (type=="betabinomial"){
    
    x<-sum(gamma)+a
    y<-length(gamma)-sum(gamma)+b
    return<-log(beta(x,y))-log(beta(a,b))
    
    # Stirling approximation
    if (return=="-Inf"){
      return=0.5*log(2*pi)+(x-0.5)*log(x)+(y-0.5)*log(y)-(x+y-0.5)*log(x+y)
    }	
  }
  if (type=="MRF"){
    return<-t(gamma)%*%Sigma%*%gamma*beta(a+sum(gamma),b-sum(gamma))/beta(a,b)
  }
  return
  
}


# logarithm of the g-function

log_g<-function(gamma,X,Y,eta,lambda,v0=0,v1,type,probs,a,b,Sigma){

	q<-sum(gamma)
	if(q>0){
		n<-length(Y)
		X_gamma<-data.matrix(X[,gamma==1])
		diag<-rep(v1,sum(gamma))
		X_tilde<-rbind(X_gamma,sqrt(1/diag)*diag(q))
		Ssq<-t(Y)%*%Y-t(Y)%*%X_gamma%*%solve(t(X_gamma)%*%X_gamma+1/diag*diag(q))%*%t(X_gamma)%*%Y
		log<- -0.5*log(det(t(X_tilde)%*%X_tilde))-
			q/2*log(v1)-
		(n+eta)/2*log(eta*lambda+Ssq)+
		log_prior(gamma=gamma,type=type,probs=probs,a=a,b=b,Sigma=Sigma)
 		log} else {
		Ssq<-t(Y)%*%Y
		log<- -(n+eta)/2*log(eta*lambda+Ssq)+
		log_prior(gamma=gamma,type=type,probs=probs,a=a,b=b,Sigma=Sigma)
 		log
		}
}

##########################################################################
# **************** EMVS PROCEDURE ************************************** #
##########################################################################

# Y          ... n x 1 vector of centered responses
# X          ... n x p matrix of standardized predictors
# v0s        ... sequence of spike variance parameters
# v1         ... slab variance parameter; if missing then imposed a prior and estimated
# beta_init  ... starting values for the beta vector; if missing then the limiting
#		     case of deterministic annealing used
# sigma_init ... a starting value for sigma
# epsilon    ... convergence margin parameter
# type 	 ... type of the prior on the model space
#                type="betabinomial"  for the betabinomial prior with hyperparameters a and b
#	           type="MRF"           for the MRF prior  with hyperparameters theta and Sigma
#                type="fixed"         for the fixed inclusion probability p
#		     type="logistic"	  for the logistic regression prior with a grouping matrix Z
# temperature... inverse temperature parameter for annealing; if missing then
#		     t=1 without annealing
# mu         ... sparsity MRF hyper-parameter, where mu(1,...1)' is the mean vector of the MRF prior
# Sigma      ... smoothness matrix in MRF prior
# p          ... fixed prior probability of inclusion, if type="fixed"			
# a,b        ... hyperparameters for the betabinomial prior, if type="betabinomial"
# a_v1,b_v1  ... parameters of the prior for v1, if v1 missing
# v1_g       ... v1 value for the model evaluation by the g-function
# Z          ... group identification matrix for the logistic prior
# MRF implemented here with fixed sparsity parameter mu


EMVS<-function(Y,X,v0s,v1,type,
               beta_init,sigma_init,epsilon,
               temperature,Z,
               mu,Sigma,p,a,b,a_v1,b_v1,v1_g){

	expit<-function(x){exp(x)/(1+exp(x))}
	L<-length(v0s)            # number of steps in the regularization
	dim<-dim(X)[2]		  # number of predictors
    n<-nrow(X) 

	intersects<-numeric(L)    # intersection points between posterior weighted spike and slab
	sigmas<-numeric(L)	  # vector of sigma MAP estimates for each spike	
	betas<-matrix(0,L,dim)    # L x p matrix of MAP beta estimates for each spike	
	posts<-matrix(0,L,dim)    # L x p matrix of conditional posterior inclusion probabilities
	
      log_post<-numeric(L)      # logarithm of the g-function models associated with v0s
      thetas<-numeric(L)
	if (missing(v1)){v1s<-numeric(L)}else {v1s=rep(v1,L)}    

	if (missing(temperature)){temperature=1}
	
	niters<-numeric(L)	  # number of iterations



	tXY<-t(X)%*%Y
	tXX<-t(X)%*%X

	for (i in (1:L)){         

		v0<-v0s[i]

		# initialization with starting values

		if (missing(beta_init)){
                 beta_k<-M_beta(tXY,X,tXX,rep((v0s[i]+v1)/(2*v0s[i]*v1),dim))}else
		            {beta_k<-beta_init}

		beta_new<-beta_k
		sigma_k<-sigma_init

	      if (missing(v1)){v_k<-1000} else {v_k<-v1}    # if v1 missing, starting value for v1 is 1000
		if (type=="betabinomial"){p_k<-0.5}            # starting value 0.5 for the inclusion probability 
		if (type=="fixed"){p_k<-p}
		if (type=="logistic"){theta_k<-rep(0,ncol(Z))}
	      eps<-epsilon+1
		niter<-1
	
	while(eps>epsilon){

		print(paste("v0=",round(v0,2),"; iter=",niter,sep=""))
	
		
		if (type=="betabinomial"|type=="fixed"){
			E_step<-E_v0v1(beta_k,sigma_k,v0,v_k,p_k,temperature)			
						}
		if (type=="MRF"){
			E_step<-E_v0v1_MRF(beta_k,sigma_k,v0,v_k,mu=rep(mu,ncol(Sigma)),Sigma,temperature)
					}
	
		if (type=="logistic"){
			linpred<-Z%*%theta_k
			E_step<-E_v0v1(beta_k,sigma_k,v0,v_k,expit(linpred),temperature)	
					   }


		inv_var<-E_step$inv_var 		# conditional posterior weighted spike and slab precisions 
		post<-E_step$prob				# conditional inclusion probabilities
		beta_k<-M_beta(tXY,X,tXX,inv_var)

		sigma_k<-M_sigma(Y,X,beta_k,inv_var,1,1)

		if (type=="betabinomial"){p_k<-M_p(post,a=a,b=b)}
		if (type=="logistic"){theta_k<-M_theta(theta_k,Z,post,a=a,b=b)}
		
		if (missing(v1)){v_k<-max(10*v0,M_v1(v_k,beta_k,sigma_k,post,a=a_v1,b=b_v1))} # v1 has to be at least 10times greater than v0
		
		eps<-max(abs(beta_new-beta_k))
		beta_new<-beta_k
		print(eps)
		niter<-niter+1
				}


	posts[i,]<-post
	betas[i,]<-beta_new
	sigmas[i]<-sigma_k
    v1s[i]<-v_k
	c<-sqrt(v1s[i]/v0s[i])
	if (type=="betabinomial"|type=="fixed"){
		thetas[i]<-p_k
	    w<-(1-p_k)/p_k	
     		if (w!=0){
		intersects[i]<-sigmas[i]*sqrt(v0s[i])*sqrt(2*log(w*c)*c^2/(c^2-1))}else{
		intersects[i]=0}
			  
								}
	index<-post>0.5
	
#	log_post[i]<-as.numeric(log_g(as.numeric(index),X,Y,1,1,v0=0,v1=v1_g,type="betabinomial",a=a,b=b))

      niters[i]<-niter	      	
    
	if (type=="logistic"){

      list<-list(betas=betas,log_post=log_post,intersects=intersects,sigmas=sigmas,v1s=v1s,v0s=v0s,
	           niters=niters,posts=posts,inv_var=inv_var,type=type,theta=theta_k)
				
				   }

	if (type=="betabinomial"){

      list<-list(betas=betas,log_post=log_post,intersects=intersects,sigmas=sigmas,v1s=v1s,v0s=v0s,
	           niters=niters,posts=posts,inv_var=inv_var,type=type,thetas=thetas)
				
				   }

	if (type=="fixed"){

      list<-list(betas=betas,log_post=log_post,intersects=intersects,sigmas=sigmas,v1s=v1s,v0s=v0s,
	           niters=niters,posts=posts,inv_var=inv_var,type=type,p_k=p_k)
				
				   }

	if (type=="MRF"){

      list<-list(betas=betas,log_post=log_post,intersects=intersects,sigmas=sigmas,v1s=v1s,v0s=v0s,
	           niters=niters,posts=posts,inv_var=inv_var,type=type)
				
				   }




    
}

	list

}

##########################################################################
# ****************** REGULARIZATION PLOT ******************************* #		
##########################################################################

reg_plot<-function(result,plot_type){

	betas<-result$betas
	intersects<-result$intersects
	log_post<-result$log_post
	v0s<-result$v0s

	if (plot_type=="both"){
	par(mfrow=c(1,2))
	
	color<-apply(result$posts,2,function(x){as.numeric(x>0.5)})

					    
	plot(betas[,1]~v0s,type="l",pch=19,ylim=c(min(betas),max(betas)),
		xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2)
	for (i in (2:dim(betas)[2])){
		points(betas[,i]~v0s,type="l",pch=19,ylim=c(min(betas),max(betas)),
		col="grey",lty=2,lwd=1)	
					    }
	points(betas[,1]~v0s,pch=19,ylim=c(min(betas),max(betas)),
		xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,
		col=3*color[,1]+1,cex=0.8)
	for (i in (2:dim(betas)[2])){
		points(betas[,i]~v0s,pch=19,ylim=c(min(betas),max(betas)),
		col=3*color[,i]+1,lwd=1,cex=0.8)	
						}
      if (result$type=="betabinomial"|result$type=="fixed"){
	points(intersects~v0s,type="l",col=2,lwd=2,lty=1)
	points(-intersects~v0s,type="l",col=2,lwd=2,lty=1)
				   }
	abline(h=0,lwd=2,lty=2,col="lightgrey")	
	title("EMVS Regularization Plot")
	text(max(v0s),betas[length(v0s),],labels=paste("x",1:ncol(betas),sep=""),col=4)


	plot(as.numeric(log_post)~v0s,pch=4,type="b",col=2,lwd=1,xlab=expression(v[0]),ylab=
	expression(log(g(gamma))))

	title(expression(Log(g(gamma))))}
	
	if (plot_type=="reg"){
	color<-apply(result$posts,2,function(x){as.numeric(x>0.5)})

					    
	plot(betas[,1]~v0s,type="l",pch=19,ylim=c(min(betas),max(betas)),
		xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2)
	for (i in (2:dim(betas)[2])){
		points(betas[,i]~v0s,type="l",pch=19,ylim=c(min(betas),max(betas)),
		col="grey",lty=2,lwd=1)	
					    }
	points(betas[,1]~v0s,pch=19,ylim=c(min(betas),max(betas)),
		xlab=expression(v[0]),ylab=expression(hat(beta)),lwd=1,
		col=3*color[,1]+1,cex=0.8)
	for (i in (2:dim(betas)[2])){
		points(betas[,i]~v0s,pch=19,ylim=c(min(betas),max(betas)),
		col=3*color[,i]+1,lwd=1,cex=0.8)	
						}
      if (result$type=="betabinomial"|result$type=="fixed"){
	points(intersects~v0s,type="l",col=2,lwd=2,lty=1)
	points(-intersects~v0s,type="l",col=2,lwd=2,lty=1)
				   }
	abline(h=0,lwd=2,lty=2,col="lightgrey")	
	title("EMVS Regularization Plot")
	text(max(v0s),betas[length(v0s),],labels=paste("x",1:ncol(betas),sep=""),col=3,cex=1.2)


}
	if (plot_type=="gf"){


	plot(as.numeric(log_post)~v0s,pch=4,type="b",col=2,lwd=1,xlab=expression(v[0]),ylab=
	expression(log(g(gamma))))

	title(expression(Log(g(gamma))))

}



}

# returns indices of variables included in the local median probability model with the highest g-function

best_model<-function(result){
	
	which<-which.max(result$log_post)
	(1:dim(result$betas)[2])[result$posts[which,]>0.5]


}




