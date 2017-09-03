#####################################################
# Network-based linear regression by GBLasso.
# Wei Pan, 3/27/08
#####################################################

## squared error loss:
Xb<-function(x, b){
sum(x*b)
}

L<-function(Y,X,b){
r<-Y-apply(X, 1, Xb, b)
sum(r*r)
}

## network-based penalty value:
T<-function(b, netwk, gamma, wt){
T0<-0
K<-nrow(netwk)
if (gamma>=1){
  for(k in 1:K)
    T0<-T0 + ( (abs(b[netwk[k,1]])^gamma)/wt[netwk[k,1]] +
               (abs(b[netwk[k,2]])^gamma)/wt[netwk[k,2]] )^(1/gamma)
  } else{
  for(k in 1:K)
    T0<-T0 + max( abs(b[netwk[k,1]]/wt[netwk[k,1]]),  
                  abs(b[netwk[k,2]]/wt[netwk[k,2]])) 
  }
T0
}

deltaT<-function(b, dbj, j, netwk, gamma, wt){
dT0<-0
K<-nrow(netwk)
K1<-(1:K)[netwk[,1]==j]
K2<-(1:K)[netwk[,2]==j]
if (gamma>=1){
  for(k in K1){
    dT0<-dT0 + ( ((abs(b[netwk[k,1]]+dbj))^gamma)/wt[netwk[k,1]] +
               (abs(b[netwk[k,2]])^gamma)/wt[netwk[k,2]] )^(1/gamma) -
               ( (abs(b[netwk[k,1]])^gamma)/wt[netwk[k,1]] +
               (abs(b[netwk[k,2]])^gamma)/wt[netwk[k,2]] )^(1/gamma)
    }
  for(k in K2){
    dT0<-dT0 + ( (abs(b[netwk[k,1]])^gamma)/wt[netwk[k,1]] +
                ((abs(b[netwk[k,2]]+dbj))^gamma)/wt[netwk[k,2]] )^(1/gamma) -
               ( (abs(b[netwk[k,1]])^gamma)/wt[netwk[k,1]] +
               (abs(b[netwk[k,2]])^gamma)/wt[netwk[k,2]] )^(1/gamma)
    }
  } else{
  for(k in K1)
    dT0<-dT0 + max( abs((b[netwk[k,1]]+dbj)/wt[netwk[k,1]]),  
                    abs(b[netwk[k,2]]/wt[netwk[k,2]])) -
             max( abs(b[netwk[k,1]]/wt[netwk[k,1]]),
                  abs(b[netwk[k,2]]/wt[netwk[k,2]]))
  for(k in K2)
    dT0<-dT0 + max( abs(b[netwk[k,1]]/wt[netwk[k,1]]),  
                    abs((b[netwk[k,2]]+dbj)/wt[netwk[k,2]])) -
             max( abs(b[netwk[k,1]]/wt[netwk[k,1]]),
                  abs(b[netwk[k,2]]/wt[netwk[k,2]]))
  }
dT0
}

##GBLasso for least squares loss and our network-based penalty:
##Input:
##       Y: response; X: covariates;
##       netwk: adjacency matrix based on a network; only two columns;
##       gamma: gamma-norm used in the penalty; ANY gamma < 1 for L_infty!
##       wt: weight used for the norm;
GBLasso<-function(Y,X, netwk, gamma=2, wt, epsilon=0.1, ksi=0.01, MAXITER=1000){

n<-length(Y);
p<-length(X[1,]);

b0<-rep(0, p)
etaX<-rep(0, p)
for(j in 1:p)
  etaX[j]<-sum(Y*X[,j])

#step 1:
b1<-b0
while(T(b1, netwk, gamma, wt)==T(b0, netwk, gamma, wt)){
  b0<-b1
  jhat<-which.max(abs(etaX))
  if (length(jhat)>1) jhat<-jhat[1]
  shat<-sign(etaX[jhat])*epsilon
  b1[jhat]<-b0[jhat] + shat
  for(j in 1:p)
    etaX[j]<-etaX[j]- shat*sum(X[,jhat]*X[,j])
  }
lambda0<-(L(Y,X, b0) - L(Y,X, b1))/
         (T(b1, netwk, gamma, wt) - T(b0, netwk, gamma, wt))

betas<-matrix(0, nrow=MAXITER, ncol=p)
lambdas<-rep(0,MAXITER) 

   nsol<-1
   betas[nsol,]<-b0
   lambdas[nsol]<-lambda0

#steps 2-3:
lambda1<-lambda0
while (lambda1>0 && nsol<MAXITER){
b0<-b1
lambda0<-lambda1

s<-epsilon
dT<-rep(0, p)
for(j in 1:p){
  b01<-b0
  b01[j]<-b01[j]+s
  dT[j]<-deltaT(b0, s, j, netwk, gamma, wt) 
#cat("j=", j, " deltaT=", dT[j], " but dT=", dT[j]<-T(b01, netwk, gamma, wt) - T(b0, netwk, gamma, wt), "\n")
  }
dGamma1<-0-2*s*etaX + s*s*(n-1) + lambda0*dT

s<-0-epsilon
dT<-rep(0, p)
b01<-b0
for(j in 1:p){
  b01<-b0
  b01[j]<-b01[j]+s
  dT[j]<-deltaT(b0, s, j, netwk, gamma, wt) 
#cat("j=", j," deltaT=", dT[j], " But dT=", dT[j]<-T(b01, netwk, gamma, wt) - T(b0, netwk, gamma, wt), "\n")
  }
dGamma2<-0-2*s*etaX + s*s*(n-1) + lambda0*dT

j1hat<-which.min(dGamma1)
j2hat<-which.min(dGamma2)
if (dGamma1[j1hat]<=dGamma2[j2hat]) {
   jhat<-j1hat
   dGamma<-dGamma1[j1hat]
   shat<-epsilon
   } else {
   jhat<-j2hat
   dGamma<-dGamma2[j2hat]
   shat<-0-epsilon
   }

if (dGamma < 0-ksi){
   b1<-b0
   b1[jhat]<-b1[jhat] + shat
   lambda1<-lambda0
   for(j in 1:p)
     etaX[j]<-etaX[j]- shat*sum(X[,jhat]*X[,j])
#cat("F: ", jhat, shat, dGamma, lambda0,"\n")
   
   } else{
   #record current sol:
   if (lambda0<lambdas[nsol]){
     nsol<-nsol+1
     betas[nsol,]<-b0
     lambdas[nsol]<-lambda0
     }
 
   jhat<-which.max(abs(etaX))
   if (length(jhat)>1) jhat<-jhat[1]
   # below could result a dead loop b/w F<->B; e.g. when sim=3, wt=1
   #shat<-sign(etaX[jhat])*epsilon
   # thus modify the step size so that F & B will NOT overlap:
   shat<-sign(etaX[jhat])*epsilon*sqrt(2)
   b1<-b0
   b1[jhat]<-b0[jhat] + shat
   for(j in 1:p)
     etaX[j]<-etaX[j]- shat*sum(X[,jhat]*X[,j])
   
   lambda0<-(L(Y,X, b0) - L(Y,X, b1))/
            (T(b1, netwk, gamma, wt) - T(b0, netwk, gamma, wt))
   lambda1<-min(lambda1, lambda0)
#cat("B: ", jhat, shat, lambda0, lambda1, "\n")
   }
} #while
list(betas=betas[1:nsol,], lambdas=lambdas[1:nsol], nsol=nsol)

}
   
