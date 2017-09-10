# seed: seed value for random numbers
# p: # of genes
# g: # of pathways
# gsm: expected size of a pathway
# ntrain: # of samples for training
# ntune: # of samples for tuning
# ntest: # of samples for testing
# sigma2: error variance
# edii: edge density between important genes in G0
# ediu: edge density between important and unimportant genes in G0
# eduu: edge density between unimportant genes in G0
# Gmode: the mode of G
#    0 - G=G0
#    1 - completely random graph
#    2 - G=G0 but edges between important and unimportant genes are thresholded
# Gthres: the threshold when Gmode=2
# savefile: file to save data

DG_batch <- function(R,head,seed,p,g,gsm,ntrain,ntune,ntest,sigma2,edii=0.5,ediu=0.05,eduu=0.2,Gmode=0,Gthres=0.1,batch=0)
{
  for ( i in 1:R )
  {
    fname = sprintf("%s/data%03d",head,batch+i)
    if ( !file.exists(fname) )
      DG(seed+batch+i,p,g,gsm,ntrain,ntune,ntest,sigma2,edii,ediu,eduu,Gmode,Gthres,fname)
  }
}

DG <- function(seed,p,g,gsm,ntrain,ntune,ntest,sigma2,edii=0.5,ediu=0.05,eduu=0.2,Gmode=0,Gthres=0.05,savefile=NULL)
{
  set.seed(seed)
  
  realp = p
  if ( p>10000 )
    p = 10000
  q = sample(4:8,1)
  n = ntrain
  pathway = matrix(FALSE,p,g)
  ngene = rpois(g,gsm)
  pathway[1:ngene[1],1] = TRUE
  for ( i in 2:g )
    pathway[sample(1:p,ngene[i]),i] = TRUE
    
  E0 = matrix(0,0,2)
    
  for ( i in 1:g )
  {
    Eg = matrix(0,0,2)
    mbr = which(pathway[,i])
    for ( j in mbr )
      for ( k in mbr )
      {
        if ( j>=k )
          next
        if ( k<=q )
        {
          if ( runif(1) < edii )
            Eg = rbind(Eg,matrix(c(j,k,k,j),2,2))
        }
        else if ( j<=q )
        {
          if ( runif(1) < ediu )
            Eg = rbind(Eg,matrix(c(j,k,k,j),2,2))
        }
        else
        {
          if ( runif(1) < eduu )
            Eg = rbind(Eg,matrix(c(j,k,k,j),2,2))
        }
      }
    E0 = rbind(E0,Eg)
  }
  E0 = unique(E0[order(E0[,1],E0[,2]),])
  nE0 = nrow(E0)

  idx = E0[,1]+p*(E0[,2]-1)

  A = matrix(0,p,p)
  A[idx] = sample(c(-1,1),nE0,TRUE)
  A = A*t(A)
  A[1:q,1:q] = abs(A[1:q,1:q])
 
  sumA = apply(abs(A),1,sum)
  maxsum = pmax(rep(sumA,p),rep(sumA,each=p))*1.1+0.1
  A = -A/maxsum
  diag(A) = 1

  pcor = -A[idx]/sqrt(A[E0[,1]+p*(E0[,1]-1)]*A[E0[,2]+p*(E0[,2]-1)])
 
  A = chol(A)
  A = backsolve(A,diag(p))
  dA = apply(A^2,1,sum)
  A = A / sqrt(dA)

  beta = c(runif(q,1,3),rep(0,p-q))

  X = matrix(rnorm(n*p),n) %*% t(A)
  y = drop(X%*%beta + sqrt(sigma2)*rnorm(n))

  Xtune = matrix(rnorm(ntune*p),ntune) %*% t(A)
  ytune = drop(Xtune%*%beta + sqrt(sigma2)*rnorm(ntune))

  Xtest = matrix(rnorm(ntest*p),ntest) %*% t(A)
  ytest = drop(Xtest%*%beta + sqrt(sigma2)*rnorm(ntest))

  if ( realp > p )
  {
    beta = c(beta,rep(0,realp-p))
    X = cbind(X,matrix(rnorm(n*(realp-p)),n))
    Xtune = cbind(Xtune,matrix(rnorm(n*(realp-p)),n))
    Xtest = cbind(Xtest,matrix(rnorm(n*(realp-p)),n))
    p = realp
  }

  if ( Gmode==0 )
    E = E0
  else if ( Gmode==1 )
  {
    e = nE0/2
    E = matrix(0,nE0,2)
    E[1:e,1] = sample(1:p,e,T)
    E[e+1:e,1] = sample(1:p,e,T)
    E[1:e,2] = E[e+1:e,1]
    E[e+1:e,2] = E[1:e,1]
    neq = which(E[,1]!=E[,2])
    E = E[neq,]  
    E = unique(E[order(E[,1],E[,2]),])
  }
  else if ( Gmode==2 )
  {
    set = abs(pcor) > Gthres
    E = E0[set,]
  }
  nE = nrow(E)
    
  data = list(p=p,q=q,g=g,n=n,ntune=ntune,ntest=ntest,sigma2=sigma2,edii=edii,ediu=ediu,eduu=eduu,Gmode=Gmode,Gthres=Gthres)
  data = c(data,list(pathway=pathway,ngene=ngene,pcor=pcor,E0=E0,nE0=nE0,E=E,nE=nE))
  data = c(data,list(beta=beta,X=X,y=y,Xtune=Xtune,ytune=ytune,Xtest=Xtest,ytest=ytest))

  if ( !is.null(savefile) )
    save(data,file=savefile)

  data
}



