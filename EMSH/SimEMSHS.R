<<<<<<< HEAD
source("EMSH/EMSHS.R")
=======
source("/longgroup/changgee/Prj2/EMSH/EMSHS.R")
>>>>>>> 7d4f4242f163c24fce8446d9b871c4c6c5607a1f

SimEMShrink <- function(r,mu,nu,useE=TRUE)
{
  l = length(mu)
  FN = rep(0,l)
  FP = rep(0,l)
  MSPE = rep(0,l)
  
  before = proc.time()
  for ( i in 1:r )
  {
    source("SimulSet.R")
    for ( s in 1:l )
    {
      if ( useE )
        res = EMShrink(y,X,mu[s],nu,t(E))
      else
        res = EMShrink(y,X,mu[s],nu)      
      
      FN[s] = FN[s] + mean(res$beta[1:10] == 0)
      FP[s] = FP[s] + mean(res$beta[11:100] != 0)
      
      yhat = Xtest %*% res$beta
      MSPE[s] = MSPE[s] + mean((ytest-yhat)^2)
    }
  }
  after = proc.time() - before
  list(mu=mu,nu=nu,FN=FN/r,FP=FP/r,MSPE=MSPE/r,r=r,timePerDataset=after[1]/r) 
}

SimEMSHS <- function(r,mu,nu,useE=TRUE,a_omega=2,head,offset=0)
{
  M = length(mu)
  FNtune = matrix(0,M,r)
  FPtune = matrix(0,M,r)
  MSPEtune = matrix(0,M,r)
  FN = rep(0,r)
  FP = rep(0,r)
  MSPE = rep(0,r)
  
  before = proc.time()
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/dataset%03d",head,offset+i))
    if ( i==1 )
      beta = matrix(0,dataset$p,M)
    
    for ( j in 1:M )
    {
      print(j)
      if ( useE )
        fit = EMSHS(dataset$y,dataset$X,mu[j],nu,dataset$E,a_omega=a_omega,descent="diagonal")
      else
        fit = EMSHS(dataset$y,dataset$X,mu[j],nu,descent="diagonal") 
      beta[,j] = fit$beta
      FNtune[j,i] = sum(fit$beta[1:dataset$q,]==0)
      FPtune[j,i] = sum(fit$beta[(dataset$q+1):dataset$p,]!=0)
    }
    
    yhattune = dataset$Xtune %*% beta
    MSPEtune[,i] = apply((dataset$ytune-yhattune)^2,2,mean)
    j = which.min(MSPEtune[,i])
    FN[i] = FNtune[j,i]
    FP[i] = FPtune[j,i]
    yhat = dataset$Xtest %*% beta[,j]
    MSPE[i] = mean((dataset$ytest-yhat)^2)
  }  
  
  after = proc.time() - before
  list(r=r,mu=mu,nu=nu,c=c,pred=T,vs=T,FNtune=FNtune,FPtune=FPtune,MSPEtune=MSPEtune,FN=FN,FP=FP,MSPE=MSPE,time=after[1])
}

