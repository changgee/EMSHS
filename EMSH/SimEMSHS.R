if ( file.exists("~/project/EMSHS/EMSH/EMSHS.R") )
  source("~/project/EMSHS/EMSH/EMSHS.R")
if ( file.exists("EMSH/EMSHS.R") )
  source("EMSH/EMSHS.R")


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

