if ( file.exists("~/project/EMSHS/EMVS/EMVS_Rcode_original.R") )
  source("~/project/EMSHS/EMVS/EMVS_Rcode_original.R")
if ( file.exists("EMVS/EMVS_Rcode_original.R") )
  source("EMVS/EMVS_Rcode_original.R")
#library(EMVS)

SimEMVS <- function(r,v0,v1,datapath,batch=0)
{
  M = length(v0)

  FNrate = matrix(0,M,r)
  FPrate = matrix(0,M,r)
  MSTE = matrix(0,M,r)
  MSPE = matrix(0,M,r)
  time = matrix(0,M,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))

    fit <- EMVS(data$y,data$X,v0=v0,v1=v1,type="betabinomial",sigma_init=1,epsilon=1e-5,a=1,b=1)
    beta = t(fit$betas * (fit$p>0.5))
    FNtune[,i] = apply(fit$p[,1:dataset$q]<0.5,1,sum)
    FPtune[,i] = apply(fit$p[,(dataset$q+1):dataset$p]>0.5,1,sum)
    
    yhattune = dataset$Xtune %*% beta
    MSPEtune[,i] = apply((dataset$ytune-yhattune)^2,2,mean)
    j = which.min(MSPEtune[,i])
    FN[i] = FNtune[j,i]
    FP[i] = FPtune[j,i]
    yhat = dataset$Xtest %*% beta[,j]
    MSPE[i] = mean((dataset$ytest-yhat)^2)
  }  
  
  after = proc.time() - before
  list(r=r,v0=v0,v1=v1,pred=T,vs=T,FNtune=FNtune,FPtune=FPtune,MSPEtune=MSPEtune,FN=FN,FP=FP,MSPE=MSPE,time=after[1])
}


SimEMVSS <- function(r,v0,v1,eta,datapath,batch=0)
{
  M = length(v0)
  FNtune = matrix(0,M,r)
  FPtune = matrix(0,M,r)
  MSPEtune = matrix(0,M,r)
  FN = rep(0,r)  
  FP = rep(0,r)
  MSPE = rep(0,r)
  
  before = proc.time()
  for ( i in 1:r )
  {
    load(sprintf("%s/dataset%03d",head,offset+i))
    if ( i==1 )
      G = diag(dataset$p)

    G[dataset$E[,1]+(dataset$E[,2]-1)*dataset$p] = eta
    
#    for ( j in 1:M )
#    {
#      fit <- EMVS(dataset$y,dataset$X,v0=v0[j],v1=v1,type="MRF",mu=0,Sigma=G,temperature=temp,sigma_init=1,epsilon=1e-5,v1_g=v1)
#      beta[,j] = fit$betas * (fit$p>0.5)
#      FNtune[j,i] = sum(fit$p[1:dataset$q]<0.5)
#      FPtune[j,i] = sum(fit$p[(dataset$q+1):dataset$p]>0.5)
#    }
    
    fit <- EMVS(dataset$y,dataset$X,v0=v0,v1=v1,type="MRF",mu=0,Sigma=G,temperature=temp,sigma_init=1,epsilon=1e-5,v1_g=v1)
    beta = t(fit$betas * (fit$p>0.5))
    FNtune[,i] = apply(fit$p[,1:dataset$q]<0.5,1,sum)
    FPtune[,i] = apply(fit$p[,(dataset$q+1):dataset$p]>0.5,1,sum)

    G[dataset$E[,1]+(dataset$E[,2]-1)*dataset$p] = 0
    
    yhattune = dataset$Xtune %*% beta
    MSPEtune[,i] = apply((dataset$ytune-yhattune)^2,2,mean)
    j = which.min(MSPEtune[,i])
    FN[i] = FNtune[j,i]
    FP[i] = FPtune[j,i]
    yhat = dataset$Xtest %*% beta[,j]
    MSPE[i] = mean((dataset$ytest-yhat)^2)
  }  
  
  after = proc.time() - before
  list(r=r,v0=v0,v1=v1,eta=eta,pred=T,vs=T,FNtune=FNtune,FPtune=FPtune,MSPEtune=MSPEtune,FN=FN,FP=FP,MSPE=MSPE,time=after[1])
}
