
if ( !require(glmnet) )
{
  install.packages("glmnet")
  library(glmnet)
}

if ( file.exists("~/project/EMSHS/RG/Ridge.R") )
  source("~/project/EMSHS/RG/Ridge.R")


SimLasso <- function(r,s,head,offset=0)
{
  M = length(s)
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
    fit = glmnet(dataset$X,dataset$y,intercept=FALSE,lambda=s)
    beta = matrix(coef(fit)[-1,],dataset$p)
    FNtune[,i] = apply(beta[1:dataset$q,]==0,2,sum)
    FPtune[,i] = apply(beta[(dataset$q+1):dataset$p,]!=0,2,sum)
    yhattune = dataset$Xtune %*% beta
    MSPEtune[,i] = apply((yhattune-as.vector(dataset$ytune))^2,2,mean)
    j = which.min(MSPEtune[,i])
    FN[i] = FNtune[j,i]
    FP[i] = FPtune[j,i]
    yhat = dataset$Xtest %*% beta[,j]
    MSPE[i] = mean((yhat-as.vector(dataset$ytest))^2)
  }
  after = proc.time() - before
  list(r=r,s=s,pred=T,vs=T,FNtune=FNtune,FPtune=FPtune,MSPEtune=MSPEtune,FN=FN,FP=FP,MSPE=MSPE,time=after[1])
}



SimALasso <- function(r,s1,s2,head,offset=0)
{
  M = length(s1)
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
    fit = Ridge(dataset$X,dataset$y,s2)
    w = 1/abs(fit)
    fit = glmnet(dataset$X,dataset$y,intercept=FALSE,penalty.factor=w,lambda=s1)
    beta = matrix(coef(fit)[-1,],dataset$p)
    FNtune[,i] = apply(beta[1:dataset$q,]==0,2,sum)
    FPtune[,i] = apply(beta[(dataset$q+1):dataset$p,]!=0,2,sum)
    yhattune = dataset$Xtune %*% beta
    MSPEtune[,i] = apply((yhattune-as.vector(dataset$ytune))^2,2,mean)
    j = which.min(MSPEtune[,i])
    FN[i] = FNtune[j,i]
    FP[i] = FPtune[j,i]
    yhat = dataset$Xtest %*% beta[,j]
    MSPE[i] = mean((yhat-as.vector(dataset$ytest))^2)
  }
  after = proc.time() - before
  list(r=r,s1=s1,s2=s2,pred=T,vs=T,FNtune=FNtune,FPtune=FPtune,MSPEtune=MSPEtune,FN=FN,FP=FP,MSPE=MSPE,time=after[1])
}
