if ( !require(glmnet) )
{
  install.packages("glmnet")
  library(glmnet)
}

if ( !require(grpregOverlap) )
{
  install.packages("grpregOverlap")
  library(grpregOverlap)
}

source("Ridge/Ridge.R")

DataLasso <- function(y,X,s,fold,k)
{
  r = ncol(fold)
  K = length(k)
  N = nrow(X)
  p = ncol(X)
  
  M = length(s)
  SSPECV = matrix(0,M,r)
  beta = array(0,c(p,K,M,r))
  L = array(0,c(K,M,r))
  
  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      fit = glmnet(X[-ik,],y[-ik],intercept=FALSE,lambda=s)

      beta[,j,,i] = matrix(coef(fit)[-1,],p)
      L[j,,i] = apply(beta[,j,,i]!=0,2,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,N=N,K=K,fold=fold,k=k,s=s,beta=beta,L=L,SSPECV=SSPECV,time=after)
}


DataALasso <- function(y,X,s1,s2,fold,k)
{
  r = ncol(fold)
  K = length(k)
  N = nrow(X)
  p = ncol(X)
  
  M = length(s1)
  SSPECV = matrix(0,M,r)
  beta = array(0,c(p,K,M,r))
  L = array(0,c(K,M,r))
  
  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      fit = Ridge(X[-ik,],y[-ik],s2)
      w = 1/abs(fit)
      fit = glmnet(X[-ik,],y[-ik],intercept=FALSE,penalty.factor=w,lambda=s1)

      beta[,j,,i] = matrix(coef(fit)[-1,],p)
      L[j,,i] = apply(beta[,j,,i]!=0,2,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,N=N,K=K,fold=fold,k=k,s1=s1,s2=s2,beta=beta,L=L,SSPECV=SSPECV,time=after)
}

DataGLasso <- function(y,X,grp,lam,fold,k)
{
  r = ncol(fold)
  K = length(k)
  N = nrow(X)
  p = ncol(X)
  
  M = length(lam)
  SSPECV = matrix(0,M,r)
  beta = array(0,c(p,K,M,r))
  L = array(0,c(K,M,r))
  
  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      fit <- grpregOverlap(data$X[-ik,],data$y[-ik],grp,penalty="grLasso",lambda=lam)
      
      beta[,j,,i] = matrix(coef(fit)[-1,],p)
      L[j,,i] = apply(beta[,j,,i]!=0,2,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,N=N,K=K,fold=fold,k=k,lam=lam,beta=beta,L=L,SSPECV=SSPECV,time=after)
}