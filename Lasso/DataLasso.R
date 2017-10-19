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

DataLasso <- function(y,X,w,s,fold,k)
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
      
      fit = glmnet(X[-ik,]*w[-ik],y[-ik]*w[-ik],intercept=FALSE,lambda=s)

      beta[,j,,i] = matrix(coef(fit)[-1,],p)
      L[j,,i] = apply(beta[,j,,i]!=0,2,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,N=N,K=K,fold=fold,k=k,s=s,beta=beta,L=L,SSPECV=SSPECV,time=after)
}


DataALasso <- function(y,X,w,s1,s2,fold,k)
{
  r = ncol(fold)
  K = length(k)
  N = nrow(X)
  p = ncol(X)
  
  D1 = length(s1)
  D2 = length(s2)
  SSPECV = array(0,c(D1,D2,r))
  beta = array(0,c(p,K,D1,D2,r))
  L = array(0,c(K,D1,D2,r))
  
  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])

      for ( d2 in 1:D2 )
      {
        fit = Ridge(X[-ik,]*w[-ik],y[-ik]*w[-ik],s2[d2])
        w = 1/abs(fit)
        fit = glmnet(X[-ik,]*w[-ik],y[-ik]*w[-ik],intercept=FALSE,penalty.factor=w,lambda=s1)
  
        beta[,j,,d2,i] = matrix(coef(fit)[-1,],p)
        L[j,,d2,i] = apply(beta[,j,,d2,i]!=0,2,sum)
        yhat = X[ik,] %*% beta[,j,,d2,i]
        SSPECV[,d2,i] = SSPECV[,d2,i] + apply((y[ik]-yhat)^2,2,sum)
      }
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,N=N,K=K,fold=fold,k=k,s1=s1,s2=s2,beta=beta,L=L,SSPECV=SSPECV,time=after)
}

DataGLasso <- function(y,X,w,pathway,lam,fold,k)
{
  r = ncol(fold)
  K = length(k)
  N = nrow(X)
  p = ncol(X)
  g = ncol(pathway)
  
  M = length(lam)
  SSPECV = matrix(0,M,r)
  beta = array(0,c(p,K,M,r))
  L = array(0,c(K,M,r))

  colnames(X) = 1:p
  grp = list()
  for ( j in 1:g )
    eval(parse(text=paste("grp = c(grp,list(gr",j,"=which(pathway[,j])))",sep="")))
  rem = which(apply(pathway,1,sum)==0)
  cnt = 0
  for ( j in rem )
  {
    cnt = cnt + 1
    eval(parse(text=paste("grp = c(grp,list(gr",g+cnt,"=j))",sep="")))
  }
  
  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      fit <- grpregOverlap(X[-ik,]*w[-ik],y[-ik]*w[-ik],grp,penalty="grLasso",lambda=lam)
      
      beta[,j,,i] = matrix(coef(fit)[-1,],p)
      L[j,,i] = apply(beta[,j,,i]!=0,2,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,N=N,K=K,fold=fold,k=k,lam=lam,beta=beta,L=L,SSPECV=SSPECV,time=after)
}