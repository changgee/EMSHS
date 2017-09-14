if ( file.exists("~/project/EMSHS/EMSH/EMSHS.R") )
  source("~/project/EMSHS/EMSH/EMSHS.R")
if ( file.exists("EMSH/EMSHS.R") )
  source("EMSH/EMSHS.R")

DataEMSHS_CV <- function(y,X,mu,nu,w=NULL,E=NULL,a_omega=NULL,b_omega=NULL,fold,k)
{
  r = ncol(fold)
  K = length(k)
  n = nrow(X)
  p = ncol(X)
  
  M = length(mu)
  SSPECV = matrix(0,M,r)
  beta = array(0,c(p,K,M,r))
  lambda = array(0,c(p,K,M,r))
  L = array(0,c(K,M,r))
  
  if ( is.null(w) )
    w = rep(1,n)

  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      for ( m in 1:M )
      {
        if ( is.null(E) )
          fit = EMSHS(y[-ik],X[-ik,],mu[m],nu,w=w[-ik],descent="diagonal") 
        else
          fit = EMSHS(y[-ik],X[-ik,],mu[m],nu,E,w=w[-ik],a_omega=a_omega,b_omega=b_omega,descent="diagonal")
        beta[,j,m,i] = fit$beta
        lambda[,j,m,i] = fit$lambda
      }
      
      L[j,,i] = apply(beta[,j,,i]!=0,2,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,n=n,K=K,fold=fold,k=k,mu=mu,nu=nu,a_omega=a_omega,b_omega=b_omega,beta=beta,L=L,lambda=lambda,SSPECV=SSPECV,time=after)
}


