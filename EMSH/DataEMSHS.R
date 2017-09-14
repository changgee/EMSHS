if ( file.exists("~/project/EMSHS/EMSH/EMSHS.R") )
  source("~/project/EMSHS/EMSH/EMSHS.R")
if ( file.exists("EMSH/EMSHS.R") )
  source("EMSH/EMSHS.R")

DataEMSHS <- function(y,X,mu,nu,w=NULL,E=NULL,a_omega=NULL,b_omega=NULL,fold,k)
{
  r = ncol(fold)
  K = length(k)
  n = nrow(X)
  p = ncol(X)
  
  if ( is.null(w) )
    w = rep(1,n)

  if ( is.null(a_omega) )
    a_omega = 0

  D1 = length(mu)
  D2 = length(nu)
  D3 = length(a_omega)

  SSPECV = array(0,c(D1,D2,D3,r))
  beta = array(0,c(p,K,D1,D2,D3,r))
  lambda = array(0,c(p,K,D1,D2,D3,r))
  L = array(0,c(K,D1,D2,D3,r))
  
  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      for ( d1 in 1:D1 )
        for ( d2 in 1:D2 )
          for ( d3 in 1:D3 )
          {
            if ( is.null(E) | a_omega[d3] == 0 )
              fit = EMSHS(y[-ik],X[-ik,],mu[d1],nu[d2],w=w[-ik],descent="diagonal") 
            else
              fit = EMSHS(y[-ik],X[-ik,],mu[d1],nu[d2],E,w=w[-ik],a_omega=a_omega[d3],b_omega=b_omega,descent="diagonal")
            beta[,j,d1,d2,d3,i] = fit$beta
            lambda[,j,d1,d2,d3,i] = fit$lambda
            yhat = X[ik,] %*% beta[,j,d1,d2,d3,i]
            SSPECV[d1,d2,d3,i] = SSPECV[d1,d2,d3,i] + sum((y[ik]-yhat)^2)
          }
      
      L[j,,,,i] = apply(beta[,j,,,,i]!=0,2:4,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,n=n,K=K,fold=fold,k=k,mu=mu,nu=nu,a_omega=a_omega,b_omega=b_omega,beta=beta,L=L,lambda=lambda,SSPECV=SSPECV,time=after)
}


