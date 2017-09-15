if ( file.exists("~/project/EMSHS/EMVS/EMVS_Rcode_original.R") )
  source("~/project/EMSHS/EMVS/EMVS_Rcode_original.R")
if ( file.exists("EMVS/EMVS_Rcode_original.R") )
  source("EMVS/EMVS_Rcode_original.R")

DataEMVSS <- function(y,X,v0,v1,eta=0,E=NULL,fold,k)
{
  r = ncol(fold)
  K = length(k)
  n = nrow(X)
  p = ncol(X)
  
  M = length(v0)
  SSPECV = matrix(0,M,r)
  beta = array(0,c(p,K,M,r))
  L = array(0,c(K,M,r))  

  before = proc.time()
  
  if ( !is.null(E) )
  {
    G = matrix(0,p,p)
    G[E[,1]+(E[,2]-1)*p] = eta
  }
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])

      if ( !is.null(E) )
        fit <- EMVS(y[-ik],X[-ik,],v0=v0,v1=v1,type="MRF",mu=0,Sigma=G,sigma_init=1,epsilon=1e-4,v1_g=v1)
      else
        fit <- EMVS(y[-ik],X[-ik,],v0=v0,v1=v1,type="betabinomial",sigma_init=1,epsilon=1e-4,a=1,b=1)
    
      beta[,j,,i] = t(fit$betas*(fit$p>0.5))
      L[j,,i] = apply(fit$p>0.5,1,sum)
      yhat = X[ik,] %*% beta[,j,,i]
      SSPECV[,i] = SSPECV[,i] + apply((y[ik]-yhat)^2,2,sum)
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,n=n,K=K,fold=fold,k=k,v0=v0,v1=v1,eta=eta,beta=beta,L=L,SSPECV=SSPECV,time=after)
}

