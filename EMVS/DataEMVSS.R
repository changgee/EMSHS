  if ( file.exists("~/project/EMSHS/EMVS/EMVS_Rcode_original.R") )
  source("~/project/EMSHS/EMVS/EMVS_Rcode_original.R")
if ( file.exists("EMVS/EMVS_Rcode_original.R") )
  source("EMVS/EMVS_Rcode_original.R")

DataEMVSS <- function(y,X,w,v0,v1,eta=0,E=NULL,fold,k)
{
  r = ncol(fold)
  K = length(k)
  n = nrow(X)
  p = ncol(X)
  
  D1 = length(v0)
  D2 = length(eta)
  SSPECV = array(0,c(D1,D2,r))
  beta = array(0,c(p,K,D1,D2,r))
  L = array(0,c(K,D1,D2,r))  

  before = proc.time()
  
  if ( !is.null(E))
  {
    G = matrix(0,p,p)
    G[E[,1]+(E[,2]-1)*p] = 1
  }
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      for ( d2 in 1:D2 )
      {
        if ( !is.null(E) & eta[d2] > 0 )
          fit <- EMVS(y[-ik]*w[-ik],X[-ik,]*w[-ik],v0=v0,v1=v1,type="MRF",mu=0,Sigma=G*eta[d2],sigma_init=1,epsilon=1e-4,v1_g=v1)
        else
          fit <- EMVS(y[-ik]*w[-ik],X[-ik,]*w[-ik],v0=v0,v1=v1,type="betabinomial",sigma_init=1,epsilon=1e-4,a=1,b=1)
      
        beta[,j,,d2,i] = t(fit$betas*(fit$p>0.5))
        L[j,,d2,i] = apply(fit$p>0.5,1,sum)
        yhat = X[ik,] %*% beta[,j,,d2,i]
        SSPECV[,d2,i] = SSPECV[,d2,i] + apply((y[ik]-yhat)^2,2,sum)
      }
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,n=n,K=K,fold=fold,k=k,v0=v0,v1=v1,eta=eta,beta=beta,L=L,SSPECV=SSPECV,time=after)
}

