if ( !require(glmgraph) )
{
  install.packages("glmgraph")
  library(glmgraph)
}

if ( file.exists("~/project/EMSHS/Net/netreg1.2.R") )
  source("~/project/EMSHS/Net/netreg1.2.R")
if ( file.exists("Net/netreg1.2.R") )
  source("Net/netreg1.2.R")

DataNet1 <- function(y,X,E,lam1,fold,k)
{
  r = ncol(fold)
  K = length(k)
  n = nrow(X)
  p = ncol(X)
  
  D1 = length(lam1)
  D2 = 9
  SSPECV = array(0,c(D1,D2,r))
  beta = array(0,c(p+1,K,D1,D2,r))
  L = array(0,c(K,D1,D2,r))
  
  LM = LapMat(E,p)

  before = proc.time()
  
  for ( i in 1:r )
  {
    for ( j in 1:K )
    {
      ik = which(fold[,i]==k[j])
      
      fit = glmgraph(X[-ik,],y[-ik],LM,family="gaussian",penalty="lasso",lambda1=lam1)

      for ( d2 in 1:D2 )
      {
        beta[-1,j,,d2,i] = fit$betas[[d2]]
      
        L[j,,d2,i] = apply(beta[-1,j,,d2,i]!=0,2,sum)
        yhat = cbind(1,X[ik,]) %*% beta[,j,,d2,i]
        SSPECV[,d2,i] = SSPECV[,d2,i] + apply((y[ik]-yhat)^2,2,sum)
      }
    }
  }
  
  after = proc.time()-before
  
  list(r=r,p=p,n=n,K=K,fold=fold,k=k,lam1=lam1,beta=beta,L=L,SSPECV=SSPECV,time=after)
}


LapMat <- function(E,p)
{
  L = matrix(0,p,p)
  e = nrow(E)
  for ( i in 1:e )
  {
    L[E[i,1],E[i,1]] = L[E[i,1],E[i,1]] + 1
    L[E[i,1],E[i,2]] = -1
  }
  for ( i in 1:p )
  {
    if ( L[i,i] == 0 )
      next
    d = sqrt(L[i,i])
    L[,i] = L[,i]/d
    L[i,] = L[i,]/d
  }
  L
}

