if ( !require(glmgraph) )
{
#  install.packages("glmgraph")
#  library(glmgraph)
}

if ( file.exists("~/project/EMSHS/Net/netreg1.2.R") )
  source("~/project/EMSHS/Net/netreg1.2.R")
if ( file.exists("Net/netreg1.2.R") )
  source("Net/netreg1.2.R")

# Pan et al (2010)
SimNet2 <- function(r,datapath,batch=0)
{
  FNrate = matrix(0,1000,r)
  FPrate = matrix(0,1000,r)
  MSTE = matrix(0,1000,r)
  MSPE = matrix(0,1000,r)
  nsol = rep(0,r)
  lam = matrix(0,1000,r)
  time = rep(0,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    
    time[i] = system.time(fit <- GBLasso(data$y,data$X,data$E,wt=rep(1,data$p)))[1]

    nsol[i] = fit$nsol
    lam[1:nsol[i],i] = fit$lambdas
    beta = t(fit$betas)
    FNrate[1:nsol[i],i] = apply(beta[1:data$q,]==0,2,mean)
    FPrate[1:nsol[i],i] = apply(beta[(data$q+1):data$p]!=0,2,mean)
    MSTE[1:nsol[i],i] = apply((data$ytune-data$Xtune%*%beta)^2,2,mean)
    MSPE[1:nsol[i],i] = apply((data$ytest-data$Xtest%*%beta)^2,2,mean)
  }  
  
  list(r=r,batch=batch,nsol=nsol,lam=lam,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}


# Li and Li (2008)
SimNet1 <- function(r,lam1,lam2,datapath,batch=0)
{
  D1 = length(lam1)
  D2 = length(lam2)

  FNrate = array(0,c(D1,D2,r))
  FPrate = array(0,c(D1,D2,r))
  MSTE = array(0,c(D1,D2,r))
  MSPE = array(0,c(D1,D2,r))
  time = array(0,c(D1,D2,r))
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    
    L = LapMat(data$E,data$p)

    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
      {
        time[d1,d2,i] = system.time(fit <- glmgraph(data$X,data$y,L,family="gaussian",penalty="lasso",lambda1=lam1[d1],lambda2=lam2[d2]))[1]

        beta = as.vector(fit$betas)
        FNrate[d1,d2,i] = mean(beta[1:data$q]==0)
        FPrate[d1,d2,i] = mean(beta[(data$q+1):data$p]!=0)
        MSTE[d1,d2,i] = mean((data$ytune-data$Xtune%*%beta)^2)
        MSPE[d1,d2,i] = mean((data$ytest-data$Xtest%*%beta)^2)
      }
  }  
  
  list(r=r,batch=batch,lam1=lam1,lam2=lam2,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
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
