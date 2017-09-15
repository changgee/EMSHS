if ( file.exists("~/project/EMSHS/EMVS/EMVS_Rcode_original.R") )
  source("~/project/EMSHS/EMVS/EMVS_Rcode_original.R")
if ( file.exists("EMVS/EMVS_Rcode_original.R") )
  source("EMVS/EMVS_Rcode_original.R")
#library(EMVS)

SimEMVSS <- function(r,v0,v1,eta=0,datapath,batch=0)
{
  D1 = length(v0)
  D2 = length(eta)

  FNrate = array(0,c(D1,D2,r))
  FPrate = array(0,c(D1,D2,r))
  MSTE = array(0,c(D1,D2,r))
  MSPE = array(0,c(D1,D2,r))
  time = array(0,c(D1,D2,r))
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    if ( data$p <= 10000 )
    {
      G = diag(0,data$p)
      G[data$E[,1]+(data$E[,2]-1)*data$p] = 1
    }

    for ( d2 in 1:D2 )
    {
      if ( eta[d2]==0 )
        time[,d2,i] = system.time(fit <- EMVS(data$y,data$X,v0=v0,v1=v1,type="betabinomial",sigma_init=1,epsilon=1e-3,a=1,b=1))[1]
      else
        time[,d2,i] = system.time(fit <- EMVS(data$y,data$X,v0=v0,v1=v1,type="MRF",mu=0,Sigma=G*eta[d2],sigma_init=1,epsilon=1e-3,v1_g=v1))[1]

      beta = t(fit$betas*(fit$p>0.5))

      FNrate[,d2,i] = apply(fit$p[,1:data$q]<0.5,1,mean)
      FPrate[,d2,i] = apply(fit$p[,(data$q+1):data$p]>0.5,1,mean)
    
      yhattune = data$Xtune %*% beta
      MSTE[,d2,i] = apply((data$ytune-yhattune)^2,2,mean)
      yhattest = data$Xtest %*% beta
      MSPE[,d2,i] = apply((data$ytest-yhattest)^2,2,mean)
    }
  }
  
  list(r=r,batch=batch,v0=v0,v1=v1,eta=eta,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}

SimEMVSS_old <- function(r,v0,v1,eta=0,datapath,batch=0)
{
  D1 = length(v0)
  D2 = length(eta)

  FNrate = array(0,c(D1,D2,r))
  FPrate = array(0,c(D1,D2,r))
  MSTE = array(0,c(D1,D2,r))
  MSPE = array(0,c(D1,D2,r))
  time = array(0,c(D1,D2,r))
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    G = diag(0,data$p)

    for ( d2 in 1:D2 )
    {
      G[data$E[,1]+(data$E[,2]-1)*data$p] = eta[d2]
      for ( d1 in 1:D1 )
      {
        if ( eta[d2]==0 )
          time[d1,d2,i] = system.time(fit <- EMVS(data$y,data$X,v0=v0[d1],v1=v1,type="betabinomial",sigma_init=1,epsilon=1e-3,a=1,b=1))[1]
        else
          time[d1,d2,i] = system.time(fit <- EMVS(data$y,data$X,v0=v0[d1],v1=v1,type="MRF",mu=0,Sigma=G,sigma_init=1,epsilon=1e-3,v1_g=v1))[1]

        beta = as.vector(fit$betas*(fit$p>0.5))

        FNrate[d1,d2,i] = mean(fit$p[1,1:data$q]<0.5)
        FPrate[d1,d2,i] = mean(fit$p[1,(data$q+1):data$p]>0.5)
    
        yhattune = data$Xtune %*% beta
        MSTE[d1,d2,i] = mean((data$ytune-yhattune)^2)
        yhattest = data$Xtest %*% beta
        MSPE[d1,d2,i] = mean((data$ytest-yhattest)^2)
      }
    }
  }
  
  list(r=r,batch=batch,v0=v0,v1=v1,eta=eta,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}

