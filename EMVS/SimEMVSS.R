if ( file.exists("~/project/EMSHS/EMVS/EMVS_Rcode_original.R") )
  source("~/project/EMSHS/EMVS/EMVS_Rcode_original.R")
if ( file.exists("EMVS/EMVS_Rcode_original.R") )
  source("EMVS/EMVS_Rcode_original.R")
#library(EMVS)

SimEMVSS <- function(r,v0,v1,eta=0,datapath,batch=0)
{
  D1 = length(v0)

  FNrate = matrix(0,D1,r)
  FPrate = matrix(0,D1,r)
  MSTE = matrix(0,D1,r)
  MSPE = matrix(0,D1,r)
  time = matrix(0,D1,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    if ( eta>0 )
    {
      G = diag(data$p)
      G[data$E[,1]+(data$E[,2]-1)*data$p] = eta
    }

    for ( d1 in 1:D1 )
    {
      if ( eta==0 )
        time[d1,i] = system.time(fit <- EMVS(data$y,data$X,v0=v0[d1],v1=v1,type="betabinomial",sigma_init=1,epsilon=1e-5,a=1,b=1))[1]
      else
        time[d1,i] = system.time(fit <- EMVS(data$y,data$X,v0=v0[d1],v1=v1,type="MRF",mu=0,Sigma=G,sigma_init=1,epsilon=1e-5,v1_g=v1))[1]

      beta = as.vector(fit$betas*(fit$p>0.5))

      FNrate[d1,i] = mean(fit$p[1,1:data$q]<0.5)
      FPrate[d1,i] = mean(fit$p[1,(data$q+1):data$p]>0.5)
    
      yhattune = data$Xtune %*% beta
      MSTE[d1,i] = mean((data$ytune-yhattune)^2)
      yhattest = data$Xtest %*% beta
      MSPE[d1,i] = mean((data$ytest-yhattest)^2)
    }
  }
  
  list(r=r,batch=batch,v0=v0,v1=v1,eta=eta,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}

