if ( file.exists("~/project/EMSHS/EMSH/EMSHS.R") )
  source("~/project/EMSHS/EMSH/EMSHS.R")
if ( file.exists("EMSH/EMSHS.R") )
  source("EMSH/EMSHS.R")


SimEMSHS <- function(r,munu,nu,a_omega,datapath,batch=0)
{
  D1 = length(munu)
  D2 = length(nu)
  D3 = length(a_omega)

  FNrate = array(0,c(D1,D2,D3,r))
  FPrate = array(0,c(D1,D2,D3,r))
  MSTE = array(0,c(D1,D2,D3,r))
  MSPE = array(0,c(D1,D2,D3,r))
  omegaii = array(0,c(D1,D2,D3,r))
  omegaiu = array(0,c(D1,D2,D3,r))
  omegauu = array(0,c(D1,D2,D3,r))
  time = array(0,c(D1,D2,D3,r))

  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    i1 = data$E[,1]<=data$q
    i2 = data$E[,2]<=data$q
    ii = i1 & i2
    iu = xor(i1,i2)
    uu = !i1 & !i2
    
    for ( d1 in 1:D1 )
      for ( d2 in 1:D2 )
        for ( d3 in 1:D3 )
        {
          if ( a_omega[d3]>0 )
            time[d1,d2,d3,i] = system.time(fit <- EMSHS(data$y,data$X,munu[d1]-nu[d2],nu[d2],data$E,a_omega=a_omega[d3]))[1]
          else
            time[d1,d2,d3,i] = system.time(fit <- EMSHS(data$y,data$X,munu[d1]-nu[d2],nu[d2]))[1]
          FNrate[d1,d2,d3,i] = mean(fit$beta[1:data$q,1]==0)
          FPrate[d1,d2,d3,i] = mean(fit$beta[(data$q+1):data$p,1]!=0)
          MSTE[d1,d2,d3,i] = mean((data$ytune-data$Xtune%*%fit$beta[,1])^2)
          MSPE[d1,d2,d3,i] = mean((data$ytest-data$Xtest%*%fit$beta[,1])^2)
          if ( a_omega[d3]>0 )
          {
            omegaii[d1,d2,d3,i] = mean(fit$omega[ii,1])
            omegaiu[d1,d2,d3,i] = mean(fit$omega[iu,1])
            omegauu[d1,d2,d3,i] = mean(fit$omega[uu,1])
          }
        }
  }  
  
  list(r=r,batch=batch,munu=munu,nu=nu,a_omega=a_omega,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,omegaii=omegaii,omegaiu=omegaiu,omegauu=omegauu,time=time)
}

