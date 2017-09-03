
if ( !require(glmnet) )
{
  install.packages("glmnet")
  library(glmnet)
}

if ( file.exists("~/project/EMSHS/Ridge/Ridge.R") )
  source("~/project/EMSHS/Ridge/Ridge.R")
if ( file.exists("~Ridge/Ridge.R") )
  source("Ridge/Ridge.R")


SimLasso <- function(r,s,datapath,batch=0)
{
  D1 = length(s)

  FNrate = matrix(0,D1,r)
  FPrate = matrix(0,D1,r)
  MSTE = matrix(0,D1,r)
  MSPE = matrix(0,D1,r)
  time = matrix(0,D1,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))

    for ( d1 in 1:D1 )
    {
      time[d1,i] = system.time(fit <- glmnet(data$X,data$y,intercept=FALSE,lambda=s[d1]))[1]
      beta = coef(fit)[-1,]
      FNrate[,i] = mean(beta[1:data$q]==0)
      FPrate[,i] = mean(beta[(data$q+1):data$p]!=0)
      yhattune = data$Xtune %*% beta
      MSTE[,i] = mean((yhattune-data$ytune)^2)
      yhattest = data$Xtest %*% beta
      MSPE[,i] = mean((yhattest-data$ytest)^2)
    }
  }
  list(r=r,batch=batch,s=s,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}


SimALasso <- function(r,s1,s2,datapath,batch=0)
{
  D1 = length(s1)

  FNrate = matrix(0,D1,r)
  FPrate = matrix(0,D1,r)
  MSTE = matrix(0,D1,r)
  MSPE = matrix(0,D1,r)
  time = matrix(0,D1,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))

    fit = Ridge(data$X,data$y,s2)
    w = 1/abs(fit)

    for ( d1 in 1:D1 )
    {
      time[d1,i] = system.time(fit <- glmnet(data$X,data$y,intercept=FALSE,penalty.factor=w,lambda=s1[d1]))[1]
      beta = coef(fit)[-1,]
      FNrate[,i] = mean(beta[1:data$q]==0)
      FPrate[,i] = mean(beta[(data$q+1):data$p]!=0)
      yhattune = data$Xtune %*% beta
      MSTE[,i] = mean((yhattune-data$ytune)^2)
      yhattest = data$Xtest %*% beta
      MSPE[,i] = mean((yhattest-data$ytest)^2)
    }
  }
  list(r=r,batch=batch,s1=s1,s2=s2,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}



