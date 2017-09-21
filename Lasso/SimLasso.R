
if ( !require(glmnet) )
{
  install.packages("glmnet")
  library(glmnet)
}

if ( !require(grpregOverlap) )
{
  install.packages("grpregOverlap")
  library(grpregOverlap)
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
      FNrate[d1,i] = mean(beta[1:data$q]==0)
      FPrate[d1,i] = mean(beta[(data$q+1):data$p]!=0)
      yhattune = data$Xtune %*% beta
      MSTE[d1,i] = mean((yhattune-data$ytune)^2)
      yhattest = data$Xtest %*% beta
      MSPE[d1,i] = mean((yhattest-data$ytest)^2)
    }
  }
  list(r=r,batch=batch,s=s,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}


SimGLasso <- function(r,lam,datapath,batch=0)
{
  D1 = length(lam)

  FNrate = matrix(0,D1,r)
  FPrate = matrix(0,D1,r)
  MSTE = matrix(0,D1,r)
  MSPE = matrix(0,D1,r)
  time = matrix(0,D1,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))

    grp = list()
    for ( j in 1:data$g )
      eval(parse(text=paste("grp = c(grp,list(g",j,"=which(data$pathway[,j])))",sep="")))
    rem = which(apply(data$pathway,1,sum)==0)
    cnt = 0
    for ( j in rem )
    {
      cnt = cnt + 1
      eval(parse(text=paste("grp = c(grp,list(g",data$g+cnt,"=j))",sep="")))
    }

    for ( d1 in 1:D1 )
    {
      time[d1,i] = system.time(fit <- grpregOverlap(data$X,data$y,grp,penalty="grLasso",lambda=lam[d2]))[1]
      beta = coef(fit)[-1]
      FNrate[d1,i] = mean(beta[1:data$q]==0)
      FPrate[d1,i] = mean(beta[(data$q+1):data$p]!=0)
      yhattune = data$Xtune %*% beta
      MSTE[d1,i] = mean((yhattune-data$ytune)^2)
      yhattest = data$Xtest %*% beta
      MSPE[d1,i] = mean((yhattest-data$ytest)^2)
    }
  }
  list(r=r,batch=batch,s=s,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}


SimALasso <- function(r,s1,s2,datapath,batch=0)
{
  D1 = length(s1)
  D2 = length(s2)

  FNrate = array(0,c(D1,D2,r))
  FPrate = array(0,c(D1,D2,r))
  MSTE = array(0,c(D1,D2,r))
  MSPE = array(0,c(D1,D2,r))
  time = array(0,c(D1,D2,r))
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))

    for ( d2 in 1:D2 )
    {
      t2 = system.time(fit <- Ridge(data$X,data$y,s2[d2]))[1]
      w = 1/abs(fit)

      for ( d1 in 1:D1 )
      {
        time[d1,d2,i] = t2 + system.time(fit <- glmnet(data$X,data$y,intercept=FALSE,penalty.factor=w,lambda=s1[d1]))[1]
        beta = coef(fit)[-1,]
        FNrate[d1,d2,i] = mean(beta[1:data$q]==0)
        FPrate[d1,d2,i] = mean(beta[(data$q+1):data$p]!=0)
        yhattune = data$Xtune %*% beta
        MSTE[d1,d2,i] = mean((yhattune-data$ytune)^2)
        yhattest = data$Xtest %*% beta
        MSPE[d1,d2,i] = mean((yhattest-data$ytest)^2)
      }
    }
  }
  list(r=r,batch=batch,s1=s1,s2=s2,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}



