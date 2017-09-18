
library(R.matlab)

OpenMatlab <- function(port=9999,wait=10)
{
  if ( !exists("matlab") )
  {
    Matlab$startServer(port=port)
    matlab <<- Matlab(port=port,remote=TRUE)
    Sys.sleep(wait)
    open(matlab)
    setVerbose(matlab, threshold=-2)
    setOption(matlab, "readResult/interval", 10); # Default is 1 second
    setOption(matlab, "readResult/maxTries",30*(60/10)); # ~30 minutes

    evaluate(matlab,"cd Vannucci;")
  }
}

CloseMatlab <- function()
{
  if ( exists("matlab") )
    close(matlab)
}


VBSMRFPCA <- function(X,y,Xtest,ytest,Mpk,E,thres,niter,bin)
{
  setVariable(matlab,data=cbind(exp(y),1),X=X,Mpk=Mpk,E=E,niter=niter,bi=bin)
  
  evaluate(matlab,"[n p] = size(X);")
  evaluate(matlab,"[p K] = size(Mpk);")
  
  evaluate(matlab,"MRF=zeros(p);")
  evaluate(matlab,"MRF(E(:,1)+p*(E(:,2)-1)) = 1;")
  
  evaluate(matlab,"Pt = [(1:K)' sum(Mpk)'];")
  
  evaluate(matlab,"h = 0.1; Hg = h*eye(K);")
  evaluate(matlab,"h0 = 10^6;")
  evaluate(matlab,"phi = 0.05;")
  evaluate(matlab,"optPC = 1;")
  evaluate(matlab,"optPLS = 0;")
  evaluate(matlab,"mu = -2;")
  evaluate(matlab,"eta = 0.04;")
  evaluate(matlab,"a = 3; b= 0.5;")
  
  evaluate(matlab,"init = 5;")
  
  #  v = getVariable(matlab,c("X","data","Pt","K","Hg","h0","phi","init","niter","optPC","optPLS","MRF","eta","mu","Mpk","a","b"))
  evaluate(matlab,"[Theta,Gamma,logProb,numpath,numvar,propmove,move,Ylatent] = bvspathMRF_C(X,data,Pt,K,Hg,h0,phi,init,niter,optPC,optPLS,MRF,eta,mu,Mpk,a,b);")
  
  evaluate(matlab,"Theta = Theta((bi+1):niter);")
  evaluate(matlab,"Gamma = Gamma((bi+1):niter);")
  
  evaluate(matlab,"FreqTheta = varfreqMRF_C(Theta,K);")
  evaluate(matlab,"ThetaSel = find(FreqTheta>0.5);")
  
  evaluate(matlab,"pathOut = varfreqcondMioMRF_C(Gamma,Theta,p,ThetaSel,Mpk);")
  setVariable(matlab,ThresholdGamma=thres)
  evaluate(matlab,"GammaSel = [];")
  evaluate(matlab,"MSPE = [];")
  setVariable(matlab,Z=y,ZV=ytest,Xf=Xtest,flag="LS",optBeta=1)
  
  evaluate(matlab,"for j=1:length(ThresholdGamma)
                     Sel = pathOut.freq>ThresholdGamma(j);
                     GammaSel = [GammaSel,Sel];
                     if ( sum(Sel) == 0 )
                       MSE = mean((ZV-mean(Z)).^2);
                     else
                       [Yfpred,Yhat.R2,MSE,I] = predicopathMRF_Cnuovo(flag,X,Z,Xf,ZV,find(Sel),ThetaSel,Hg,[],[],optPC,optPLS,Mpk,optBeta,h);
                     end
                     MSPE = [MSPE MSE];
                   end")
  
  res <- getVariable(matlab,c("GammaSel","MSPE"))
  
  list(Gamma=res$GammaSel,MSPE=drop(res$MSPE))
}  



SimVBSMRF <- function(r,thres,niter,bin,datapath,batch=0)
{
  M = length(thres)

  FNrate = matrix(0,M,r)
  FPrate = matrix(0,M,r)
  MSTE = matrix(0,M,r)
  MSPE = matrix(0,M,r)
  time = rep(0,r)
  
  for ( i in 1:r )
  {
    print(i)
    load(sprintf("%s/data%03d",datapath,batch+i))
    
    before = as.vector(Sys.time())
    fit <- VBSMRFPCA(data$X,data$y,data$Xtune,data$ytune,data$pathway,data$E,thres,niter,bin)
    after = as.vector(Sys.time()) - before
    
    FNrate[,i] = apply(fit$Gamma[1:data$q,]==0,2,mean)
    FPrate[,i] = apply(fit$Gamma[(data$q+1):data$p,]==1,2,mean)
    MSTE[,i] = fit$MSPE
    
    fit <- VBSMRFPCA(data$X,data$y,data$Xtest,data$ytest,data$pathway,data$E,thres,niter,bin)
    MSPE[,i] = fit$MSPE
    time[i] = after[1]
  }  
  
  list(r=r,batch=batch,thres=thres,niter=niter,bin=bin,FNrate=FNrate,FPrate=FPrate,MSTE=MSTE,MSPE=MSPE,time=time)
}

