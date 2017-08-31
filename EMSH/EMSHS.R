if ( file.exists("/longgroup/changgee/Prj2/EMSH/DWL.R") )
  source("/longgroup/changgee/Prj2/EMSH/DWL.R")
if ( file.exists("EMSH/DWL.R") )
  source("EMSH/DWL.R")

# Function EMSHS
# y: n by 1 response
# X: n by p predictors
# mus: vector of mu
# nu: nu
# E: e by 2 matrix with edges, must be sorted by E[,1] and then E[,2]
#    a single edge (j,k) must be duplicated to (k,j) in E
#    NULL if no edge

EMSHS <- function(y,X,mus,nu,E=NULL,a_sigma=1,b_sigma=1,a_omega=2,b_omega=1,sigma_init=NULL,alpha_init=NULL,descent="diagonal",eps=1e-5)
{
  n = nrow(X)
  p = ncol(X)
  Xy = drop(t(X) %*% y)
  
  # Initialize graph information
  # nadj[i] has the number of neighboring variables of variable i
  # idx[i] keeps the beginning and ending index in E[,2] for neighboring variables of variable i
  nadj = rep(0,p)
  if ( is.null(E) )
  {
    e = 0
    descent = "diagonal"
  }
  else
  {
    e = nrow(E)
    for ( i in 1:e )
      nadj[E[i,1]] = nadj[E[i,1]] + 1
  }
  idx = diffinv(nadj)  

  
  # initialize for DWL
  DWL.n <<- n
  DWL.p <<- p
  
  DWL.X <<- as.matrix(X)
  DWL.y <<- as.vector(y)

  DWL.Init()  

  
  
  # Initialize beta
  beta = rep(0,p)

    
  
  # Initialize sigma
  c1 = DWL.yy/2 + b_sigma
  c2 = 0
  c3 = n+p+2*a_sigma+2
  if ( is.null(sigma_init) )
    sigma = (c2+sqrt(c2^2+8*c1*c3))/(2*c3)
  else
    sigma = sigma_init
  
  
  # Initialize storages
  M = length(mus)
  
  niters = rep(0,M)
  Qs = rep(0,M)
  betas = matrix(0,p,M)
  sigmas = rep(0,M)
  alphas = matrix(0,p,M)
  

  
  for ( mm in 1:M )
  {
    mu = mus[mm]
    
    # Initialize alpha
    if ( is.null(alpha_init) )
      alpha = rep(mus[mm],p)
    else
      alpha = alpha_init

    niter = 0

    repeat
    {
      niter = niter + 1
      
      ## E-Steps
      
      # E-Step: omega
      Eomega = 2*nu*a_omega / (2*nu*b_omega + (alpha[E[,1]]-alpha[E[,2]])^2)
      
  
      pQ = c3*log(sigma) + c1/sigma^2 + c2/sigma - sum(alpha) + sum((alpha-mu)^2)/2/nu + sum(Eomega*(alpha[E[,1]]-alpha[E[,2]])^2)/4/nu
      
      
      ## M-Step
      
      # M-Step: beta 
      res = DWL(sigma*exp(alpha))
      beta = res$coef
      c1 = (DWL.yy - sum((DWL.Xy+DWL.C)*beta))/2 + b_sigma
      c2 = sum(exp(alpha)*abs(beta))

      
      # M-Step: sigma
      sigma = (c2+sqrt(c2^2+8*c1*c3))/(2*c3)


  
      # M-Step: alpha
      if ( descent == "newton" )
      {
        EOmega = matrix(0,p,p)
        EOmega[E[,1]+(E[,2]-1)*p] = -Eomega
        diag(EOmega) = 1-apply(EOmega,1,sum)
        
        H = sigma*EOmega + nu*diag(exp(alpha)*abs(beta),p)
        tmp = drop(EOmega%*%(alpha-mu))
        g = sigma*tmp - nu*sigma*rep(1,p) + nu*exp(alpha)*abs(beta)
        f = sigma*sum(tmp*(alpha-mu))/2 - nu*sigma*sum(alpha) + nu*sum(exp(alpha)*abs(beta)) 
        
        cH = chol(H)
        dir = backsolve(cH,forwardsolve(t(cH),g))
        maxdir = max(abs(dir))
        m = sum(g*dir)
        
        ss = 1
        repeat
        {
          nalpha = alpha - ss*dir
          nc2 = sum(exp(nalpha)*abs(beta))
          tmp = drop(EOmega%*%(nalpha-mu))
          nf = sigma*sum(tmp*(nalpha-mu))/2 - nu*sigma*sum(nalpha) + nu*nc2 
          
          if ( ss*maxdir < eps | f-nf > ss*m*0.05 )
          {
            alpha = nalpha
            c2 = nc2
            break
          }
          ss = ss/2
        }      
      }
      else if ( descent == "diagonal" )
      {
        h = rep(0,p)
        g = rep(0,p)
        
        for ( i in 1:p )
        {
          if ( idx[i] < idx[i+1] )
          {
            adjidx = (idx[i]+1):idx[i+1]
            adjvar = E[adjidx,2]
            sEomega = sum(Eomega[adjidx])
            sEomegaalpha = sum(Eomega[adjidx]*alpha[adjvar])
          }
          else
          {
            sEomega = 0
            sEomegaalpha = 0
          }
          
          h[i] = sigma*(1+sEomega) + nu*exp(alpha[i])*abs(beta[i])
          g[i] = sigma*((1+sEomega)*alpha[i] - mu - sEomegaalpha) - nu*sigma + nu*exp(alpha[i])*abs(beta[i])
        }
        f = -sigma*(mu+nu)*sum(alpha) + nu*c2 + sigma*sum(alpha^2)/2 + sigma*sum(Eomega*(alpha[E[,1]]-alpha[E[,2]])^2)/4
        
        dir = g/h
        maxdir = max(abs(dir))
        m = sum(g*dir)
        
        ss = 1
        repeat
        {
          nalpha = alpha - ss*dir
          nc2 = sum(exp(nalpha)*abs(beta))
          nf = -sigma*(mu+nu)*sum(nalpha) + nu*nc2 + sigma*sum(nalpha^2)/2 + sigma*sum(Eomega*(nalpha[E[,1]]-nalpha[E[,2]])^2)/4
          
          if ( ss*maxdir < eps | f-nf > ss*m*0.05 )
          {
            alpha = nalpha
            c2 = nc2
            break
          }
          ss = ss/2
        }      
      }


      Q = c3*log(sigma) + c1/sigma^2 + c2/sigma - sum(alpha) + sum((alpha-mu)^2)/2/nu + sum(Eomega*(alpha[E[,1]]-alpha[E[,2]])^2)/4/nu
      
      if ( (pQ-Q)/abs(Q) < eps )
        break
      
    }
  
    niters[mm] = niter
    Qs[mm] = Q
    betas[,mm] = beta
    sigmas[mm] = sigma
    alphas[,mm] = alpha
  
  }
  
  list(niter=niters,Q=Qs,beta=betas,sigma=sigmas,alpha=alphas,lambda=exp(alphas),y=y,X=X,E=E,mu=mus,nu=nu,a_sigma=a_sigma,b_sigma=b_sigma,a_omega=a_omega,b_omega=b_omega,descent=descent,eps=eps)
}




