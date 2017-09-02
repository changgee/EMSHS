
Ridge <- function(X,y,lambda)
{
  n = length(y)
  Xy = drop(t(X)%*%y)
  beta = c()
  for ( l in lambda )
  {
    tmp1 = Xy / l
    tmp2 = X %*% tmp1
    tmp3 = diag(n) + X %*% t(X) / l
    ctmp3 = chol(tmp3)
    Xbeta = backsolve(ctmp3,forwardsolve(t(ctmp3),tmp2))
    tmp2 = (t(X) %*% Xbeta) / l
    beta = cbind(beta,tmp1-tmp2,deparse.level=0)
  }
  beta
}


