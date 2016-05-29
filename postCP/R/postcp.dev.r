postcp <- function(formula, data=numeric(), bp=integer(), family = gaussian, homoscedastic=TRUE, par=numeric(), sigma=1.0) {
  X=model.matrix(formula,data)
  n=nrow(X)
  p=ncol(X)
  beta=matrix(par,nrow=p)
  Xbeta=X%*%beta
  K=ncol(beta)
  # compute the response variable y
  cp=bp
  CP=c(0,cp,n)
  y=rep(NA,n)
  for (k in 1:K) {
    idx=(CP[k]+1):(CP[k+1])
    y[idx]=rnorm(diff(range(idx))+1,Xbeta[idx,k],sigma)
  }
  # compute the log-evidence le
  # this is basically the only model-specific part !
  le=-(Xbeta-y)^2/2/sigma^2-log(sigma)-0.5*log(2*pi)
  # build workspace for forward-backward
  lFw=matrix(0,n,K)
  lFw[]=le[]
  lBk=matrix(0,n,K)
  # call forward-backward
  loglik=FwBk(lFw,lBk)-lchoose(n-1,K-1)
  # compute posterior distribution
  post=lFw+lBk
  post=exp(post-apply(post,1,logsumexp))
  # compute posterior cp distribution
  post.cp=matrix(0,n,K-1)
  for (k in 1:(K-1)) {
    aux=lFw[-n,k]+le[-1,k+1]+lBk[-1,k+1]
    post.cp[-n,k]=exp(aux-max(aux))
  }
  return(list(loglik=loglik,post=post,post.cp=post.cp))
}

logsumexp=function(l) {
  i=which.max(l);
  res=l[i]+log1p(sum(exp(l[-i]-l[i])));
  if (is.nan(res)) res=-Inf;
  return(res);
}
