plot.postCP.dev<-function(x, main=character(),xlab="index",ylab="probability",p.col="blue",pch=16,p.cex=NA,m.col="brown",m.lty=1,m.lwd=1,l.col=NA,l.lty=NA,l.lwd=NA){
  # x: results from function postCP
  if(missing(x)){
    stop("Need to use results of postCP first")
  }
  n = x$n

  plot(c(1,n),c(0,max(x$post.cp)),t="n",xlab=xlab,ylab=ylab,col=p.col,pch=pch,cex=p.cex)
  for (k in 1:ncol(x$post.cp))
    points(x$post.cp[,k],t='l',col=k,lty=k)
  title(main=main)
}