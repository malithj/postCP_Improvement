plot.postCP <- function(x, line.prob="cp",line.mean="seg",rm.outliers=0,main=character(),xlab=NA,ylab=NA,p.col="blue",pch=16,p.cex=NA,m.col="brown",m.lty=1,m.lwd=1,l.col=NA,l.lty=NA,l.lwd=NA, ...){
# x: results from function postCP
  if(missing(x)){
     stop("Need to use results of postCP first")
  }
  
  if (is.element("data",names(x))){
    data=as.double(x$data);
    mu=as.double(x$means);
    n=length(data);
  } else line.mean="post";


   if (!is.element(line.mean,c("seg","post")))  stop("\nSpecify means as seg (segment) or post (posterior) means"); 
   if (!is.element("post.cp",names(x))) {
    stop("\nRerun postCP with keep=TRUE");
  }
   if (line.mean=="seg"&!is.element("data",names(x))) stop("\n Segment means can only be calculated if data is entered");
  level.based=is.element("post.level",names(x))  
  if (!level.based){
    if (!is.element(line.prob,c("cp","state"))) {stop("\nSpecify posterior probabilities as cp (change-point) or state (state)");} else{
    if (!is.element(line.prob,c("cp","state","level"))) {stop("\nSpecify posterior probabilities as cp (change-point), state (state) or level (level)");} }
  } 
 
  if (is.na(xlab)) xlab="Index";
  if (is.na(ylab)) ylab="data";
  range.100=range(data);
  if (line.prob=="cp") pprob=x$post.cp else{
  if (line.prob=="state") pprob=x$post.state;
  if (line.prob=="level")  pprob=x$post.level;}
  if (is.element("data",names(x))){
    if (rm.outliers!=0){
      ylim=quantile(data,probs=c(rm.outliers/2,(1-rm.outliers/2)))
     } else  ylim=range(data);

    if (pch==16){
      if (is.na(p.cex)) p.cex=1/log(n+1,100);
    }
    plot(data,ylim=ylim,main=main,xlab=xlab,ylab=ylab,col=p.col,pch=pch,cex=p.cex);
    if (line.mean=="seg"){
      index=diff(c(0,x$seg,n));
      means=mu[rep(1:length(mu),index)]; 
    }else{ 
      index=diff(c(0,x$bestcp,n));     
      means=as.vector(x$post.state%*%mu); 
    }
    lines(means,col=m.col,lty=m.lty,lwd=m.lwd);
    k1=1;
    crj=pretty(pprob);
    range1=diff(range(pprob));
    pprob=(pprob-min(pprob))/range1*diff(ylim)+ylim[1];
    mid1=round(length(crj)/2);
  } else{
      ylim=range(pprob);
      k1=2;
   }
  line.vec <- function(y){
    if (is.na(y[1]))  out=1:ncol(pprob) else{ 
      if (length(y)==1) out=rep(y,ncol(pprob)) else{
        if (length(y)>=ncol(pprob)) out=y else stop("\nSpecify line vector of at least same length as number of displayed lines");
      }
    }
    return(out);
  }
  l.col=line.vec(l.col);
  l.lwd=line.vec(l.lwd);
  l.lty=line.vec(l.lty);

  if (is.na(l.lty[1])) l.lty=1:ncol(pprob);
  if (is.na(l.lwd[1])) l.lwd=1;  
  plab=paste("P(",line.prob,")",sep="");
  if (!is.element("data",names(x))) {      
     plot(pprob[,1],type="l",main=main,xlab=xlab,ylim=ylim,col=l.col[1],lty=l.lty[1],lwd=l.lwd,ylab=plab);
  } else {  
     axis(side=4,labels=format(crj),at=crj/range1*diff(ylim)+ylim[1]);
     mtext(plab,at=mean(crj[mid1:(mid1+1)])/range1*diff(ylim)+ylim[1],side=4);
   }
  for (i in k1:ncol(pprob)) lines(pprob[,i],col=l.col[i],lty=l.lty[i],lwd=l.lwd);
  if (!is.element("data",names(x))) xlims=c(1,nrow(pprob)) else xlims=c(1,length(data));
  out=list(x=xlims,y=ylim);
  return(out);

}

