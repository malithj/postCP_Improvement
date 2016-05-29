postCPsample <- function(postCP.res, nsamples=100, gen.data="n", prior=0.5, prior.type="n", verbose=TRUE,debug=FALSE) UseMethod("postCPsample")

postCPsample.default <- function(postCP.res, nsamples=100, gen.data="n", prior=0.5, prior.type="n",verbose=TRUE,debug=FALSE){
# postCP.res: results from function postCP
# nsamples: number of sets of generated changepoints
# gen.data: generate a matrix of data, one row of length n for each of nsamples replicates
  if(missing(postCP.res)){
     stop("Need to use results of postCP first")
  }
  data=as.double(postCP.res$data)
  mu=as.double(postCP.res$means)
  sigma=as.double(postCP.res$sds)
  prior=as.double(postCP.res$prior)
  ptype=as.integer(postCP.res$prior.type)
  lprob=numeric()
  if (!is.element("lbackward",names(postCP.res))) {
    stop("\nRe-run postCP with keep=TRUE");
  }
  level.based=is.element("post.level",names(postCP.res)) 
  if (!level.based) lbackward=as.double(as.vector(postCP.res$lbackward)) else {
     warning("It is advisable to run postCPsample directly within postCP due to roundoff error");
  lbackward=as.double(unlist(postCP.res$lbackward))}
  
  if (!is.element("lprob",names(postCP.res))){
  n=postCP.res$n
  J=length(mu)
  if (!is.element("nseg",names(postCP.res))) {
    stop("\nRe-run postCP with keep=TRUE");
  } else   nseg=postCP.res$nseg
  model=1
  if (postCP.res$model=="normal") model=2} else{
    n=nrow(postCP.res$lprob);
    J=ncol(postCP.res$lprob);
    lprob=as.vector(postCP.res$lprob);
    if (gen.data=="p") print("Please choose non-parametric model");
  }

 if (!is.element("nseg",names(postCP.res))) nseg=postCP.res$nseg;
 probs=(length(lprob)>0)
 if (prior.type=="n") ptype=1 else
 if (prior.type=="o") ptype=2 else
 if (prior.type=="s") ptype=3 
  cpsize=nsamples*(nseg-1)

  cpvector=rep(0,cpsize)
  levvector=rep(0,cpsize)

 # if (!is.element("lbackward",names(postCP.res))) {
  .C("postCPsample",rdata=data,rmu=mu,rsigma=sigma,lprob=as.double(lprob),rprobs=as.integer(probs),rlbackward=as.double(lbackward),cpvector=as.double(cpvector),levvector=as.double(levvector),rnsamples=as.integer(nsamples),nn=as.integer(n),JJ=as.integer(J),rnseg=as.integer(nseg),rmodel=as.integer(model),rprior=as.double(prior),rpriortype=as.integer(ptype),rlevelbased=as.integer(level.based),rverbose=as.integer(verbose),rlevelbased=as.integer(level.based),rdebug=as.integer(debug), DUP=FALSE,PACKAGE="postCP")
  cpmatrix=matrix(cpvector,ncol=nseg-1)
  levmatrix=matrix(levvector,ncol=nseg-1)
#} else{
#  .C("postCPsampleext",lprob=as.double(lprob),rlbackward=as.double(lbackward),as.double(cpvector),rnsamples=as.integer(nsamples),nn=as.integer(n),JJ=as.integer(J),rprior=as.double(prior),rpriortype=as.integer(ptype),rverbose=as.integer(verbose),rdebug=as.integer(debug), DUP=FALSE,PACKAGE="postCP")
 # cpmatrix=matrix(cpvector,ncol=J-1)}
  rmu=postCP.res$means

  if (level.based) level.ind=postCP.res$level.ind

  if (!level.based) rmu1=rmu else  rmu1=rmu[level.ind[diff(c(0,level.ind))!=0]]
  fbsample <- list(changepoints=cpmatrix)
  if (level.based) {
	levmatrix=matrix(levvector,ncol=nseg-1)
	cpmatrix=cbind(cpmatrix,levmatrix)
      }	
     if (gen.data=="p"|gen.data=="np"){
	getSegment <- function(cpv,n,level.based=FALSE) {
		if (!level.based) out=rep(1:(length(cpv)+1),diff(c(0,cpv,n))) else {  
			lv = c(1,cpv[(length(cpv)/2+1):length(cpv)])
			cpv=cpv[1:(length(cpv)/2)]
			out=rep(lv,diff(c(0,cpv,n)))
		}
  		return(out)}
		if (gen.data=="p"){
			seglabs <- apply(cpmatrix,1,getSegment,n=n,level.based=level.based) 
		if (model==1) x=matrix(rpois(length(seglabs),rmu1[seglabs]),ncol=n,byrow=T)
		if (model==2) x=matrix(rnorm(length(seglabs),rmu1[seglabs],sigma),ncol=n,byrow=T)		
		} else{
			getIndex <- function(i,cpv,n,level.based) which(getSegment(cpv,n,level.based)==i)
			getSample <- function (cpv,n,level.based) {
				if (!level.based) levels=1:(length(cpv)+1) else levels=unique(c(1,cpv[(length(cpv)/2+1):length(cpv)]))
				ind.list=lapply(levels,getIndex,cpv,n=n,level.based=level.based) 
				ind.order=order(unlist(ind.list))
				getSample2 <- function(ind) {
					if (length(ind)>1) out=sample(ind,length(ind),replace=TRUE) else out=ind;
					return(out);
				}
				unlist(lapply(ind.list,getSample2))[ind.order]
			}	
			seglabs=matrix(apply(cpmatrix,1,getSample,n=n,level.based=level.based),ncol=nsamples)
			x=matrix(data[seglabs],ncol=n,byrow=T)
		}
	fbsample$x=x
	names(fbsample)[2] <- "data"
	}
     if (!is.element(gen.data,c("n","np","p"))) print("Please choose np (non-parametric resampling), p (parametric resampling) or n (no resampling)");

  class(fbsample) <- "postCPsample"
  fbsample
}

print.postCPsample <- function(x,...)
{ cat("\nGenerated change-points ($changepoints):\n")
  str(x[[1]])
  if (length(x)==2){
    cat("\nMatrix of generated data ($data):\n")
    str(x[[2]])
  }
}

postCP <- function(data=numeric(),seg=integer(),model=1,lprob=numeric(),level.ind=numeric(),keep=TRUE,ci=0.9,viterbi=TRUE,initsegci=TRUE,nsamples=0,gen.data="n",prior=0.5,prior.type="n",epsilon=1e-9,disp.equal=TRUE,eps.nb=1e-8,verbose=TRUE,debug=FALSE) UseMethod("postCP")

postCP.default <- function(data=numeric(),seg=integer(),model=1,lprob=numeric(),level.ind=numeric(),keep=TRUE,ci=0.9,viterbi=TRUE,initsegci=TRUE,nsamples=0,gen.data="n",prior=0.5,prior.type="n",epsilon=1e-9,disp.equal=TRUE,eps.nb=1e-8,verbose=TRUE,debug=FALSE) {
  if ((model!=1)&(model!=2)&(model!=3)){
    stop("Choose model=1 (Poisson) or 2 (normal) or 3 (negative binomial)")
  }
  n=length(data);

  nseg=length(seg)+1
  
  table.seg=table(seg)
  level.based=FALSE
  if (length(level.ind)>0) level.based=TRUE
  if (level.based) J=length(unique(level.ind))
  if (sum(table.seg>1)>0) {    
    stop(paste("Change-point",names(table.seg)[table.seg>1],"is repeated at least twice, please include only once"));
  }
  seg=sort(seg)
  if (length(lprob)>0 & length(data)>0){
    warning("Both data/segmentation/model and log-densities entered, log-densities will be used in HMM calculations");
  }
  if (length(lprob)==0){
  if (n<2) {
    stop("Not enough data to segment");
  }
  if (nseg>n) {
    stop("Number of segments larger than data");
  }
  if (length(seg)>0){
    if (max(seg)>=n) stop("At least one change-point initialized to larger than n");
	if (min(seg)<1)    stop("Location of first change-point must be at least 1");
	if (sum(round(seg)!=seg)>0) stop("Enter all change-points as integers");
  }
  if (viterbi==FALSE) initsegci=TRUE;
  if (initsegci==FALSE) viterbi=TRUE;
  if ((sum(data<0)>0)&model==1) {
    stop("Negative numbers with Poisson distribution specified, choose model=2");
  }
  if ((sum(round(data)!=data)>0)&model==1) {
    stop("Non-integer numbers with Poisson distribution specified, choose model=2");
  }
  } else{
  if (!is.matrix(lprob)) {
    stop("Enter lprob as a matrix, with n rows (one for each observation) and J columns (one for each hidden state)");
  }
  zero.prob=rowSums(lprob!=-Inf)==0
  if (sum(zero.prob)>0) {stop(paste("Observation in row(s)",paste(which(zero.prob),collapse=","),"in matrix lprob have zero probabilities for all states, please re-enter lprob matrix"))}

    n=nrow(lprob);
    J=ncol(lprob);
    nseg=J;
    lprob=as.vector(lprob);
    J1=length(unique(level.ind))
    if (level.based&(J!=J1)){
       stop("Number of different levels in level.ind does not match number of columns in lprob")
    }

  }
  if (!level.based) J=nseg else {
    
    nseg1=sum(diff(level.ind)!=0)+1
    if ((nseg1!=nseg)){
       if (length(lprob)==0) stop("Number of level changes in level.ind does not match number of change-points in seg")     
    }
  }
  if (model==3){
    if(disp.equal) {
	  if (level.based) stop("For negative binomial model with equal disperson, please select level.based=FALSE.");
          print("For negative binomial model with equal dispersion, segments will be chosen by PDPA method.");
	  segout=Segmentor(data,model=model,Kmax=max(nseg,2));
	  if (nseg>1) seg=(getBreaks(segout))[nseg,1:(nseg-1)];
	  param.nb=list(prob=getParameters(segout)[nseg,1:nseg],disp=getOverdispersion(segout));
      param.nb$prob=param.nb$prob+eps.nb;
	  param.nb$size = param.nb$disp/(1 - param.nb$prob+eps.nb)
      param.nb$mu = param.nb$disp * param.nb$prob/(1 - param.nb$prob+eps.nb)
	  lprob=matrix(dnbinom(rep(data,nseg),size=param.nb$size,mu=param.nb$mu[rep(1:nseg,each=n)],log=TRUE),ncol=nseg)
    }else{
      if (!level.based) state.temp=rep(1:nseg,diff(c(0,seg,n))) else state.temp=level.ind;
      out.mle=lapply(split(data,state.temp),mleNB,eps.nb);
      out.sizes=sapply(out.mle,function(x) x[1]);
      out.means=sapply(out.mle,function(x) x[2]);
      # log-density of neg binomial, 1 column for each possible state
      lprob=matrix(dnbinom(rep(data,J),size=out.sizes[rep(1:J,each=n)],mu=out.means[rep(1:J,each=n)],log=TRUE),ncol=J);
   #    level.based=FALSE;
    }
  }
  probs=(length(lprob)>0);
  if (sum(prior<0|prior>1)>0) {
    stop("Choose all priors between 0 and 1");
  }
  if (level.based){
         is.wholenumber <-
         function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

    if (length(level.ind)!=n) {
      stop("Enter a vector of level indices of length n");
    }
    if (min(level.ind)!=1) {
      stop("Level indices must start at one");
    }
    if (level.ind[1]!=1) {
      stop("Level index of first observation must be 1");
    }

    if (sum(!is.wholenumber(level.ind))>0) {
      stop("All level indices must be integers");
    }
    if (max(diff(sort(unique(level.ind))))>1) {
      stop("Set of level indices must be consecutive");
    }
  }
  if ((length(prior)!=n)&(prior.type=="o")) {
    stop("n transitions (1 for each observation) need to be chosen for prior vector");
  }
  if ((length(prior)!=nseg)&(prior.type=="s")) {
    stop(" transitions (1 for each state) need to be chosen for prior vector");
  }

 if (prior.type=="n") ptype=1 else
 if (prior.type=="o") ptype=2 else
 if (prior.type=="s") ptype=3 else
   stop("if prior.type specified, must be by (o)bservation or (s)tate");
 if (!level.based) level.ind=rep(1:nseg,diff(c(0,seg,n))) else{
    seg=which(diff(level.ind)!=0)
    nseg=length(seg)+1
 }
  lforward=0
   lbackward=0
 if (!keep){
   bestcp=rep(0,nseg-1)
   post.cp=0
 } else{
     if (!level.based){
       lforward=rep(0,n*J)
       lbackward=rep(0,n*J)
     } else{
       rm(lforward)
       rm(lbackward)
       lforward=rep(0,n*J*nseg)
       lbackward=rep(0,n*J*nseg)

     }
       post.cp=rep(0,(n-1)*(nseg-1))
       bestcp=rep(0,nseg-1)
    
 } 
if (level.based) best.level=rep(0,nseg)
 if (length(data)>0) {
   rmu=aggregate(data,list(level.ind),mean)$x
   rsigma=sqrt(sum((data-rmu[level.ind])^2)/n)
   if (model==1) rmu=rmu+epsilon
   if (model==2) rsigma=rsigma+epsilon
   } else{
      rmu=0; rsigma=0;rseg=0;rinitsegci=0;if (!level.based) rmu1=rmu else  rmu1=rmu[level.ind[diff(c(0,level.ind))!=0]]
}
 cpconfint=rep(0,(nseg-1)*3)
 evidence=0

 #  if ((sum(rmu==0)>0)&(model==1)&!probs) {
#    stop("One state consists entirely of zeroes, FB algorithm will not work with Poisson dist");
#  }
 

 cpvector=rep(0,nsamples*(nseg-1))
 if (level.based) levvector=rep(0,nsamples*(nseg-1))
 if (level.based){
 .C("postCPlev",data=as.double(data),rseg=as.integer(seg),nn=as.integer(n),JJ=as.integer(J),rnseg=as.integer(nseg), lforward=as.double(lforward),lbackward=as.double(lbackward),lprob=as.double(lprob),rprobs=as.integer(probs),cp=as.double(post.cp),bestcp=as.double(bestcp),bestlevel=as.double(best.level),rmu=as.double(rmu),rsigma=as.double(rsigma),cpconfint=as.double(cpconfint),cpvector=as.double(cpvector),levvector=as.double(levvector),rci=as.double(ci),rviterbi=as.integer(viterbi),rinitsegci=as.integer(initsegci),rnsamples=as.integer(nsamples),rmodel=as.integer(model),rprior=as.double(prior),rpriortype=as.integer(ptype),rmout=as.integer(keep),rlevelbased=as.integer(1),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")} else{
     .C("postCP",data=as.double(data),rseg=as.integer(seg),nn=as.integer(n),JJ=as.integer(J),lforward=as.double(lforward),lbackward=as.double(lbackward),cp=as.double(post.cp),bestcp=as.double(bestcp),rmu=as.double(rmu),rsigma=as.double(rsigma),lprob=as.double(lprob),rprobs=as.integer(probs),cpconfint=as.double(cpconfint),cpvector=as.double(cpvector),rci=as.double(ci),rviterbi=as.integer(viterbi),rinitsegci=as.integer(initsegci),rnsamples=as.integer(nsamples),rmodel=as.integer(model),rprior=as.double(prior),rpriortype=as.integer(ptype),rmout=as.integer(keep),rlevelbased=as.integer(0),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")} 
#else{
    #  .C("postCPext",lprob=as.double(lprob),nn=as.integer(n),JJ=as.integer(J),lforward=as.double(lforward),lbackward=as.double(lbackward),cp=as.double(post.cp),bestcp=as.double(bestcp),cpconfint=as.double(cpconfint),cpvector=as.double(cpvector),rci=as.double(ci),rnsamples=as.integer(nsamples),rprior=as.double(prior),rpriortype=as.integer(ptype),rmout=as.integer(keep),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")
  #  } 
# }
 if (!level.based) rmu1=rmu else  rmu1=rmu[level.ind[diff(c(0,level.ind))!=0]]
   if (nsamples>0) {
     cpmatrix=matrix(cpvector,ncol=nseg-1)
     fbsample <- list(changepoints=cpmatrix)
     names(fbsample) <- "changepoints"
     if (level.based) {
	levmatrix=matrix(levvector,ncol=nseg-1)
	cpmatrix=cbind(cpmatrix,levmatrix)
      }	
     if (gen.data=="p"|gen.data=="np"){
	getSegment <- function(cpv,n,level.based=FALSE) {
		if (!level.based) out=rep(1:(length(cpv)+1),diff(c(0,cpv,n))) else {  
			lv = c(1,cpv[(length(cpv)/2+1):length(cpv)])
			cpv=cpv[1:(length(cpv)/2)]
			out=rep(lv,diff(c(0,cpv,n)))
		}
  		return(out)}
		if (gen.data=="p"){
			seglabs <- apply(cpmatrix,1,getSegment,n=n,level.based=level.based) 
		if (model==1) x=matrix(rpois(length(seglabs),rmu1[seglabs]),ncol=n,byrow=T)
		if (model==2) x=matrix(rnorm(length(seglabs),rmu1[seglabs],rsigma),ncol=n,byrow=T)		
		} else{
			getIndex <- function(i,cpv,n,level.based) which(getSegment(cpv,n,level.based)==i)
			getSample <- function (cpv,n,level.based) {
				if (!level.based) levels=1:(length(cpv)+1) else levels=unique(c(1,cpv[(length(cpv)/2+1):length(cpv)]))
				ind.list=lapply(levels,getIndex,cpv,n=n,level.based=level.based) 
				ind.order=order(unlist(ind.list))
				getSample2 <- function(ind) {
					if (length(ind)>1) out=sample(ind,length(ind),replace=TRUE) else out=ind;
					return(out);
				}
				unlist(lapply(ind.list,getSample2))[ind.order]
			}	
			seglabs=matrix(apply(cpmatrix,1,getSample,n=n,level.based=level.based),ncol=nsamples)
			x=matrix(data[seglabs],ncol=n,byrow=T)
		}
	fbsample$x=x
	names(fbsample)[2] <- "data"
	}
     if (!is.element(gen.data,c("n","np","p"))) print("Please choose np (non-parametric resampling), p (parametric resampling) or n (no resampling)");
    
     cp.out=matrix(cpconfint,ncol=3)
     colnames(cp.out)=c("est",paste("lo.",ci,sep=""),paste("hi.",ci,sep=""))
     if (model==2) {
        model.dist="normal"
        postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu1,sds=rsigma,prior=prior,prior.type=ptype,fbsample=fbsample)
     } else{
       model.dist="";
       if (model==1) model.dist="Poisson";
       if (model==3) model.dist="negative Binomial";
       postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu1,prior=prior,prior.type=ptype,fbsample=fbsample)    
     }

   } else{ 
     cp.out=matrix(cpconfint,ncol=3)
     colnames(cp.out)=c("est",paste("lo.",ci,sep=""),paste("hi.",ci,sep=""))
     if (model==2) {
       model.dist="normal"
      postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu1,sds=rsigma,prior=prior,prior.type=ptype)
     }else{
         model.dist="";
        if (model==1)  model.dist="Poisson";
        if (model==3)  model.dist="negative Binomial";             
        postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu1,prior=prior,prior.type=ptype)
     }     
  } 
  postCP.res$nseg=nseg
  if (viterbi) {postCP.res$bestcp=bestcp
     if (level.based)    {
       best.level[1]=1
       postCP.res$best.level= rep(best.level,diff(c(0,bestcp,n)))
    }
  }
  
  if (keep){
     if (length(data)>0) { 
        postCP.res$data=data
        postCP.res$seg=seg
     }
     if (probs){
     postCP.res$lprob=matrix(lprob,ncol=J)
     }

   if (!level.based) {
        postCP.res$lforward=matrix(lforward,ncol=nseg)
        postCP.res$lbackward=matrix(lbackward,ncol=nseg)
   }else {
        segindmat=rep(rep(1:nseg,each=n),J)
        levindmat=rep(1:J,each=n*nseg)
        indmat=rep(rep(1:n),nseg*J)
	postCP.res$level.ind=level.ind
        postCP.res$lforward=lapply(split(lforward,list(levindmat)),matrix,ncol=nseg)
        postCP.res$lbackward=lapply(split(lbackward,list(levindmat)),matrix,ncol=nseg)
   }
     if (nseg>1) postCP.res$post.cp=matrix(post.cp,ncol=nseg-1) else postCP.res$post.cp=numeric()
     if (!level.based) postCP.res$post.state=exp(postCP.res$lforward+postCP.res$lbackward-postCP.res$lforward[1,1]-postCP.res$lbackward[1,1]) else{
     lsum<-function(lx) {
        max.l=which.max(lx)
         lx[max.l]+log(sum(exp(-abs(lx-lx[max.l]))))
     }
     levidence=lforward[1]+lbackward[1]
     lforbacsum=exp(lforward+lbackward-levidence)
     
     post.state=sapply(split(lforbacsum,list(indmat,segindmat)),sum) # posterior probability of state
     postCP.res$post.state= matrix(post.state,ncol=nseg)

     lforbacsum=exp(lforward+lbackward-levidence)
     post.level=sapply(split(lforbacsum,list(indmat,levindmat)),sum) # posterior probability of level
     postCP.res$post.level= matrix(post.level,ncol=J)
     }
  } else  postCP.res$seg=seg
  if (model==3) postCP.res$eps.nb=eps.nb;
  if (probs&model.dist=="") postCP.res$model="blank"
  postCP.res$call <- match.call()
  class(postCP.res) <- "postCP"
  postCP.res
}

print.postCP <- function(x,...)
{ cat("Call:\n")
  print(x$call)
  data.ent=is.element("data",names(x)) # T if data entered instead of log-densities
  if (!is.element("best.level",names(x))) {
     if (data.ent)  out.table=data.frame(seg=1:length(x$means),start=c(1,x$seg+1),end=c(x$seg,x$n),size=c(x$seg,x$n)-c(1,x$seg+1)+1,mean=x$means) else out.table=data.frame(seg=1:length(x$means),start=c(1,x$cp.est[,1]+1),end=c(x$cp.est[,1],x$n),size=c(x$cp.est[,1],x$n)-c(1,x$cp.est[,1]+1)+1)
  } else{
     if (data.ent)  out.table=data.frame(seg=1:length(x$means),level=x$level.ind[c(x$seg,x$n)],start=c(1,x$seg+1),end=c(x$seg,x$n),size=c(x$seg,x$n)-c(1,x$seg+1)+1,mean=x$means) else out.table=data.frame(seg=1:length(x$means),level=x$level.ind[c(x$cp.est,x$n)],start=c(1,x$cp.est[,1]+1),end=c(x$cp.est[,1],x$n),size=c(x$cp.est[,1],x$n)-c(1,x$cp.est[,1]+1)+1)
  }

  if (data.ent) {cat("\nDistribution:")
  print(x$model)}
  cat("\nSegment information:\n")
  print(out.table)
  if (!is.element("fbsample",names(x))) print(x[!is.element(names(x),c("call","data","lforward","lbackward","lprob","seg","means","post.cp","means","post.state","post.level","model","prior","prior.type","best.level"))])
  else {
    cat("\nEstimated change points (with ci):\n")
    print(x[!is.element(names(x),c("fbsample","data","level.ind","lprob","seg","means","lforward","lbackward","post.cp","post.state","post.level","model","call","prior","prior.type","best.level"))])
    cat("Generated change-points: ($fbsample$changepoints)\n")
    str(x$fbsample[[1]])
    if (length(x$fbsample)==2){
      cat("Generated data: ($fbsample$data)\n")
      str(x$fbsample[[2]])
    }
  }
  if (is.element("best.level",names(x))){
   cat("\nVector of indices of levels: ($best.level)")
      str(x$best.level)}
  if (is.element("post.cp",names(x))){
    if (data.ent){
      cat("\nData: ($data) ")
      str(x$data)
    } else{
      cat("\nLog-densities ($lprob) ")
      str(x$lprob)
    }    
    if (is.element("level.ind",names(x))){
      cat("\nLevel indices: ($level.ind)")
      str(x$level.ind)   
    } 
    cat("\nBackward matrix: ($lforward)")
    str(x$lbackward)
    cat("\nForward matrix: ($lbackward)")
    str(x$lforward)
    cat("\nChange-point matrix: ($post.cp)")
    str(x$post.cp)
    cat("\nState probability matrix: ($post.state)")
    str(x$post.state)
    if (is.element("post.level",names(x))){
     cat("\nLevel probability matrix: ($post.level)")
      str(x$post.level)   
    }
  }
  cat("\nPrior:\n")
  str(x$prior)
  cat("\nPrior type: ")
  if (x$prior.type==1) cat("homogeneous HMM\n")
  if (x$prior.type==2) cat("by observation\n")
  if (x$prior.type==3) cat("by hidden state\n")
}


viterbi <- function(data=numeric(),seg=integer(),model=1,lprob=numeric(),level.ind=numeric(),prior=0.5,prior.type="n",epsilon=1e-9,verbose=TRUE,debug=FALSE) {
  if (length(lprob)>0) model=0;
  if ((model!=1)&(model!=2)&(model!=3)&(length(lprob)==0)){
    stop("Choose model=1 (Poisson), 2 (normal) or 3 (negative binomial)")
  }
  n=length(data);

  nseg=length(seg)+1
  
  table.seg=table(seg)
  level.based=FALSE
  if (length(level.ind)>0) level.based=TRUE
  if (level.based) J=length(unique(level.ind))
  if (sum(table.seg>1)>0) {    
    stop(paste("Change-point",names(table.seg)[table.seg>1],"is repeated at least twice, please include only once"));
  }
  seg=sort(seg)
  if (length(lprob)>0 & length(data)>0){
    stop("Both data/segmentation/model and log-densities entered, specify only one");
  }
  if (length(lprob)==0){
  if (n<2) {
    stop("Not enough data to segment");
  }
  if (nseg>n) {
    stop("Number of segments larger than data");
  }
  if (length(seg)>0){
    if (max(seg)>=n) stop("At least one change-point initialized to larger than n");
	if (min(seg)<1)    stop("Location of first change-point must be at least 1");
	if (sum(round(seg)!=seg)>0) stop("Enter all change-points as integers");
  }
  
  if ((sum(data<0)>0)&model==1) {
    stop("Negative numbers with Poisson distribution specified, choose model=2");
  }
  if ((sum(round(data)!=data)>0)&model==1) {
    stop("Non-integer numbers with Poisson distribution specified, choose model=2");
  }
  } else{
  if (!is.matrix(lprob)) {
    stop("Enter lprob as a matrix, with n rows (one for each observation) and J columns (one for each hidden state)");
  }
  zero.prob=rowSums(exp(lprob),na.rm=TRUE)==0
 # if (sum(zero.prob)>0) {stop(paste("Observation in row(s)",paste(which(zero.prob),collapse=","),"in matrix lprob have zero probabilities for all states, please re-enter lprob matrix"))}

    n=nrow(lprob);
    J=ncol(lprob);
    nseg=J;
    lprob=as.vector(lprob);
    J1=length(unique(level.ind))
    if (level.based&(J!=J1)){
       stop("Number of different levels in level.ind does not match number of columns in lprob")
    }

  }

  level.based=FALSE
  if (length(level.ind)>0) level.based=TRUE
  if (level.based) J=length(unique(level.ind))

  if (!level.based) {
    J=nseg; 
  }        else {
    
    nseg1=sum(diff(level.ind)!=0)+1
    if ((nseg1!=nseg)){
       if (length(lprob)==0) stop("Number of level changes in level.ind does not match number of change-points in seg")     
    }
  }
    probs=(length(lprob)>0)
  if (sum(prior<0|prior>1)>0) {
    stop("Choose all priors between 0 and 1");
  }
  if (level.based){
         is.wholenumber <-
         function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

    if (length(level.ind)!=n) {
      stop("Enter a vector of level indices of length n");
    }
    if (min(level.ind)!=1) {
      stop("Level indices must start at one");
    }
    if (level.ind[1]!=1) {
      stop("Level index of first observation must be 1");
    }

    if (sum(!is.wholenumber(level.ind))>0) {
      stop("All level indices must be integers");
    }
    if (max(diff(sort(unique(level.ind))))>1) {
      stop("Set of level indices must be consecutive");
    }
  }
  if ((length(prior)!=n)&(prior.type=="o")) {
    stop("n transitions (1 for each observation) need to be chosen for prior vector");
  }
  if ((length(prior)!=nseg)&(prior.type=="s")) {
    stop(" transitions (1 for each state) need to be chosen for prior vector");
  }

 if (prior.type=="n") ptype=1 else
 if (prior.type=="o") ptype=2 else
 if (prior.type=="s") ptype=3 else
   stop("if prior.type specified, must be by (o)bservation or (s)tate");
 if (!level.based) level.ind=rep(1:nseg,diff(c(0,seg,n))) else{
    seg=which(diff(level.ind)!=0)
    nseg=length(seg)+1
 }
   bestcp=rep(0,nseg-1)
if (level.based) best.level=rep(0,nseg) else best.level=0
 if (!probs) {
   rmu=aggregate(data,list(level.ind),mean)$x
   rsigma=sqrt(sum((data-rmu[level.ind])^2)/n)
   if (model==1) rmu=rmu+epsilon
   if (model==2) rsigma=rsigma+epsilon
   } else{
      rmu=0; rsigma=0;rseg=0;
}

 evidence=0

 #  if ((sum(rmu==0)>0)&(model==1)&!probs) {
#    stop("One state consists entirely of zeroes, FB algorithm will not work with Poisson dist");
#  }
  probs=(length(lprob)>0);
 .C("viterbi2",data=as.double(data),lprob=as.double(lprob), rprobs=as.integer(probs),bestcp=as.double(bestcp),bestlevel=as.double(best.level),rprior=as.double(prior),rpriortype=as.integer(ptype),nn=as.integer(n),JJ=as.integer(J),rnseg=as.integer(nseg),rmu=as.double(rmu),rsigma=as.double(rsigma),rmodel=as.integer(model), level=as.integer(level.based),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")

#else{
    #  .C("postCPext",lprob=as.double(lprob),nn=as.integer(n),JJ=as.integer(J),lforward=as.double(lforward),lbackward=as.double(lbackward),cp=as.double(post.cp),bestcp=as.double(bestcp),cpconfint=as.double(cpconfint),cpvector=as.double(cpvector),rci=as.double(ci),rnsamples=as.integer(nsamples),rprior=as.double(prior),rpriortype=as.integer(ptype),rmout=as.integer(keep),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")
  #  } 
# }

  out=list()   
  out$bestcp=bestcp
   if (level.based)    {
       best.level[1]=1
       out$best.level= rep(best.level,diff(c(0,bestcp,n)))
    }
 out
}

