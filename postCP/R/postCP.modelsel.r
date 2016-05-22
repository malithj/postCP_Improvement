
lesum<-function(lx) { # finds log of sum of exponentials
        max.l=which.max(lx);
        out=lx[max.l]+log(sum(exp(-abs(lx-lx[max.l]))));
        return(out);
     }

mleNB <- function(x,eps.nb=1e-8){
if ((length(x)>1)&(mean(x)>0)) out=try(fitdistr(x,"Negative Binomial",control=list(reltol=eps.nb))$estimate) else out=c(mean(x),mean(x))+eps.nb;
if (class(out)=="try-error") out=c(mean(x),mean(x))+eps.nb;
return(out);
}

postCPcrit<-function(data,seg=numeric(),model,disp.equal=TRUE,eps.nb=1e-8, param.nb=NULL, prior=0.5, prior.type="n"){

  k=length(seg)+1;
  print(paste("Number of segments:",k));
  
  n=length(data);
  bestcp=numeric();
  table.seg=table(seg)
  if (sum(table.seg>1)>0) {    
    stop(paste("Change-point",names(table.seg)[table.seg>1],"is repeated at least twice, please include only once"));
  }
  seg=sort(seg)
  if (n<2) {
    stop("Not enough data to segment");
  }
  if (k>n) {
    stop("Number of segments larger than data");
  }
  if (length(seg)>0){
    if (max(seg)>=n) stop("At least one change-point initialized to larger than n");
	if (min(seg)<1)    stop("Location of first change-point must be at least 1");
	if (sum(round(seg)!=seg)>0) stop("Enter all change-points as integers");
  }
  if (!is.element(model,1:3)) stop("Please enter model= 1 (Poisson), 2 (normal), 3 (negative binomial)");
  if (model==3&disp.equal) {
    if(class(param.nb)=="NULL"){
      print("For negative binomial model with equal dispersion, segments will be chosen by PDPA method.");
	  segout=Segmentor(data,model=model,Kmax=max(k,2));
	  if (k>1) seg=(getBreaks(segout))[k,1:(k-1)];
  	  param.nb=list(prob=getParameters(segout)[k,1:k],disp=getOverdispersion(segout));
	}
	param.nb$prob=param.nb$prob+eps.nb;
	 param.nb$size = param.nb$disp/(1 - param.nb$prob+eps.nb)
     param.nb$mu = param.nb$disp * param.nb$prob/(1 - param.nb$prob+eps.nb)
  }  
  if (k==1){
    out=postCP(data,seg,model,keep=TRUE,verbose=FALSE,initsegci=FALSE,ci=0,disp.equal=disp.equal,prior=prior,prior.type=prior.type);
    if (model==2){
      mBIC=0;
      BIC=-sum(dnorm(data,mean(data),sd(data)*sqrt((n-1)/n),TRUE))+(2)*log(n);
      entropy=0;
      AIC=-sum(dnorm(data,mean(data),sd(data)*sqrt((n-1)/n),TRUE))+(2)*2;
    } 
    if (model==1){
      mBIC=-n*mean(data)*log(mean(data))+log(n);
      BIC=-sum(dpois(data,mean(data),TRUE))+log(n);
      AIC=-sum(dpois(data,mean(data),TRUE))+2;
      entropy=0;
    } 
    if (model==3){
      if (!disp.equal){
        negbin.est=mleNB(data,eps.nb=eps.nb);
        BIC=-sum(dnbinom(data,mu=negbin.est[2],size=negbin.est[1],log=TRUE))+1*log(n);
        AIC=-sum(dnbinom(data,mu=negbin.est[2],size=negbin.est[1],log=TRUE))+1*2;
      }else{
	    BIC=-sum(dnbinom(data,size=param.nb$size,mu=param.nb$mu,log=TRUE))+(1)*log(n);
        AIC=-sum(dnbinom(data,size=param.nb$size,mu=param.nb$mu,log=TRUE))+1*2; 
        lprob.matrix=matrix(dnbinom(rep(data,1),size=param.nb$size,mu=param.nb$mu,log=TRUE),ncol=1);	   
        out=postCP(lprob=lprob.matrix,keep=TRUE,verbose=FALSE,viterbi=FALSE,initsegci=FALSE,ci=0,prior=prior,prior.type=prior.type);
      }
      entropy=0;
    } 
  }else{
   if (model==3){
    state.temp=rep(1:k,diff(c(0,seg,n)));
    if (!disp.equal){
      out.mle=lapply(split(data,state.temp),mleNB,eps.nb=eps.nb);
      out.sizes=sapply(out.mle,function(x) x[1]);
      out.means=sapply(out.mle,function(x) x[2]);
      # log-density of neg binomial, 1 column for each possible state
      lprob.matrix=matrix(dnbinom(rep(data,k),size=out.sizes[rep(1:k,each=n)],mu=out.means[rep(1:k,each=n)],log=TRUE),ncol=k);
      # initial run of postCP
      out=viterbi(lprob=lprob.matrix,verbose=FALSE,prior=prior,prior.type=prior.type);
      bestcp=out$bestcp;
      state.temp=rep(1:k,diff(c(0,out$bestcp,n)));
      out.mle=lapply(split(data,state.temp),mleNB,eps.nb=eps.nb);
      out.sizes=sapply(out.mle,function(x) x[1]);
      out.means=sapply(out.mle,function(x) x[2]);
      lprob.matrix=matrix(dnbinom(rep(data,k),size=out.sizes[rep(1:k,each=n)],mu=out.means[rep(1:k,each=n)],log=TRUE),ncol=k);
      rm(out);
      }else{
	    lprob.matrix=matrix(dnbinom(rep(data,k),size=param.nb$size,mu=param.nb$mu[rep(1:k,each=n)],log=TRUE),ncol=k);	           
	    bestcp=seg;
    }    
    out=postCP(lprob=lprob.matrix,keep=TRUE,verbose=FALSE,viterbi=FALSE,initsegci=FALSE,ci=0,prior=prior,prior.type=prior.type);
    if (!disp.equal) BIC=-sum(dnbinom(data,size=out.sizes[state.temp],mu=out.means[state.temp],log=TRUE))+(k*1)*log(n) else BIC=-sum(dnbinom(data,size=param.nb$size,mu=param.nb$mu[state.temp],log=TRUE))+(k*1)*log(n);
    if (!disp.equal) AIC=-sum(dnbinom(data,size=out.sizes[state.temp],mu=out.means[state.temp],log=TRUE))+(k*1)*2 else AIC=-sum(dnbinom(data,size=param.nb$size,mu=param.nb$mu[state.temp],log=TRUE))+(k*1)*2;
   }else{
      out=viterbi(data,seg,model,verbose=FALSE,prior=prior,prior.type=prior.type);
      bestcp=out$bestcp;
     # redo postCP using a posteriori bestCP (viterbi) from previous postCP run
      state.temp=rep(1:k,diff(c(0,bestcp,n)));
      out.means=aggregate(data,list(state.temp),mean)$x;
      out=postCP(data,bestcp,model,keep=TRUE,verbose=FALSE,debug=FALSE,ci=0,prior=prior,prior.type=prior.type);
      if (model==2) {
        # use for normal data
	    SSB=sum(table(state.temp)*(out.means-mean(data))^2);
        SSA=sum((data-mean(data))^2);
        mBIC=-(0.5*(n-k+2)*log(1+SSB/(SSA-SSB))+lgamma(0.5*(n-k+2))-lgamma(0.5*(n+1))+0.5*(k-1)*log(SSA)-0.5*sum(log(table(state.temp)))+(1.5-k)*log(n));
        BIC=-sum(dnorm(data,out.means[state.temp],out$sds,TRUE))+(k+1)*log(n);
        AIC=-sum(dnorm(data,out.means[state.temp],out$sds,TRUE))+(k+1)*2;
      }
      if (model==1){
        # use for Poisson data
        mBIC=-(sum(table(state.temp)*out.means*log(out.means))-0.5*sum(log(table(state.temp)))+(0.5-k)*log(n));
        BIC=-sum(dpois(data,out.means[state.temp],TRUE))+(k)*log(n);
        AIC=-sum(dpois(data,out.means[state.temp],TRUE))+(k)*2;
      }
    }
    if (is.element("lprob",names(out))){
	ldensmatrix=out$lprob;
    } else{
    if (model==1){
       ldensmatrix=matrix(dpois(rep(data,k),out$means[rep(1:k,each=n)],log=TRUE),ncol=k);
    }
    if (model==2){
       ldensmatrix=matrix(dnorm(rep(data,k),out$means[rep(1:k,each=n)],out$sds,log=TRUE),ncol=k);
    } 
    if (model==3){
       if (!disp.equal) ldensmatrix=matrix(dnbinom(rep(data,k),size=out.sizes[rep(1:k,each=n)],mu=out.means[rep(1:k,each=n)],log=TRUE),ncol=k)else 
       ldensmatrix=matrix(dnbinom(rep(data,k),size=param.nb$size,mu=param.nb$mu[rep(1:k,each=n)],log=TRUE),ncol=k);	  
    }    
   }	 
   if (out$prior.type==1) prior=out$prior;
   if (out$prior.type==2) prior=out$prior;
   if (out$prior.type==3) prior=matrix(rep(out$prior,each=n-1),ncol=k-1);
   lpost.trans=ldensmatrix[-1,-1]+out$lbackward[-1,-1]-out$lbackward[-n,-k]+log(prior);   
   lpost.transx=ldensmatrix[-1,-k]+out$lbackward[-1,-k]-out$lbackward[-n,-k]+log(1-prior);
   post.ncp=(out$post.state[-n,-k]-out$post.cp);
   post.ncp[post.ncp<1e-8]=0;
   lpost.trans[lpost.trans<(-1000)]=NA;
   lpost.transx[lpost.transx<(-1000)]=NA;
   entropy=-sum(out$post.state[1,1]*log(out$post.state[1,1]),na.rm=TRUE)-sum(out$post.cp*lpost.trans,na.rm=TRUE)-sum(post.ncp*lpost.transx,na.rm=TRUE);
  }
 
  # postCPS=postCP(rep(0,n),seg=seg,model=1,viterbi=FALSE,prior=prior,prior.type=prior.type);
   postCPK=postCP(lprob=matrix(rep(0,n*k),ncol=k),viterbi=FALSE,prior=prior,prior.type=prior.type);
  ICL=-out$lforward[1,1]-out$lbackward[1,1]+postCPK$lforward[1,1]+postCPK$lbackward[1,1]+log(choose(n-1,k-1))+entropy;
  if (model!=3) { scores=c(ICL,AIC,BIC,mBIC); names(scores)=c("ICL","AIC","BIC","mBIC")} else {scores=c(ICL,AIC,BIC); names(scores)=c("ICL","AIC","BIC")};
  results<-list(scores=scores,cp.loc=bestcp);
  results;
}

postCPmodelsel <- function(data,K.range,model=1,greedy=FALSE,disp.equal=TRUE,eps.nb=1e-8,prior=0.5,prior.type="n"){

  n=length(data);
  K.min=K.range[1];
  K.max=K.range[length(K.range)];
  entropy=numeric();
  BIC=numeric();
  mBIC=numeric();
  AIC=numeric();
  bestcp=list();
  k=K.min;
  if (length(K.range)<2) stop("Please choose models with at least two different values of K.");
  if (!greedy) {
	  segout=Segmentor(data,model=model,Kmax=K.max);
	  seg.matrix=getBreaks(segout);
  }  else segout=NULL;

  modelsel1<- function(k,data,segout, model, disp.equal,eps.nb = eps.nb, prior, prior.type){
    errorinseg=FALSE;
  	 if (k>1) {
	    if (class(segout)=="NULL") seg=GreedySegmente(data,k)[2:k]-1 else  seg=(getBreaks(segout))[k,1:(k-1)];
		seg=sort(seg);
		table.seg=table(seg);
        if ((max(table(seg))>1)| (max(seg)>=n)|(min(seg)<1)|(sum(round(seg)!=seg)>0)){
           errorinseg=TRUE;
		   if (model!=3) { scores=rep(NA,4); names(scores)=c("ICL","AIC","BIC","mBIC")} else {scores=rep(NA,3); names(scores)=c("ICL","AIC","BIC")};
            out=list(scores=scores,cp.loc=rep(NA,k-1));
        }
	} else seg=numeric();
	if (!errorinseg){	  
	  if ((class(segout)!="NULL"&(model==3))&disp.equal) param.nb=list(prob=getParameters(segout)[k,1:k],disp=getOverdispersion(segout)) else param.nb=NULL;
	  out=postCPcrit(data,seg,model,disp.equal,eps.nb, param.nb, prior, prior.type);
	 }
	return(out);
  }
  out=sapply(K.range,modelsel1,data=data,segout=segout,model=model,disp.equal=disp.equal,eps.nb=eps.nb, prior=prior, prior.type=prior.type);
  scores=unlist(out[1,]);
  ICL=scores[names(scores)=="ICL"];
  BIC=scores[names(scores)=="BIC"];
  AIC=scores[names(scores)=="AIC"];
  names(ICL)=K.range;
  names(AIC)=K.range;
  names(BIC)=K.range;
  if (model!=3) {
    mBIC=scores[names(scores)=="mBIC"];    names(mBIC)=K.range;
  }
  
  best.K.ICL <- which.min(ICL);
  best.K.AIC <- which.min(AIC);
  best.K.BIC <- which.min(BIC);
  bestcp=out[2,];
  cp.loc=list();
  cp.loc$ICL = bestcp[[best.K.ICL]];
  cp.loc$BIC = bestcp[[best.K.BIC]];
  cp.loc$AIC = bestcp[[best.K.AIC]];
  if (model!=3) { 
	best.K.mBIC=which.min(mBIC); 
	} else best.K.mBIC=NA;
  if (model!=3) cp.loc$mBIC = bestcp[[best.K.mBIC]];
   
  if (model!=3) { scores=cbind(ICL,AIC,BIC,mBIC); colnames(scores)=c("ICL","AIC","BIC","mBIC");rownames=K.range;} else {scores=cbind(ICL,AIC,BIC); colnames(scores)=c("ICL","AIC","BIC");rownames=K.range;};
  if (model!=3) results<-list(ICL=K.range[best.K.ICL],AIC=K.range[best.K.AIC],BIC=K.range[best.K.BIC],mBIC=K.range[best.K.mBIC],scores=scores,cp.loc=cp.loc)else{
    results<-list(ICL=K.range[best.K.ICL],AIC=K.range[best.K.AIC],BIC=K.range[best.K.BIC],scores=scores,cp.loc=cp.loc,eps.nb=eps.nb)
  }
  results;
}
