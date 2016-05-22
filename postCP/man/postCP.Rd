\name{postCP}
\alias{postCP}
\alias{print.postCP}
\alias{postCP.default}
\title{
  Characterize uncertainty in change-point estimates.
}
\description{
  The functions are used for change-point problems, after an initial set of change-points within the data has already been obtained. The function postCP obtains estimates of posterior probabilities of change-point and hidden states for each observation, and confidence intervals for the positions of the change-point.
}
\usage{
postCP(data=numeric(),seg=integer(),model=1,lprob=numeric(),level.ind=numeric(),
keep=TRUE,ci=0.9,viterbi=TRUE,initsegci=TRUE,nsamples=0,gen.data="n",prior=0.5,
prior.type="n",epsilon=1e-9,disp.equal=TRUE,eps.nb=1e-8,verbose=TRUE,debug=FALSE)
}
\arguments{
\item{data}{A vector of observations to be segmented. Must have no missing values.}
\item{seg}{For K segments, a vector of K-1 initial estimates of change-point locations, corresponding to the index of the last observation of the first K-1 segments.}
\item{model}{Emission distribution of observed data, 1: Poisson (default), 2: normal, 3: negative binomial.}
\item{level.ind}{Optional for level-based model, where distributions may be shared within non-adjacent segments. Vector of indices of length n, each index corresponding to observations with the same mean. Index of first observation must be 1; indices must be integers and range from 1:J for J different levels (observations in different segments may share same levels). Observations in adjacent segments may not share the same mean.}
\item{lprob}{Instead of data, seg, and model, enter n*K matrix of log-densities for emission distribution, each row i corresponds to each observations, and each column k is the log-density of observation in row i in state k}
\item{keep}{Saves raw data, forward and backward matrices, and matrix of posterior change-point probabilities. Default=TRUE}
\item{ci}{Confidence interval displayed for each change-point location.}
\item{viterbi}{Obtain estimate of best set of hidden states by Viterbi algorithm, default=TRUE (always true if initsegci=FALSE)}
\item{initsegci}{If true, estimate change-point confidence interval around initial segmentation (seg), if false estimate confidence intervals around observations with highest posterior probabilities (found by viterbi algorithm). Default=TRUE, always FALSE if lprob matrix is specified.}
\item{nsamples}{Generated number of replicates of K-1 change-points.}
\item{gen.data}{Generate matrix of segmented data, each row corresponding to a different replicate. "p" for parametric resampling, "np" for non-parametric re-sampling, "n" for no resampling.}
\item{prior}{Transition probabilities chosen as priors for HMM model. For homogeneous HMM, enter a scalar between 0 and 1 (default: 0.5). For heterogeneous HMM, enter a vector of probabilities, of length n for n observations or K for K hidden states.}
\item{prior.type}{Type of prior entered. For "n": homogeneous HMM (default), "o": vector of probabilities for each observation, "s": vector of probabilities for each hidden state.}
\item{epsilon}{Value added to mu (for Poisson data) or sigma (for normal data) to enable forward-backward in case of variance estimates of 0.}
\item{disp.equal}{Equal dispersion parameter across all segments for negative binomial (defaults to TRUE). If disp.equal=TRUE, seg will automatically be chosen by the PDPA algorithm from Segmentor function.}
\item{eps.nb}{Relative tolerance of BFGS method for finding MLE of negative Binomial data, increase if there are convergence problems.}
\item{verbose}{Additional screen output, default=TRUE.}
\item{debug}{Additional screen output for debugging purposes, default=FALSE.}
}

\value{
\item{model}{Emission distribution (Poisson or normal).} 
\item{n}{Length of segmented data.}
\item{nseg}{Number of segments in data.}
\item{cp.est}{Matrix of estimated change-points, with lower and upper bounds of confidence intervals (given in column names).}
\item{means}{Mean of data in each of K segments.}
\item{sds}{For normal data, pooled standard deviation.}
\item{fbsample}{If nsamples>0, information for randomly generated change-points and data.}
\item{fbsample$changepoints}{Matrix of generated changepoint locations of nsamples rows, each row containing K-1 changepoints.}
\item{fbsample$data}{Matrix of generated data of nsamples rows, each row containing vector of n observations.}
\item{data}{If keep=TRUE, the segmented data (if entered).}
\item{lprob}{If keep=TRUE, the log-densities (if entered).}
\item{lforward}{If keep=TRUE, matrix of log of forward probabilities, n rows corresponding to observations and K columns corresponding to states. If using level-based model, lforward will be a list with one matrix corresponding to each level.}
\item{lbackward}{If keep=TRUE, matrix of log of backward probabilities, n rows corresponding to observations and K columns corresponding to states. If using level-based model, lbackward will be a list with one matrix corresponding to each level.}
\item{post.cp}{If keep=TRUE, matrix of posterior change-point probabilities, prob at row i column k is the probability of the kth change-point occurring at observation i.}
\item{post.state}{If keep=TRUE, matrix of posterior hidden state probabilities, prob at row i column k is the probability of observation i being in state k.}
\item{post.level}{If keep=TRUE, matrix of posterior mean probabilities, prob at row i column j is the probability of observation i having mean j.}
\item{prior}{Transition probabilities chosen as priors for HMM mode.}
\item{prior.type}{Type of prior entered.}
\item{bestcp}{Vector of set of $K-1$ change-points with the highest posterior probability (found by Viterbi algorithm).}
\item{best.level}{Vector of set of $n$ indices corresponding to levels with the highest posterior probability for each observation (found by Viterbi algorithm).}
\item{eps.nb}{Relative tolerance of BFGS method used for finding MLE of negative Binomial data.}
}

\details{
\tabular{ll}{
Package: \tab postCP\cr
Type: \tab Package\cr
Version: \tab 1.7.2\cr
Date: \tab 2013-11-25\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
Typical usage is todo!

Estimates posterior probabilities of change in distribution and state probabilities in segmentation problems through the forward-backward algorithm using a restricted hidden Markov model.
}

\author{
Gregory Nuel and The Minh Luong

Maintainer: Who to complain to <the-minh.luong@parisdescartes.fr>
}

\references{
  Luong, T.M., Rozenholc, Y. & Nuel, G. (2012). Fast estimation of posterior probabilities in change-point models through a constrained hidden Markov model. \url{http://arxiv.org/pdf/1203.4394}\cr
  Rabiner, L.R. (1989) A tutorial on hidden Markov models and selected applications in speech recognition. 77(2):257-286.\cr
  Viterbi, A.J. (1967) Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. IEEE Transactions on Information Theory. 13(2): 260-269.
}


\keyword{ package }

\seealso{
\code{\link{postCPsample}}
\code{\link{postCPcrit}}
\code{\link{postCPmodelsel}}
}
\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
# obtain 85% confidence intervals for changepoints, keep forward and backward matrices
res=postCP(data=x,seg=c(N,3*N),keep=TRUE, ci=0.85);  


## alternate method, use matrix of log-densities, n*J or 8000*3 matrix 
lprob.matrix=matrix(dpois(rep(x,3),res$means[rep(1:3,each=8000)],log=TRUE),ncol=3)
res1=postCP(lprob=lprob.matrix,keep=TRUE, ci=0.85);

## select observations in 1st and 3rd segments to have the same mean
reslev=postCP(data=x,seg=c(N,3*N),keep=TRUE, ci=0.85, level.ind=c(rep(1,N),rep(2,2*N),rep(1,N)));  

## Use results from forward-backward algorithm to obtain 
## posterior probabities of being a change-point at each location
par(mfrow=c(2,2));
plot(x,main="original data");
plot(res$post.cp[,1],col="blue",t='l',lwd=2,main="change point distributions"); 
points(res$post.cp[,2],col="red",t="l",lwd=2);
## Use results to obtain probability of state membership at each point
marg.prob=exp(res$lbackward+res$lforward-res$lbackward[1,1]-res$lforward[1,1]);
plot(marg.prob[,1],col="blue",t='l',lwd=2,main="marginal distributions"); 
points(marg.prob[,2],col="red",t="l",lwd=2); points(res$marginal[,3],col="green",t="l",lwd=2);

## draw 100 samples of change-points using results from postCP, 
## and also a matrix with 100 rows of resampled data of length 8000
res.samples=postCPsample(res,100,gen.data="np");
## draw 1000 samples in same function call as data segmentation
# postCP(data=x,seg=c(N,3*N),nsamples=1000,gen.data="np");

# Use heterogeneous with transition probabilities 0.4 and 0.7 for 1 to 2 and 2 to 3 transitions, 
# respectively, use 0.5 for transition to junk state K+1 from constrained HMM
postCP(data=x,seg=c(N,3*N),keep=TRUE, ci=0.85,prior=c(0.4,0.7,0.5),prior.type="s");

# Use posterior probability of change-points at each observation as prior transition probabilities, 
# use 0.5 for transition to junk state from constrained HMM
postCP(data=x,seg=c(N,3*N),keep=TRUE, ci=0.85,prior=c(rowSums(res$post.cp),0.5),prior.type="o");

}


