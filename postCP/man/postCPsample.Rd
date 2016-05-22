\name{postCPsample}
\alias{postCPsample}
\alias{print.postCPsample}
\alias{postCPsample.default}
\title{
  Resampling of change-point data.
}
\description{
  The functions are used for change-point problems, after an initial set of change-points within the data has already been obtained. The function postCPsample obtains random samples of sets of change-points using the output of the postCP function.
}
\usage{
postCPsample(postCP.res, nsamples=100, gen.data="n", 
prior=0.5,prior.type="n",verbose=TRUE,debug=FALSE) 
}
\arguments{
\item{postCP.res}{Results of postCP with keep=TRUE, use as input for postCPsample.}\item{nsamples}{Generated number of replicates of J-1 change-points.}
\item{gen.data}{Generate matrix of segmented data, each row corresponding to a different replicate. "p" for parametric resampling, "np" for non-parametric re-sampling, "n" for no resampling. Parametric resampling is only done for Poisson or normally distributed data.}
\item{prior}{Transition probabilities chosen as priors for HMM model. For homogeneous HMM, enter a scalar between 0 and 1 (default: 0.5). For heterogeneous HMM, enter a vector of probabilities, of length n for n observations or J for J hidden states.}
\item{prior.type}{Type of prior entered. For "n":homogeneous HMM (default), "o": vector of probabilities for each observation, "s": vector of probabilities for each hidden state.}
\item{verbose}{Additional screen output, default=TRUE.}
\item{debug}{Additional screen output for debugging purposes, default=FALSE.}
}

\value{
\item{changepoints}{Matrix of generated changepoint locations of nsamples rows, each row containing J-1 changepoints.}
\item{data}{Matrix of generated data of nsamples rows, each row containing vector of n observations.}
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

}

\author{
Gregory Nuel and The Minh Luong

Maintainer: Who to complain to <the-minh.luong@parisdescartes.fr>
}

\references{
  Luong, T.M., Rozenholc, Y. & Nuel, G. (2012). Fast estimation of posterior probabilities in change-point models through a constrained hidden Markov model. \url{http://arxiv.org/pdf/1203.4394}
}

\keyword{ package }
\seealso{
\code{\link{postCP}}
}

\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
# obtain 85% confidence intervals for changepoints
res=postCP(data=x,seg=c(N,3*N),keep=TRUE, ci=0.85);  

## Use results from forward-backward to get 
## posterior probabities of a change-point at each location
## draw 100 samples of change-points using results from postCP
## also a matrix with 100 rows of generated data of length 8000
res.samples=postCPsample(res,100,gen.data="np");
## draw 1000 samples in same function call as segmenting data
# postCP(data=x,seg=c(N,3*N),nsamples=1000,gen.data="np");
}
