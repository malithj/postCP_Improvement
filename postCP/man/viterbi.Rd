\name{viterbi}
\alias{viterbi}
\title{
 Most probable set of change-points through Viterbi algorithm.
}
\description{ Obtains 'a posteriori' most probable set of change-points for segmentation problems. For level-based models, also obtains most 'a posteriori' set of most probable levels for each observation.
}
\usage{
viterbi(data=numeric(),seg=integer(),model=1,lprob=numeric(),level.ind=numeric(),
prior=0.5,prior.type="n",epsilon=1e-9,verbose=TRUE,debug=FALSE)
}
\arguments{
\item{data}{A vector of observations to be segmented. Must have no missing values.}
\item{seg}{For K segments, a vector of K-1 initial estimates of change-point locations, corresponding to the index of the last observation of the first K-1 segments.}
\item{model}{Emission distribution of observed data, 1: Poisson (default), 2: normal.}
\item{lprob}{Instead of data, seg, and model, enter n*K matrix of log-densities for emission distribution, each row i corresponds to each observations, and each column k is the log-density of observation in row i in state k.}
\item{level.ind}{Vector of indices of length n, each index corresponding to observations with the same mean. Index of first observation must be 1; indices must be integers and range from 1:J for J different levels (observations in different segments may share same levels). Observations in adjacent segments may not share the same mean.}
\item{prior}{Transition probabilities chosen as priors for HMM model. For homogeneous HMM, enter a scalar between 0 and 1 (default: 0.5). For heterogeneous HMM, enter a vector of probabilities, of length n for n observations or J for J hidden states.}
\item{prior.type}{Type of prior entered. For "n":homogeneous HMM (default), "o": vector of probabilities for each observation, "s": vector of probabilities for each hidden state.}
\item{epsilon}{Value added to mu (for Poisson data) or sigma (for normal data) to enable forward-backward in case of variance estimates of 0.}
\item{verbose}{Additional screen output, default=TRUE.}
\item{debug}{Additional screen output for debugging purposes, default=FALSE.}
}

\value{
\item{bestcp}{Vector of set of $K-1$ change-points with the highest posterior probability (found by Viterbi algorithm).}
\item{best.level}{Vector of set of $n$ indices corresponding to levels with the highest posterior probability for each observation (found by Viterbi algorithm).}
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
Viterbi, A.J. (1967) Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. IEEE Transactions on Information Theory. 13(2): 260-269.
}

\keyword{ package }
\seealso{
\code{\link{postCP}}
}

\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
# obtain a posteriori most probable set of change-points by Viterbi algorithm
viterbi(data=x,seg=c(N,3*N));  


}
