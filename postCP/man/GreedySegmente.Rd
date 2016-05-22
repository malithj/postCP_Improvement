\name{GreedySegmente}
\alias{GreedySegmente}
\title{
  Greedy segmentation of change-point data.
}
\description{
Uses greedy algorithm to segment data of length n, using locally most optimal change-points that minimize the least squares. Obtains a vector of first observation of all k segments, and n+1/}
\usage{
GreedySegmente(Y,binmax) 
}
\arguments{
\item{Y}{A vector of observations to be segmented. Must have no missing values.}
\item{binmax}{Number of segments (k) chosen, which is number of changepoints+1.}
}

\value{
\item{breaks}{A vector of estimated break points, of k+1 indices corresponding to first observation in all k segments, and n+1. To use with postCP, convert to change-point locations using breaks[2:k]-1.}
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
Gregory Nuel, The Minh Luong, Mael Thepaut, and Yves Rozenholc

Maintainer: Who to complain to <the-minh.luong@parisdescartes.fr>
}

\references{
  Luong, T.M., Rozenholc, Y. & Nuel, G. (2012). Fast estimation of posterior probabilities in change-point models through a constrained hidden Markov model. \url{http://arxiv.org/pdf/1203.4394}\cr
  Hartigan, J.A., Wong, M.A. (1979) Algorithm AS 136: A K-Mean Clustering Algorithm. Journal of the Royal Statistical Society, Series C (Applied Statistics) 28(1):100-108.
}

\keyword{ package }
\seealso{
\code{\link{postCP}}
\code{\link{postCPmodelsel}}
}

\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
firstobs=GreedySegmente(x,3); # obtain index of first observation of all k=3 segments, and n+1
cp=firstobs[2:3]-1; # obtain index of last observation of k-1=2 change-points to enter into postCP
postCP(x,cp);
}

