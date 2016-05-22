\name{lesum}
\alias{lesum}
\title{
Log of sum of exponentials.
}
\description{
  Calculates log of sum of exponentials. Useful in problems with underflow issues due to very low values.}
\usage{
lesum(lx)
}
\arguments{
\item{lx}{Vector of values, usually log-transformed data whose exponentials need to be summed.}
}

\value{
\item{out}{Log of sum of exponential of values.}
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
  Our manuscript.
}

\keyword{ package }
\seealso{
\code{\link{postCP}}
}

\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
lx=log(x);
lesum(lx);
log(sum(exp(lx)));
}
