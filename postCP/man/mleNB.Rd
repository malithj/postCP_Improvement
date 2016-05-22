\name{mleNB}
\alias{mleNB}
\title{
Maximum likelihood estimation of negative binomial data.
}
\description{
  Find maximum likelihood estimates of size and mean parameters assuming a negative binomial distribution. Uses optim function, using Broyden-Fletcher-Goldfarb-Shanno (BFGS) method.}
\usage{
mleNB(x,eps.nb=1e-8)
}
\arguments{
\item{x}{Vector of data with hypothesized negative binomial distribution.}
\item{eps.nb}{Relative tolerance of BFGS method, increase if there are convergence problems.}
}

\value{
\item{size}{Size parameter MLE of negative binomial.}
\item{mu}{Mu (mean) parameter MLE of negative binomial.}
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
Venables, W.N. and Ripley, B.D. (2002) Modern Applied Statistics with S. Fourth edition. Springer. 
  Our manuscript.
}

\keyword{ package }
\seealso{
\code{\link{postCP}}
\code{\link{fitdistr}}
\code{\link{dnbinom}}
}

\examples{
require(postCP);
N=2000; x=rnbinom(N,size=10.0,mu=5.0);
mleNB(x);
}

