\name{postCP}
\alias{postCP}
\alias{print.postCP}
\alias{postCP.default}
\title{
  Characterize uncertainty in change-point estimates.
}
\description{
  The functions are used for change-point problems, after an initial set of change-points within the data has already been obtained. The function postCP obtains estimates of posterior probabilities of change-points.
}
\usage{
postcp(formula, data=numeric(), bp=integer(), family = gaussian, par=numeric(), homoscedastic=TRUE,  sigma=1.0)

}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are as given in glm.}
\item{data}{A matrix of data to be segmented. Must have no missing values}
\item{bp}{break points. i.e. For K segments, a vector of K-1 initial estimates of change-point locations, corresponding to the index of the last observation of the first K-1 segments.}
\item{family}{Emission distribution of the observed data as in glm. Supports Poisson(default), Gaussian, Negative Binomial and Gamma distributions.}
\item{par}{A matrix of initial values of the regression coefficients}
\item{homoscedastic}{Applies for the homoscedastic gaussian model}
\item{sigma}{ Value added to sigma for gaussian data}
}


\value{
\item{model}{Emission distribution (Poisson or normal).} 
\item{n}{Length of segmented data}
\item{loglik}{loglikelihood matrix}
\item{post}{Posterior probability distribution}
\item{post.cp}{Posterior change point probability distribution}
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

\examples{
require(postCP)
set.seed(18)
#prepare data 
n=300
sigma=1.3
#Initial regression coefficients 
beta_ref=matrix(c(-0.7,0.3,0.2,-0.5,0.7,0.3),nrow=2,byrow=TRUE)  
cp=c(100,200)
num=runif(n)
#Generate a matrix of data
dd=data.frame(num)

#Apply postcp function 
res = postcp(~1+num,family=gaussian,data=dd,bp=c(100,200),par=beta_ref,sigma)

#Plot the results
plot.postCP.dev(res,main="Posterior Change Point Distribution")

}


