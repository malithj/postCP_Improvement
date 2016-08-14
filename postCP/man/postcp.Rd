\name{postCP}
\alias{postCP}
\alias{print.postCP}
\alias{postCP.default}
\title{
  Characterize uncertainty in change-point estimates.
}
\description{
  The functions are used for change-point problems, after an initial set of change-points within the data has already been obtained. The function postCP obtains estimates of posterior probabilities of change-points and hidden states for each observation.
}
\usage{
postcp(formula, data, bp=integer(), family = gaussian, sigma=1.0, maxFB = FALSE)

}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are as given in glm.}
\item{data}{A matrix of data to be segmented. Must have no missing values}
\item{bp}{break points. i.e. For K segments, a vector of K-1 initial estimates of change-point locations, corresponding to the index of the last observation of the first K-1 segments.}
\item{family}{Emission distribution of the observed data as in glm. Supports Poisson(default), Gaussian, Negative Binomial and Gamma distributions.}
\item{sigma}{ Value added to sigma for gaussian data}
\item{maxFB}{ If maxFB = TRUE : call maxFwBk if the most probable change point configuration is required otherwise marginal distributions for change points are returned}
}


\value{
\item{model}{Emission distribution (Poisson, Normal, Gamma, Binomial).}
\item{n}{Length of segmented data}
\item{loglik}{loglikelihood matrix}
\item{post}{Posterior probability distribution}
\item{post.cp}{Posterior change point probability distribution}
\item{le.updated}{Updated Log Evidence matrix}
\item{param.before}{Initial regression coefficients}
\item{param.updated}{Updated regression coefficients using postCP}
\item{response.variable}{Response variable of the data input}
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

The user can input a data matrix and specify the formula required. The package will estimate the initial regression coefficients through the change point model. It will then calculate the posterior change point probabilities and will update the parameter estimates and the log evidence.
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
#prepare data
sigma=1.3
#Change point estimates
bp=c(7,10)
#Obtain data from longley dataset
data = longley
#Apply postcp function
res = postcp(Employed ~ GNP + Armed.Forces,family=gaussian(),data=data,bp=c(7,10),sigma)

#Plot the results
plot.postCP(res,main="Posterior Change Point Probability Distribution")

#Apply postcp function with maxFB=TRUE to obtain marginal distribution
res = postcp(Employed ~ GNP + Armed.Forces,family=gaussian(),data=data,bp=c(7,10),sigma,maxFB=TRUE)

#Plot the results
plot.postCP(res,main="Posterior Change Point Probability Distribution")

}


