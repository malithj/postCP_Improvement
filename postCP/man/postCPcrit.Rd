\name{postCPcrit}
\alias{postCPcrit}
\title{
  Estimate model selection criteria in change-point model.
}
\description{
  Obtains model selection criteria for change-point problems, after an initial set of change-points within the data has already been obtained. Calculates ICL (Information Completed Likelihood), AIC (Akaike's Information Criterion), BIC (Bayesian Information Criterion), and mBIC (modified Bayesian Information Criterion). Calculated for segment-based model, with different levels estimated for each segment.
}
\usage{
postCPcrit(data,seg=numeric(),model,disp.equal=TRUE,
eps.nb=1e-8,param.nb=NULL,prior=0.5,prior.type="n") 
}
\arguments{
\item{data}{A vector of observations to be segmented. Must have no missing values.}
\item{seg}{For K segments, a vector of K-1 initial estimates of change-point locations, corresponding to the index of the last observation of the first K-1 segments.}
\item{model}{Emission distribution of observed data, 1: Poisson (default), 2: normal, 3: negative binomial.}
\item{disp.equal}{Equal dispersion parameter across all segments for negative binomial (defaults to TRUE). If disp.equal=TRUE, seg will automatically be chosen by the PDPA algorithm from Segmentor function.}
\item{eps.nb}{Epsilon for maximum likelihood estimation of negative binomial, if mle does not converge, set to higher value.}
\item{param.nb}{Parameters for maximum likelihood estimation of negative binomial if already fit by postCPmodelsel.}
\item{prior}{Transition probabilities chosen as priors for HMM model. For homogeneous HMM, enter a scalar between 0 and 1 (default: 0.5). For heterogeneous HMM, enter a vector of probabilities, of length n for n observations or K for K hidden states.}
\item{prior.type}{Type of prior entered. For "n":homogeneous HMM (default), "o": vector of probabilities for each observation, "s": vector of probabilities for each hidden state.}
}

\value{
\item{scores}{Vector of criteria.}
\item{cp.loc}{A posteriori most probable set of K-1 change-points, by Viterbi algorithm.}
}

\details{
BIC: Bayesian Information Criterion, expressed as -LL+p*log(n), where LL is the log-likelihood, p is the number of parameters, and n is the length of observations
AIC: Akaike Information Criterion, expressed as -LL+2*p, where LL is the log-likelihood, p is the number of parameters.
ICL: Information Completed Likelihood, includes entropy term to take into account clustering of observations within segments, obtained using posterior probability of observations from postCP algorithm.
mBIC: For normal and Poisson data, modified Bayesian Information Criterion, according to Zhang and Siegmund's estimate.
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
  Luong, T.M., Rozenholc, Y. & Nuel, G. (2012). Fast estimation of posterior probabilities in change-point models through a constrained hidden Markov model. \url{http://arxiv.org/pdf/1203.4394}\cr
  Zhang, N.R. & Siegmund, D.O. (2007). A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data. Biometrics, 63(1):22--32.\cr
  Picard, F., Robin, S., Lavielle, M., Vaisse, C., & Daudin, J.J. (2005) A statistical approach for array CGH data analysis. BMC bioinformatics, 6(1):27.\cr
  Schwarz, G.E. (1978). Estimating the dimension of a model. Annals of Statistics 6 (2): 461-464
}


\keyword{ package }
\seealso{
\code{\link{postCP}}
\code{\link{postCPmodelsel}}
}

\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
# obtain ICL, AIC, BIC, mBIC assuming Poisson model
postCPcrit(x,c(N,2*N),model=1);
}
