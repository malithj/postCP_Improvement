\name{postCPmodelsel}
\alias{postCPmodelsel}
\title{
  Model selection for change-point data.
}
\description{
  Estimates number of segments in data using model selection criteria for change-point problems. Minimizes ICL (Information Completed Likelihood), AIC (Akaike's Information Criterion), BIC (Bayesian Information Criterion), and mBIC (modified Bayesian Information Criterion) for values specified in K.range.
}
\usage{
postCPmodelsel(data,K.range,model=1,greedy=FALSE,
disp.equal=TRUE,eps.nb=1e-8,prior=0.5,prior.type="n") 
}
\arguments{
\item{data}{A vector of observations to be segmented. Must have no missing values.}
\item{K.range}{Vector of range of possible K number of segments.}
\item{model}{Emission distribution of observed data, 1: Poisson (default), 2: normal, 3: negative binomial.}
\item{greedy}{Algorithm used to obtain initial segmentation for each fixed K. TRUE: use greedy algorithm based on least-squares minimization, FALSE (default): use pruned dynamic programming algorithm (PDPA), using Segmentor function.}
\item{disp.equal}{Equal dispersion parameter across all segments for negative binomial (defaults to TRUE).}
\item{eps.nb}{Epsilon for maximum likelihood estimation of negative binomial, if mle does not converge, set to higher value.}
\item{prior}{Transition probabilities chosen as priors for HMM model. For homogeneous HMM, enter a scalar between 0 and 1 (default: 0.5). For heterogeneous HMM, enter a vector of probabilities, of length n for n observations or K for K hidden states.}
\item{prior.type}{Type of prior entered. For "n":homogeneous HMM (default), "o": vector of probabilities for each observation, "s": vector of probabilities for each hidden state.}
}

\value{
\item{ICL}{Estimated number of segments by ICL (Information Completed Likelihood) criterion.}
\item{AIC}{Estimated number of segments by AIC (Akaike Information Criterion).}
\item{BIC}{Estimated number of segments by BIC (Bayesian Information Criterion).}
\item{mBIC}{Estimated number of segments by mBIC (modified Bayesian Information Criterion).}
\item{scores}{Matrix of criteria, rows for each possible $K$ segments, columns for each criteria.}
\item{cp.loc}{A posteriori most probable set of change-points, by Viterbi algorithm, for most probable number of change-points according to each criteria.}
}

\details{
BIC: Bayesian Information Criterion, expressed as -LL+p*log(n), where LL is the log-likelihood, p is the number of parameters, and n is the length of observations
AIC: Akaike Information Criterion, expressed as -LL+2*p, where LL is the log-likelihood, p is the number of parameters.
ICL: Information Completed Likelihood, includes entropy term to take into account clustering of observations within segments, obtained using posterior probability of observations from postCP algorithm.
mBIC: For normal and Poisson data, modified Bayesian Information Criterion, according to Zhang and Siegmund's estimate.

Greedy segmentation consists of greedy algorithm using locally optimal binary segmentations through least squares minimization. Non-greedy segmentation uses pruned dynamic programming algorithm (PDPA). Algorithm assume homogeneous priors for states.
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
  Cleynen, A., Koskas, M., & Rigaill, G. (2012) A Generic Implementation of the Pruned Dynamic Programing Algorithm. \url{http://arxiv.org/pdf/1204.5564}\cr
  Zhang, N.R. & Siegmund, D.O. (2007). A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data. Biometrics, 63(1):22--32.\cr
  Picard, F., Robin, S., Lavielle, M., Vaisse, C., & Daudin, J.J. (2005) A statistical approach for array CGH data analysis. BMC bioinformatics, 6(1):27.\cr
  Schwarz, G.E. (1978). Estimating the dimension of a model. Annals of Statistics 6 (2): 461-464
}

\keyword{ package }
\seealso{
\code{\link{postCP}}
\code{\link{postCPcrit}}
}

\examples{
require(postCP);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
# estimate number of segments by ICL, AIC, BIC, mBIC assuming Poisson model
# choose number of segment (1 to 10) minimizing each criteria
postCPmodelsel(x,1:10,model=1);
}

