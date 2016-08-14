#include <Rcpp.h>
using namespace Rcpp;

// This is the core implementation of the forward-backward model
//
//
//
//
inline double logsumexp(double loga,double logb) {
  if (loga>logb)
    return loga+log1p(exp(logb-loga));
  else
    return logb+log1p(exp(loga-logb));
}

// [[Rcpp::export]]
SEXP FwBk(Rcpp::NumericMatrix ledata, Rcpp::NumericMatrix workdata){
  Rcpp::NumericMatrix le(ledata);
  Rcpp::NumericMatrix work(workdata);
  int n=le.nrow();
  int K=le.ncol();
  // set constraints
  for (int j=1; j<K; ++j) le(0,j)=-1e300;
  for (int j=0; j<K-1; ++j) le(n-1,j)=-1e300;
  // backward in work
  for (int j=0; j<K; ++j) work(n-1,j)=0.;
  for (int i=n-1; i>0; --i) {
    work(i-1,K-1)=le(i,K-1)+work(i,K-1);
    for (int j=0; j<K-1; ++j) work(i-1,j)=logsumexp(le(i,j)+work(i,j),le(i,j+1)+work(i,j+1));
  }
  // forward in le
  for (int i=1; i<n; ++i) {
    le(i,0)+=le(i-1,0);
    for (int j=1; j<K; ++j) le(i,j)+=logsumexp(le(i-1,j-1),le(i-1,j));
  }
  return(wrap(le(n-1,K-1)));
}

inline double max(double loga,double logb) {
  if (loga>logb)
    return loga;
  else
    return logb;
}

// [[Rcpp::export]]
SEXP maxFwBk(Rcpp::NumericMatrix ledata, Rcpp::NumericMatrix workdata){
  Rcpp::NumericMatrix le(ledata);
  Rcpp::NumericMatrix work(workdata);
  int n=le.nrow();
  int K=le.ncol();
  // set constraints
  for (int j=1; j<K; ++j) le(0,j)=-1e300;
  for (int j=0; j<K-1; ++j) le(n-1,j)=-1e300;
  // backward in work
  for (int j=0; j<K; ++j) work(n-1,j)=0.;
  for (int i=n-1; i>0; --i) {
    work(i-1,K-1)=le(i,K-1)+work(i,K-1);
    for (int j=0; j<K-1; ++j) work(i-1,j)=max(le(i,j)+work(i,j),le(i,j+1)+work(i,j+1));
  }
  // forward in le
  for (int i=1; i<n; ++i) {
    le(i,0)+=le(i-1,0);
    for (int j=1; j<K; ++j) le(i,j)+=max(le(i-1,j-1),le(i-1,j));
  }
  return(wrap(le(n-1,K-1)));
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
