#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double weighted_sum_kernel_cpp(double hb,
                               NumericVector S,
                               NumericVector pi1,
                               NumericVector pi0) {
  int n = S.size();
  double auc_est = 0.0;
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double diff = (S[i] - S[j]) / hb;
      auc_est +=pi1[i]* pi0[j] * R::pnorm(diff, 0.0, 1.0, 1, 0);
    }
  }
  
  return auc_est;
}

// [[Rcpp::export]]
NumericVector weighted_IF_cpp(double hb,
                               NumericVector S,
                               NumericVector pi1,
                               NumericVector pi0) {
  int n = S.size();
  NumericVector IF_est(n);
  double diff=0.0;
  double pnorm_val=0.0;
  double full_sum=0.0;

  for (int i = 0; i < n; ++i) {
    full_sum=0.0;
    for (int j = 0; j < n; ++j) {
      diff = (S[i] - S[j]) / hb;
      pnorm_val=R::pnorm(diff, 0.0, 1.0, 1, 0);
      full_sum += (pi1[i]* pi0[j] *pnorm_val + pi0[i]*pi1[j]*(1-pnorm_val) );
    }
    IF_est[i]=full_sum/n;
  }
  
  return IF_est;
}


