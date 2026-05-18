#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector der_h_weighting_XI(NumericVector pi1,
                                    NumericVector pi0,
                                   NumericVector S,
                                   NumericMatrix x,
                                   double hb) {
  int n = S.size();
  int p = x.ncol();
  double w = 0;
  double diff = 0;
  double Phi = 0;
  NumericVector res(p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      w = pi1[i] * pi0[j];
      diff = (S[i] - S[j]) / hb;
      Phi = R::pnorm(diff, 0.0, 1.0, 1, 0);
      for (int k = 0; k < p; k++) {
        res[k] += w * Phi * (x(i, k) + x(j, k));
      }
    }
  }
  return res/(n*n);
}

// [[Rcpp::export]]
NumericVector der_h_weighting_beta(NumericVector pi1,
                                    NumericVector pi0,
                                    NumericVector S,
                                    NumericVector py,
                                    NumericMatrix W,
                                    double hb) {
  int n = S.size();
  int pq = W.ncol();
  NumericVector res(pq);
  double w = 0;
  double diff = 0;
  double Phi = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      w = pi1[i] * pi0[j];
      diff = (S[i] - S[j]) / hb;
      Phi = R::pnorm(diff, 0.0, 1.0, 1, 0);
      for (int k = 0; k < pq; k++) {
        res[k] += w * Phi *( (1-py[i])*W(i, k) -py[j]* W(j, k));
      }
    }
  }
  return res/(n*n);
}
// [[Rcpp::export]]
NumericVector der_h_model_beta(NumericVector pi1,
                                   NumericVector pi0,
                                   NumericVector S,
                                   NumericVector py,
                                   NumericMatrix Predictor,
                                   double hb) {
  int n = S.size();
  int pq = Predictor.ncol();
  NumericVector res(pq);
  double w = 0;
  double diff = 0;
  double dnorm_phi = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      w = pi1[i] * pi0[j];
      diff = (S[i] - S[j]) / hb;
      dnorm_phi = R::dnorm(diff, 0.0, 1.0,  0)/hb;
      for (int k = 0; k < pq; k++) {
        res[k] += w * dnorm_phi *( Predictor(i, k) - Predictor(j, k));
      }
    }
  }
  return res/(n*n);
}