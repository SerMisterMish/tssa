#include <RcppArmadillo/Lighter>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector reconstruct_group_t3(const arma::dcube& data) {
  auto I = data.n_rows, L = data.n_cols, K = data.n_slices, N = I + L + K - 2;
  int count;
  double sum;
  Rcpp::NumericVector result(N);
  for (int C = 2; C < N + 2; ++C) {
    sum = 0; count = 0;
    for (int i = 0; i < C - 1; ++i) {
      for (int l = 0; l < C - i; ++l) {
        if (i < I && l < L && C - i - l - 2 < K) {
          sum += data(i, l, C - i - l - 2);
          ++count;
        }
      }
    }
    result[C - 2] = sum / count;
  }
  return result;
}
