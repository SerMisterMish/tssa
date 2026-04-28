#include <RcppArmadillo/Lighter>

// [[Rcpp::export]]
Rcpp::ComplexVector reconstruct_group_tn(const Rcpp::ComplexVector& data) {
  Rcpp::IntegerVector dims = data.attr("dim");
  size_t N = 0, n = dims.size();
  for (size_t i = 0; i < n; ++i) N += dims[i] - 1;
  ++N;

  for (auto it = data.begin(); it != data.end(); ++it) {
    std::cout << *it << ' ';
  }
  return Rcpp::ComplexVector(0);
  /*
  size_t count;
  Rcpp::Rcomplex sum, mean;
  Rcpp::ComplexVector result(N);
  
  for (size_t C = n - 1; C < N + n - 1; ++C) {
    sum.r = 0; sum.i = 0;
    count = 0;
    std::vector<size_t> idx(n, 0);
    
    
    for (size_t i = 0; i < C - 1; ++i) {
      for (size_t l = 0; l < C - i; ++l) {
        if (i < dims[0] && l < dims[1] && C - i - l - 2 < dims[2]) {
          sum += data(i, l, C - i - l - 2);
          ++count;
        }
      }
    }
    
    mean = sum / static_cast<double>(count);
    result[C - 2] = mean; 
  }
  return result;
   */
}