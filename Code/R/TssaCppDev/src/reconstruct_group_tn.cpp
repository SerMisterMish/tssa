#include <RcppArmadillo/Lighter>

static inline size_t flatten_arr_idx(const Rcpp::IntegerVector& dims, const std::vector<size_t>& idx) {
  size_t n = dims.size(), lin = 0, stride = 1;
  for (size_t i = 0; i < n; ++i) {
    size_t k = idx[i];
    if (k < 0 || k >= static_cast<size_t>(dims[i])) 
      throw std::out_of_range(std::format("Index `{}` out of bounds for dimension {}", k + 1, i + 1));
    lin += k * stride;
    stride *= dims[i];
  }
  return lin;
}

static inline Rcomplex ndarr_get_elem(const Rcpp::ComplexVector& data, const std::vector<size_t>& idx) {
  Rcpp::IntegerVector dims = data.attr("dim");
  size_t lin = flatten_arr_idx(dims, idx);
  return data[lin];
}

static inline bool idx_inc(std::vector<size_t>& idx, const Rcpp::IntegerVector& upper_bounds) {
  size_t n = idx.size();
  if (static_cast<size_t>(upper_bounds.size()) != n)
    throw std::length_error("Idx and Upper Bounds lengths not equal");

  bool result = false;
  for (size_t k = 0; k < n; ++k) {
    size_t i = n - k - 1;
    size_t ub = static_cast<size_t>(upper_bounds[i]);
    
    if (ub <= 1) 
      throw std::out_of_range(std::format(
          "Upper Bounds must be greater than 1. Got `{}` for dimension `{}`", 
          ub + 1, 
          i + 1));
      
    idx[i] += 1;
    if (idx[i] < ub) {
      result = true;
      break;
    }
    
    idx[i] = 0;
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::ComplexVector reconstruct_group_tn(const Rcpp::ComplexVector& data) {
  Rcpp::IntegerVector dims = data.attr("dim");
  size_t N = 0, n = dims.size();
  for (size_t i = 0; i < n; ++i) 
    N += dims[i] - 1;
  ++N;

  Rcpp::ComplexVector result(N, 0);
  std::vector<size_t> idx(n, 0), counts(N, 0);
  
  do {
    size_t idx_flat = 0;
    for (size_t i = 0; i < n; ++i) 
      idx_flat += idx[i];
    
    result[idx_flat] = result[idx_flat] + ndarr_get_elem(data, idx);
    ++counts[idx_flat];
  } while (idx_inc(idx, dims));
    
  for (size_t i = 0; i < N; ++i) {
    result[i].r /= counts[i];
    result[i].i /= counts[i];
  }
  
  return result;
}