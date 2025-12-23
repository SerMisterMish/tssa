library(rTensor)
library(Rssa)

mse <- function(true, pred) {
  err <- pred - true
  Re(mean(err * Conj(err)))
}

# true is scalar, pred is vector
rmse <- function(true, pred) {
  sqrt(mse(true, pred))
}

# true is vector, pred is matrix with predictions for each parameter in rows
rmse_params_mat <- function(true, pred) {
  apply(cbind(true, pred), 1, \(v) rmse(v[1], v[-1]))
}

# true is nd-array, pred is (n+1)d-array with n+1-mode slices as predictions
rmse_ts_nd <- function(true, pred) {
  sqrt(mean(apply(pred, length(dim(pred)), \(p) mse(true, p))))
}

mixing_rate <- function(x, r, ...)
  UseMethod("mixing_rate")

# 0 - no mixing, 1 and greater - mixed
setMethod("mixing_rate", "numeric", function(x, r, ...) {
  if (length(x) - r > 1)
    (x[r + 1] - x[r + 2]) / (x[r] - x[r + 1])
  else
    x[r + 1] / (x[r] - x[r + 1])
})

mixing_rate.ssa <- function(x, r) {
  mixing_rate(x$sigma, r)
}

setMethod("mixing_rate", "Tensor", function(x, r, dims = seq(length(x@modes))) {
  sigmas <- get_tensor_eigenvalues(x, dims)
  if (length(r) == 1)
    sapply(sigmas, mixing_rate, r = r)
  else
    sapply(seq_along(sigmas), \(i) mixing_rate(sigmas[[i]], r[i]))
})

setMethod("mixing_rate", "list", function(x, r, dims = seq(length(x$Z@modes))) {
  mixing_rate(x$Z, r, dims)
})
