polynom_eval <- function(coefs, x) {
  powers <- outer(x, seq_along(coefs) - 1, FUN = "^")
  result <- powers %*% coefs
  
  return(as.vector(result))
}

# If parameters are matrices, then a multivariate series will be constructed.
# In this parameters matrices columns correspond to summands and rows to channels
# This means, if you want to increase number of ts in a sum, you grow columns,
# and if you want to add a channel to an mv ts, you grow rows
# The result is either a 1d ts, or a matrix, where each column is a 1d ts -- channel.
construct_ts <- function(N,
                         ampl,
                         rate,
                         freq,
                         phase,
                         complex = TRUE,
                         poly_ampl_coefs = NULL,
                         from = 0,
                         by = 1) {
  if (any(c(length(ampl), length(rate), length(freq)) != length(phase))) {
    stop("Amplitudes, rates, frequencies and phases all must have the same length")
  }
  if (all(is.matrix(ampl), is.matrix(rate), is.matrix(freq), is.matrix(phase))) {
    if (any(c(dim(ampl), dim(rate), dim(freq)) != dim(phase))) {
      stop("Amplitudes, rates, frequencies and phases all must have the same dimensions")
    }
    if (is.null(poly_ampl_coefs)) {
      poly_ampl_coefs <- 1 %o% rep(1, nrow(ampl)) %o% rep(1, ncol(ampl))
    } else {
      stopifnot(
        is.array(poly_ampl_coefs) &&
          length(dim(poly_ampl_coefs)) - 1 == length(dim(ampl)) &&
          all(dim(poly_ampl_coefs)[-1] == dim(ampl))
      )
    }
  } else {
    if (any(is.matrix(ampl), is.matrix(rate), is.matrix(freq), is.matrix(phase))) {
      stop("All parameters should be either matrices, or vectors at the same time")
    }
    
    if (is.null(poly_ampl_coefs)) {
      poly_ampl_coefs <- 1 %o% rep(1, length(ampl))
    } else {
      stopifnot(
        is.matrix(poly_ampl_coefs) &&
          length(dim(poly_ampl_coefs) == 2) &&
          dim(poly_ampl_coefs)[2] == length(ampl)
      )
    }
  }
  
  if (complex) {
    period_f <- exp
    unit <- 1i
  } else {
    period_f <- cos
    unit <- 1
  }
  N_ones <- rep(1, N)
  ts_range <- seq(from = from, length.out = N, by = by)
  poly_ampl <- apply(poly_ampl_coefs, 2:length(dim(poly_ampl_coefs)), polynom_eval, x = ts_range)
  ts_arr <- poly_ampl * (N_ones %o% ampl) * exp(ts_range %o% rate) *
    period_f(N_ones %o% phase + 2 * pi * unit * ts_range %o% freq)
  
  apply(ts_arr, 1:(length(dim(ts_arr)) - 1), sum)
}

calc_rank <- function(ampl, rate, freq, phase, complex = TRUE, poly_ampl_coefs = NULL, tol = 1e-6) {
  if (any(c(length(ampl), length(rate), length(freq)) != length(phase))) {
    stop("Amplitudes, rates, frequencies and phases all must have the same length")
  }
  if (is.null(poly_ampl_coefs)) {
    pac_rank <- rep(1, length(ampl))
  } else {
    pac_rank <- apply(poly_ampl_coefs, 2:length(dim(poly_ampl_coefs)), function(v) max(which(v != 0)))
  }
  
  fr_df <- data.frame(freq = as.vector(freq), rate = as.vector(rate), pac = as.vector(pac_rank))
  fr_unique <- fr_df |>
    group_by(freq, rate) |>
    summarise(pac = max(pac))
  (2 - complex) *
    (sum(ampl != 0) + sum(fr_unique$pac) - length(ampl)) - sum(abs(fr_unique[, 1] - 0.5) < tol)
}


# TODO: This is prone to computational errors, ideally should be like `calc_rank`
calc_rank3 <- function(ts, ...) {
  if (!require("Matrix"))
    install.packages("Matrix")
  
  rankMatrix(ts)[1]
}