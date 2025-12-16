setwd(rprojroot::find_rstudio_root_file())
library(Rssa)
library(svd)
library(scales)
library(rTensor)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Rcpp)
library(RcppArmadillo)
library(future.apply)
library(progressr)

source(file = "./Source/TSSA.R")
progressr::handlers(global = TRUE)
progressr::handlers("progress")

pred_ssa_rec <- function() {
  noise <- rnorm(N, 0, sd)
  x <- signal + noise
  rec <- reconstruct(ssa(x, L = L_ssa), groups = groups_ssa)
  Reduce("+", rec)
}

pred_hossa_rec <- function() {
  noise <- rnorm(N, 0, sd)
  x <- signal + noise

  rec <- tens_ssa_reconstruct(
    x,
    # c(I_hossa, L_hossa),
    unlist(L_hossa),
    groups = groups_hossa,
    decomp = "hooi",
    trunc_dims = unlist(td_hossa),
    status = FALSE
  )

  Reduce("+", rec)
}

pred_mssa_rec <- function() {
  noise <- rnorm(N * P, 0, sd)
  x <- signal + noise
  # rec <- reconstruct(ssa(x, L = L_mssa, kind = 'mssa'), groups = groups_mssa)
  rec <- cmssa_reconstruct(x, L = L_mssa, groups = groups_mssa)
  Reduce("+", rec)
}

pred_homssa_rec <- function() {
  noise <- rnorm(N * P, 0, sd)
  x <- signal + noise

  rec <- tens_mssa_reconstruct(
    x,
    L_homssa,
    groups = groups_homssa,
    groups3 = groups3_homssa,
    decomp = "hooi",
    status = FALSE
  )

  Reduce("+", rec)
}

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

calc_errors <- function(R, pred_func, err_func, signal, progbar, seed = 5) {
  # with(plan(multicore, workers = availableCores() - 1), {
  with(plan(multisession, workers = availableCores() - 1), {
    batch_recs <- function(R) {
      # p <- progressor(R)
      future_replicate(R,
        {
          # p()
          progbar()
          pred_func()
        },
        simplify = "array",
        future.seed = seed
      )
    }
    if (sd == 0) {
      recs <- batch_recs(2)
    } # don't need R because no noise present. Set to 2, because with 1 will be simplified
    else {
      recs <- batch_recs(R)
    }
  })

  return(err_func(signal, recs))
}

# params_grid_list: list with the parameter names as names and vectors of the parameters values as values
calc_errors_by_params <- function(params_grid_list, R, comp_func, err_func, signal, seed = 5) {
  params_grid <- expand.grid(params_grid_list)
  results <- params_grid
  results$error <- numeric(nrow(params_grid))
  q <- progressor(nrow(params_grid) * R)
  # q <- progressor(nrow(params_grid))
  for (i in seq(nrow(params_grid))) {
    # print(paste0(names(params_grid), " = ", params_grid[i,], collapse = ", "))
    list2env(params_grid[i, ], envir = parent.frame())
    results[i, ]$error <- calc_errors(R, comp_func, err_func, signal, q, seed)
    # q()
  }
  return(results)
}

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
construct_periodics <- function(N,
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
  ts_arr <- poly_ampl * (N_ones %o% ampl) * (N_ones %o% exp(rate)) *
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
  fr_df <- data.frame(freq = freq, rate = rate, pac = pac_rank)
  fr_unique <- fr_df |>
    group_by(freq, rate) |>
    summarise(pac = max(pac))
  (2 - complex) *
    (sum(ampl != 0) + sum(fr_unique$pac) - length(ampl)) - sum(abs(fr_unique[, 1] - 0.5) < tol)
}

print_errors <- function(err_df) {
  err_df |>
    mutate(zapped_error = zapsmall(error)) |>
    group_by(sd) |>
    arrange(zapped_error, L_ssa, .by_group = TRUE)
}

print_errors_m <- function(err_df) {
  err_df |>
    mutate(zapped_error = zapsmall(error)) |>
    group_by(sd) |>
    arrange(zapped_error, L_mssa, .by_group = TRUE)
}

summarise_errors <- function(err_df) {
  err_df |>
    mutate(zapped_error = zapsmall(error)) |>
    group_by(sd) |>
    summarise(
      best = min(zapped_error),
      best_L = L[which.min(zapped_error)],
      avg = mean(zapped_error),
      worst = max(zapped_error)
    )
}

Ps <- c(2, 4, 8, 12, 16, 20)

for (P in Ps) {
  P_case <- sprintf("_c%d", P)
  print(P_case)
  set.seed(5)

  rates.true <- matrix(rep(0, P), nrow = P)
  rates_case <- "_nr"

  freqs.true <- matrix(rep(0.2, P), nrow = P)
  freqs_case <- "_ef"

  A <- matrix(rnorm(P), nrow = P)
  A_case <- "_da"

  poly.true <- 1 %o% rep(1, P) %o% 1
  poly_case <- ""
  # poly.true <- c(1, 1) * matrix(runif(2 * P, 0.5, 2), nrow = 2) %o% 1
  # poly_case <- "p1"

  phases.true <- matrix(rep(0, P), nrow = P)
  phases_case <- "_ep"

  # complex.signal <- TRUE
  # complex_case <- "_c" # complex
  complex.signal <- FALSE
  complex_case <- "_r" # real

  N <- 25
  n_case <- ""
  # N <- 100
  # n_case <- "_long"

  full_case <- paste0(n_case, P_case, rates_case, freqs_case, A_case, poly_case, phases_case, complex_case)
  max.r <- calc_rank(A, rates.true, freqs.true, phases.true, complex.signal, poly.true)

  stopifnot(N > max.r * 2 + 1)
  signal <- construct_periodics(N, A, rates.true, freqs.true, phases.true, complex.signal, poly.true)

  if (is.matrix(signal)) {
    groups_mssa <- list(1:max.r)

    if (n_case == "") {
      Ls_mssa <- seq(
        from = max.r + 1,
        to = N - max.r - 1,
        by = 2
      )
      Ls_homssa <- seq(
        from = max.r + 1,
        to = (N + 1) %/% 2,
        by = 2
      )
    } else {
      Ls_mssa <- seq(
        from = max.r + 1,
        to = N - max.r - 1,
        by = 10
      )
      Ls_homssa <- seq(
        from = max.r + 1,
        to = (N + 1) %/% 2,
        by = 10
      )
    }

    max.r3 <- tens_mssa_decompose(signal, N %/% 2, status = FALSE)$Z@data |>
      zapsmall() |>
      apply(3, \(x) any(x != 0)) |>
      sum()
    groups_homssa <- list(1:max.r)
    groups3_homssa <- list(1:max.r3)

    sds <- seq(0, 0.6, by = 0.2)

    params_grid_ssa <- list(sd = sds, L_mssa = Ls_mssa)
    params_grid_hossa <- list(
      sd = sds,
      L_homssa = Ls_homssa
    )
  }

  R <- 500

  if (is.matrix(signal)) {
    ssa_err <- calc_errors_by_params(params_grid_ssa, R, pred_mssa_rec, rmse_ts_nd, signal, seed = 5)
    tssa_err <- calc_errors_by_params(params_grid_hossa, R, pred_homssa_rec, rmse_ts_nd, signal, seed = 5)
    save(ssa_err, tssa_err, file = paste0("./Research/comp_results/cases_errors/one_periodic_mv", full_case, ".RData"))
  } else {
    ssa_err <- calc_errors_by_params(params_grid_ssa, R, pred_ssa_rec, rmse_ts_nd, signal, seed = 5)
    tssa_err <- calc_errors_by_params(params_grid_hossa, R, pred_hossa_rec, rmse_ts_nd, signal, seed = 5)
    save(ssa_err, tssa_err, file = paste0("./Research/comp_results/cases_errors/one_periodic", full_case, ".RData"))
  }
}
