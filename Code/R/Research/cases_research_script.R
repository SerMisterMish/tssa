setwd(rprojroot::find_rstudio_root_file())
library(Rssa)
library(svd)
library(scales)
library(rTensor)
library(ggplot2)
library(dplyr)
library(tidyr)
library(english)
library(Rcpp)
library(RcppArmadillo)
library(future.apply)
library(progressr)

source(file = "./Source/TSSA.R")
source(file = "./Source/tsgen.R")
source(file = "./Source/utils.R")

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

Ps <- c(2, 4, 8, 12, 16, 20, 40)

for (P in Ps) {
  P_case <- sprintf("_c%d", P)
  print(P_case)
  set.seed(5)
  
  # S.true <- 1
  S.true <- 2
  S_case <- sprintf("%s_periodic%s", as.english(S.true), ifelse(S.true > 1, "s", ""))

  # rates.true <- matrix(rep(0, P * S.true), nrow = P)
  # rates_case <- "_nr"
  # rates.true <- matrix(rep(-0.04, P * S.true), nrow = P)
  # rates_case <- "_Ler"
  rates.true <- matrix(rep(-0.04, P * S.true) + rep(seq(
    from = 0,
    length.out = S.true,
    by = -0.02
  ), each = P), nrow = P)
  rates_case <- "_dsr"

  # freqs.true <- matrix(rep(0.2, P * S.true), nrow = P)
  # freqs_case <- "_ef"
  freqs.true <- rep(1/8, P) %o% seq(from = 1, to = 8 / 6, length.out = S.true)
  freqs_case <- "_dsf"

  A <- matrix(rnorm(P * S.true), nrow = P)
  A_case <- "_da"

  # poly.true <- 1 %o% rep(1, P) %o% rep(1, S.true)
  # poly_case <- ""
  poly.true <- c(1, 1) * matrix(runif(2 * P, 0.5, 2), nrow = 2) %o% rnorm(S.true, mean = 1, sd = 0.25)
  poly_case <- sprintf("p%d", dim(poly.true)[1] - 1)

  phases.true <- matrix(rep(0, P * S.true), nrow = P)
  phases_case <- "_ep"

  # complex.signal <- TRUE
  # complex_case <- "_c" # complex
  complex.signal <- FALSE
  complex_case <- "_r" # real

  # N <- 25
  # n_case <- ""
  N <- 99
  n_case <- "_long"

  full_case <- paste0(n_case, P_case, rates_case, freqs_case, A_case, poly_case, phases_case, complex_case)
  max.r <- calc_rank(A, rates.true, freqs.true, phases.true, complex.signal, poly.true)

  stopifnot(N > max.r * 2 + 1)
  signal <- construct_ts(N, A, rates.true, freqs.true, phases.true, complex.signal, poly.true)

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
        by = 5
      )
      Ls_homssa <- seq(
        from = max.r + 1,
        to = (N + 1) %/% 2,
        by = 5
      )
    }

    max.r3 <- calc_rank3(signal)
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
    save(ssa_err, tssa_err, file = paste0("./Research/comp_results/cases_errors/", S_case, "_mv", full_case, ".RData"))
  } else {
    ssa_err <- calc_errors_by_params(params_grid_ssa, R, pred_ssa_rec, rmse_ts_nd, signal, seed = 5)
    tssa_err <- calc_errors_by_params(params_grid_hossa, R, pred_hossa_rec, rmse_ts_nd, signal, seed = 5)
    save(ssa_err, tssa_err, file = paste0("./Research/comp_results/cases_errors/", S_case, full_case, ".RData"))
  }
}
