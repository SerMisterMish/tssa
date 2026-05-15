setwd(rprojroot::find_rstudio_root_file())

library(abind)
library(Rssa)
library(svd)
library(scales)
library(rTensor)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(quantmod)
library(Rcpp)
library(RcppArmadillo)
library(future.apply)
library(progressr)

source(file = "./Source/TSSA.R")
source(file = "./Source/tsgen.R")
source(file = "./Source/utils.R")
source(file = "./Research/master/pred_funcs.R")

progressr::handlers(global = TRUE)
progressr::handlers("progress")

get_true_params <- function() {
  if (complex.signal) {
    true_params <- c(sort(freqs.true), sort(rates.true))
  } else {
    half_idx <- which(freqs.true %in% c(0, 1 / 2))
    if (length(half_idx) > 0) {
      true_params <- c(sort(c(
        freqs.true[half_idx], rep(freqs.true[-half_idx], each = 2) * c(1, -1)
      )), sort(c(rates.true[half_idx], rep(
        rates.true[-half_idx], 2
      ))))
    } else {
      true_params <- c(sort(rep(freqs.true, each = 2) * c(1, -1)), sort(rep(rates.true, 2)))
    }
  }
  true_params
}

get_param_grids <- function(groups_sep = NULL, est_dim = NULL, trunc_dims = NULL) {
  .groups_ssa <- list(list(1:max.r))
  .Ls_ssa <- seq(from = max.r + 1,
                to = (N + 1) %/% 2,
                by = 5)
  
  .groups_hossa <- list(list(1:max.r))
  if (!is.null(groups_sep)) {
    .groups_hossa_sep <- list(groups_sep)
    .groups_ssa <- list(groups_sep)
  } else {
    .groups_hossa_sep <- list(.groups_hossa)
  }
  
  .td_hossa <- list(NULL)
  
  Is_tmp <- seq(from = max.r + 1,
                to = N - max.r - 1,
                by = 5)
  Ls_tmp <- seq(from = max.r + 1,
                to = N - max.r - 1,
                by = 5)
  Ks_tmp <- seq(from = max.r + 1,
                to = N - max.r - 1,
                by = 5)
  
  correct.dims <- which(
    outer(Is_tmp, 
          Ls_tmp, 
          (\(i, l)ifelse((i >= l) & (N - i - l + 2 > max.r) & (N - i - l + 2 <= l), TRUE, FALSE))
    ), 
    arr.ind = TRUE)
  
  len.out <- nrow(correct.dims)
  tens.dims <- matrix(ncol = 3, nrow = len.out)
  colnames(tens.dims) <- c("I", "L", "J")
  tens.dims[, "I"] <- Is_tmp[correct.dims[, 1]]
  tens.dims[, "L"] <- Ls_tmp[correct.dims[, 2]]
  tens.dims[, "J"] <- N - tens.dims[, "I"] - tens.dims[, "L"] + 2
  
  tens.dims <- tens.dims[order(tens.dims[, 1], tens.dims[, 2], decreasing = TRUE), ]
  .tens.dims_list <- apply(tens.dims, 1, \(v) v[1:2], simplify = list)
  
  # TODO: Could come up with something better
  .tens.dims_list_4d <- list()
  for (i in Is_tmp) {
    for (l in Ls_tmp) {
      if (l > i) next
      for (k in Ks_tmp) {
        j <- N - i - l - k + 3
        if (k > l || j > k || j <= max.r) next
        .tens.dims_list_4d <- c(.tens.dims_list_4d, list(c(i, l, k)))
      }
    }
  }
  
  if (is.null(trunc_dims)) .td_hossa <- list(NULL)
  else .td_hossa <- trunc_dims
  
  .sds <- 0.6
  
  params_grid_ssa <<- list(sd = .sds, L_ssa = .Ls_ssa, groups_ssa = .groups_ssa)
  
  params_grid_hossa <<- list(sd = .sds,
                             L_hossa = .tens.dims_list,
                             use.hutd = FALSE,
                             groups_hossa = .groups_hossa,
                             groups_hossa_sep = .groups_hossa_sep,
                             td_hossa = .td_hossa)
  if (!is.null(est_dim)) params_grid_hossa <<- c(params_grid_hossa, list(est_dim_hossa = est_dim))
  
  params_grid_utssa <<- list(sd = .sds,
                            L_hossa = .tens.dims_list,
                            use.hutd = TRUE,
                            groups_hossa = .groups_hossa,
                            groups_hossa_sep = .groups_hossa_sep,
                            td_hossa = .td_hossa)
  
  params_grid_utssa_4d <<- list(sd = .sds,
                               use.hutd = TRUE,
                               L_hossa = .tens.dims_list_4d,
                               groups_hossa = .groups_hossa,
                               groups_hossa_sep = .groups_hossa_sep)
}

get_signal <- function() {
  max.r <<- calc_rank(A, rates.true, freqs.true, phases.true, complex.signal, poly.true)
  signal_sep <<- construct_ts(N, A, rates.true, freqs.true, phases.true, complex.signal, poly.true, sep = TRUE)
  signal <<- apply(signal_sep, 1, sum)
}

run_num <- function(groups_sep = NULL, est_dim = NULL, rec = FALSE, sep = FALSE, sep_nest = FALSE, est = FALSE, trunc_dims = NULL) {
  snapshot_path_ssa      <- paste0("Research/master/snapshots/", num, "_ssa.RData")
  snapshot_path_hossa    <- paste0("Research/master/snapshots/", num, "_hossa.RData")
  snapshot_path_utssa    <- paste0("Research/master/snapshots/", num, "_utssa.RData")
  snapshot_path_utssa_4d <- paste0("Research/master/snapshots/", num, "_utssa_4d.RData")
  
  get_signal()
  true_all <- list(rec = signal, sep = signal_sep, est = get_true_params())
  get_param_grids(groups_sep = groups_sep, est_dim = est_dim, trunc_dims = trunc_dims)
  
  print(num)
  if (.RECALCULATE || !file.exists(snapshot_path_ssa)) {
    print("ssa")
    ssa_res <- calc_errors_by_params(
      params_grid_ssa,
      R,
      \() pred_ssa_several(rec = rec, sep = sep, est = est),
      rmse_several,
      true_all,
      seed = 5,
      save_snap_to = snapshot_path_ssa
    )
  }
  
  if (.RECALCULATE || !file.exists(snapshot_path_hossa)) {
    print("hossa")
    hossa_res <- calc_errors_by_params(
      params_grid_hossa,
      R,
      \() pred_hossa_several(rec = rec, sep = sep, est = est),
      rmse_several,
      true_all,
      seed = 5,
      save_snap_to = snapshot_path_hossa
    )
  }
  
  if (.RECALCULATE || !file.exists(snapshot_path_utssa)) {
    print("utssa")
    utssa_res <- calc_errors_by_params(
      params_grid_utssa,
      R,
      \() pred_hossa_several(
        rec = rec,
        sep = sep,
        sep_nest = sep_nest
      ),
      rmse_several,
      true_all,
      seed = 5,
      save_snap_to = snapshot_path_utssa
    )
  }
  
  if (.RECALCULATE || !file.exists(snapshot_path_utssa_4d)) {
    print("4d utssa")
    utssa4d_res <- calc_errors_by_params(
      params_grid_utssa_4d,
      R,
      \() pred_hossa_several(rec = rec, sep = sep),
      rmse_several,
      true_all,
      seed = 5,
      save_snap_to = snapshot_path_utssa_4d
    )
  }
}

if (!exists(".RECALCULATE")) .RECALCULATE <- FALSE

R <- 500

# 1. Single channel, $N = 99$
# 1.1 Real case
# 1.1.1 s_n = 1
num <- "1_1_1"
S.true <- 1; rates.true <- 0; freqs.true <- 0; A <- 1
poly.true <- 1 %o% 1; phases.true <- 0; complex.signal <- FALSE
N <- 99
run_num(rec = TRUE)


# 1.1.2 s_n = \exp(n / 50)
num <- "1_1_2"
S.true <- 1; rates.true <- 1/50; freqs.true <- 0; A <- 1
poly.true <- 1 %o% 1; phases.true <- 0; complex.signal <- FALSE
N <- 99
run_num(rec = TRUE, est = TRUE, est_dim = as.list(1:3))

# 1.1.3 s_n = \cos(2 \pi n / 10)
num <- "1_1_3"
S.true <- 1; rates.true <- 0; freqs.true <- 1 / 10; A <- 1
poly.true <- 1 %o% 1; phases.true <- 0; complex.signal <- FALSE
N <- 99
run_num(rec = TRUE, est = TRUE, est_dim = as.list(1:3))

# 1.1.4 s_n = \exp(n / 50) \cos(2 \pi n / 10)
num <- "1_1_4"
S.true <- 1; rates.true <- 1/50; freqs.true <- 1 / 10; A <- 1
poly.true <- 1 %o% 1; phases.true <- 0; complex.signal <- FALSE
N <- 99
run_num(rec = TRUE, est = TRUE, est_dim = as.list(1:3))

# 1.1.5 s_n = 2 + \cos(2 \pi n / 10)
num <- "1_1_5"
S.true <- 2; rates.true <- c(0, 0); freqs.true <- c(0, 1 / 10); A <- c(2, 1)
poly.true <- 1 %o% c(1, 1); phases.true <- c(0, 0); complex.signal <- FALSE
N <- 99
run_num(rec = TRUE, sep = TRUE, groups_sep = list(1, 2:3))

# 1.1.6 s_n = 2 \cos(2 \pi n / 5) + \cos(2 \pi n / 10)
num <- "1_1_6"
S.true <- 2; rates.true <- c(0, 0); freqs.true <- c(1 / 5, 1 / 10); A <- c(2, 1)
poly.true <- 1 %o% c(1, 1); phases.true <- c(0, 0); complex.signal <- FALSE
N <- 99
run_num(rec = TRUE, sep = TRUE, groups_sep = list(1:2, 3:4), est = TRUE, est_dim = as.list(1:3))

# 1.1.7 s_n = 1 + n / 10
num <- "1_1_7"
S.true <- 1; rates.true <- 0; freqs.true <- 0; A <- 1
poly.true <- c(1, 1/10) %o% 1; phases.true <- 0; complex.signal <- FALSE
N <- 99
run_num(rec = TRUE)

# 1.2.1 s_n = \exp(2 \pi \iu n / 10)
num <- "1_2_1"
S.true <- 1; rates.true <- 0; freqs.true <- 1/10; A <- 1
poly.true <- 1 %o% 1; phases.true <- 0; complex.signal <- TRUE
N <- 99
run_num(rec = TRUE, est = TRUE, est_dim = as.list(1:3))

# 1.2.2 s_n = \exp(n / 50 + 2 \pi \iu n / 10)
num <- "1_2_2"
S.true <- 1; rates.true <- 1/50; freqs.true <- 1/10; A <- 1
poly.true <- 1 %o% 1; phases.true <- 0; complex.signal <- TRUE
N <- 99
run_num(rec = TRUE, est = TRUE, est_dim = as.list(1:3))

# 1.2.3 s_n = \exp(n / 50 + 2 \pi \iu n / 10) + 1
num <- "1_2_3"
S.true <- 2; rates.true <- c(1/50, 0); freqs.true <- c(1/10, 0); A <- c(1, 1)
poly.true <- 1 %o% c(1, 1); phases.true <- c(0, 0); complex.signal <- TRUE
N <- 99
run_num(rec = TRUE, sep = TRUE, groups_sep = list(1, 2))

# 1.2.4 s_n = 2 \exp(2 \pi \iu n / 5) + \exp(2 \pi \iu n / 5)
num <- "1_2_4"
S.true <- 2; rates.true <- c(0, 0); freqs.true <- c(1/5, 1/10); A <- c(2, 1)
poly.true <- 1 %o% c(1, 1); phases.true <- c(0, 0); complex.signal <- TRUE
N <- 99
run_num(rec = TRUE, est = TRUE, est_dim = as.list(1:3), sep = TRUE, groups_sep = list(1, 2))

# 1.2.4 (with truncation dimensions)
num <- "1_2_4_td"
S.true <- 2; rates.true <- c(0, 0); freqs.true <- c(1/5, 1/10); A <- c(2, 1)
poly.true <- 1 %o% c(1, 1); phases.true <- c(0, 0); complex.signal <- TRUE
N <- 99
run_num(rec = TRUE, trunc_dims = list(1, 2, 3, 1:2, 2:3))
