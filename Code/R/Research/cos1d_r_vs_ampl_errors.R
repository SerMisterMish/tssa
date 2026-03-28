source("Source/TSSA.R")
library(future.apply)
library(progressr)

progressr::handlers(global = TRUE)
progressr::handlers(progressr::handler_progress)
plan(multisession)

create_ts <- function (ampl, N, freqs, C = 1) {
  C + rowSums(ampl * cos(2 * pi * (1:N) %o% freqs))
}

calc_rec_mse <- function(ts, noise_m, noise_sd, L, r) {
  N <- length(ts)
  noise <- rnorm(N, noise_m, noise_sd)
  rec <- tens_ssa_reconstruct(
    ts + noise,
    L,
    list(1:r),
    svd.method = 'primme',
    decomp = "HOSVD",
    status = FALSE
  )$F1

  if (is.complex(ts)){
    return(mse(ts, rec))
  } else return(mse(ts, Re(rec)))
}

mse <- function(true, pred) mean((true - pred)^2)

rep_calc_mses <- function(repeats,
                         cos_num,
                         ampl,
                         L,
                         N,
                         freqs = NULL,
                         noise_m = 0,
                         noise_sd = 1,
                         C = 1,
                         seed = TRUE) {
  if (is.null(freqs))
    freqs <- seq(from = 0.2,
                 by = 0.02,
                 length.out = cos_num)
  r <- 2 * cos_num + 1
  ts <- create_ts(ampl, N, freqs, C)

  prog <- progressr::progressor(steps = repeats)
  future_replicate(repeats, {
    res <- calc_rec_mse(ts, noise_m, noise_sd, L, r)
    prog()
    res
  }, future.seed = seed)
}

rep_calc_mean_mse <- function(repeats,
                              cos_num,
                              ampl,
                              L,
                              N,
                              freqs = NULL,
                              noise_m = 0,
                              noise_sd = 1,
                              C = 1,
                              seed = TRUE) {
  mean(rep_calc_mses(repeats, cos_num, ampl, L, N, freqs, noise_m, noise_sd, C, seed))
}

repeats <- 500
N <- 499

{
  L <- c(167, 167)
  sapply(0:4,
         \(cn) rep_calc_mean_mse(
           repeats = repeats,
           cos_num = cn,
           ampl = 1,
           L = L,
           N = N,
           seed = 5
         )) |> print()
  # 0: 0.00622592 1: 0.01863985 2: 0.03106062 3: 0.04292674 4: 0.05424432

  rep_calc_mean_mse(
    repeats = repeats,
    cos_num = 1,
    ampl = 5,
    L = L,
    N = N,
    seed = 5
  ) |> print()
  # 0.01833671
}

{
  L <- c(10, 245)
  sapply(0:4,
         \(cn) rep_calc_mean_mse(
           repeats = repeats,
           cos_num = cn,
           ampl = 1,
           L = L,
           N = N,
           seed = 5
         )) |> print()
  # 0: 0.005921844 1: 0.017787906 2: 0.030039175 3: 0.041912390 4: 0.053216666

  rep_calc_mean_mse(
    repeats = repeats,
    cos_num = 1,
    ampl = 5,
    L = L,
    N = N,
    seed = 5
  ) |> print()
  # 0.0176269
}
