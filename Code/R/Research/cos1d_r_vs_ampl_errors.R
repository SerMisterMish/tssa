source("Source/TSSA.R")
library(future.apply)
library(progressr)

progressr::handlers(global = TRUE)
progressr::handlers(progressr::handler_progress)
plan(multisession)

create_ts <- function(ampl, N, freqs, C = 1) {
  C + rowSums(ampl * cos(2 * pi * (1:N) %o% freqs))
}

calc_rec_mse <- function(ts, noise_m, noise_sd, L, r, method=c("tssa", "ssa")) {
  method <- match.arg(method)
  N <- length(ts)
  noise <- rnorm(N, noise_m, noise_sd)
  if (identical(method, "tssa")) {
    rec <- tens_ssa_reconstruct(
      ts + noise,
      L,
      list(1:r),
      svd.method = 'primme',
      decomp = "HOSVD",
      status = FALSE
    )$F1
  } else if (identical(method, "ssa")) {
    rec <- reconstruct(ssa(
      ts + noise,
      L),
      groups = list(1:r),
    )$F1
  }

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
                         method = c("tssa", "ssa"),
                         seed = TRUE) {
  method <- match.arg(method)
  if (is.null(freqs))
    freqs <- seq(from = 0.2,
                 by = 0.02,
                 length.out = cos_num)
  r <- 2 * cos_num + 1
  ts <- create_ts(ampl, N, freqs, C)

  prog <- progressr::progressor(steps = repeats)
  future_replicate(repeats, {
    res <- calc_rec_mse(ts, noise_m, noise_sd, L, r, method)
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
                              method = c("tssa", "ssa"),
                              seed = TRUE) {
  method <- match.arg(method)
  mean(rep_calc_mses(repeats, cos_num, ampl, L, N, freqs, noise_m, noise_sd, C, method, seed))
}

repeats <- 500
N <- 499

# Approx equal dimensions
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

# One small dimension
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

# No signal
{
  L <- c(10, 245)
  sapply(0:4,
         \(cn) rep_calc_mean_mse(
           repeats = repeats,
           cos_num = cn,
           ampl = 0,
           L = L,
           N = N,
           C = 0,
           seed = 5
         )) |> print()
  # 0: 0.001515748, 1: 0.017708967, 2: 0.043529113, 3: 0.074340952, 4: 0.109641057  
}

# W/ signal, ssa
{
  L <- 250
  sapply(0:4,
         \(cn) rep_calc_mean_mse(
           repeats = repeats,
           cos_num = cn,
           ampl = 1,
           L = L,
           N = N,
           seed = 5,
           method = "ssa"
         )) |> print()
  # 0: 0.005974265, 1: 0.017909527, 2: 0.030242741, 3: 0.042285118, 4: 0.053698167
}


# No signal, ssa
{
  L <- 250
  sapply(0:4,
         \(cn) rep_calc_mean_mse(
           repeats = repeats,
           cos_num = cn,
           ampl = 0,
           L = L,
           N = N,
           C = 0,
           seed = 5,
           method = "ssa"
         )) |> print()
  # 0: 0.009972133, 1: 0.043918724, 2: 0.073297362, 3: 0.100151884, 4: 0.125720295 
}
