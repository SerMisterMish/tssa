library(Rssa)
library(rTensor)
library(ggplot2)
library(dplyr)
library(tidyr)
library(future.apply)
library(progressr)

setwd(rprojroot::find_rstudio_root_file())
source(file = "./Research/TSSA.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession)

# get_cp_resid_fnorm <- function(R, s, nc) {
#   p <- progressor(R)
#   future_replicate(R, {
#     p()
#     noise <- matrix(rnorm(N * Q, sd = sd), ncol = Q)
#     tens_mssa_decompose(s + noise,
#                         L, 
#                         decomp = "cp",
#                         neig = nc, 
#                         status = FALSE,
#                         max_iter = max_iter,
#                         tol = tol,
#                         start = "svd")$fnorm_resid
#   }, future.seed = 5)
# }

# get_ho_resid_fnorm <- function(R, s, r, r3) {
#   p <- progressor(R)
#   future_replicate(R, {
#     p()
#     noise <- matrix(rnorm(N * Q, sd = sd), ncol = Q)
#     tens_mssa_decompose(s + noise,
#                         L, 
#                         decomp = "hosvd",
#                         neig = r,
#                         neig3 = r3,
#                         status = FALSE)$fnorm_resid
#   }, future.seed = 5)
# }

get_resid_fnorm <- function(R) {
  p <- progressor(R)
  future_replicate(R, {
    p()
    noise <- rnorm(N * Q, sd = sd)
    x <- s + noise
    c(
      noise = noise |> fnorm_complex(),
      mssa = (x - Reduce("+", reconstruct(ssa(x, L_mssa, kind = "mssa"), groups = groups_mssa))) |> fnorm_complex(),
      cp = (x - Reduce(
        "+",
        tens_mssa_reconstruct(
          x,
          L_cp,
          decomp = "cp",
          groups = groups_cp,
          status = FALSE,
          max_iter = max_iter,
          tol = tol,
          start = "svd"
        )
      )) |> fnorm_complex(),
      ho = (x - Reduce(
        "+",
        tens_mssa_reconstruct(
          x,
          L_ho,
          decomp = "hosvd",
          groups = groups_mssa,
          groups3 = groups3,
          status = FALSE,
          max_iter = max_iter,
          tol = tol
        )
      )) |> fnorm_complex()
    )
  }, future.seed = 5)
}

get_bias <- function() {
  c(
    mssa = (s - Reduce("+", reconstruct(ssa(s, L_mssa, kind = "mssa"), groups = groups_mssa))) |> fnorm_complex(),
    cp = (s - Reduce(
      "+",
      tens_mssa_reconstruct(
        s,
        L_cp,
        decomp = "cp",
        groups = groups_cp,
        status = FALSE,
        max_iter = max_iter,
        tol = tol,
        start = "svd"
      )
    )) |> fnorm_complex(),
    ho = (s - Reduce(
      "+",
      tens_mssa_reconstruct(
        s,
        L_ho,
        decomp = "hosvd",
        groups = groups_mssa,
        groups3 = groups3,
        status = FALSE,
        max_iter = max_iter,
        tol = tol
      )
    )) |> fnorm_complex()
  )
}

N <- 71
Q <- 9
# sd = 5
sd = 0.03
R <- 500
ones <- rep(1, N)
A <- 20 * exp(ones %o% c(rep(0.05, 7), 0.03, 0.03) * 1:N)
w <- rep(1 / 12, 9)
phi <- ones %o% rep(0, 9)
max_iter = 25
tol = 1e-7

s <- A * cos(2 * pi * 1:N %o% w + phi)

groups_cp <- list(1:4)
groups_cp <- list(1:8)
# groups_mssa <- list(1:4)
groups3 <- list(1:2)
L_mssa <- 56
L_ho <- 22
L_cp <- 56

resid_fnorms <- get_resid_fnorm(R)
bias_fnorms <- get_bias()
cat(sprintf("Ranks: r_mssa = %i, r_3 = %i, r_cp = %i\n\n", 
            max(sapply(groups_mssa, max)),
            max(sapply(groups3, max)),
            max(sapply(groups_cp, max))),
  "\rResidual norms:\n", 
    sprintf("%s:\t%f\n", rownames(resid_fnorms), rowMeans(resid_fnorms)),
    "\n\rBias norms:\n",
    sprintf("%s:\t%f\n", names(bias_fnorms), bias_fnorms))

# nc <- 5
# r <- 4
# r3 <- 2

# noise_fnorm <- future_replicate(R, fnorm_complex(tens3(matrix(
#   rnorm(N * Q, sd = sd), ncol = Q
# ), L, kind = "MSSA")), future.seed = 5)
# 
# signal_cp_resid <- tens_mssa_decompose(
#   s,
#   L,
#   decomp = "cp",
#   neig = nc,
#   status = FALSE,
#   start = "svd"
# )$fnorm_resid
# signal_ho_resid <- tens_mssa_decompose(
#   s,
#   L,
#   decomp = "hosvd",
#   neig = r,
#   neig3 = r3,
#   status = FALSE
# )$fnorm_resid
# 
# noised_signal_cp_resid <- get_cp_resid_fnorm(R, s, nc)
# noised_signal_ho_resid <- get_ho_resid_fnorm(R, s, r, r3)
# 
# cat(sprintf(
#   " noise: %f\n\nCPD:\nsignal: %f\n   s+n: %f\n\nHOSVD:\nsignal: %f\n   s+n: %f\n",
#   mean(noise_fnorm),
#   mean(signal_cp_resid),
#   mean(noised_signal_cp_resid),
#   mean(signal_ho_resid),
#   mean(noised_signal_ho_resid)
# ))
