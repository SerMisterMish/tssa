library(Rssa)
library(rTensor)
library(ggplot2)
library(dplyr)
library(tidyr)

source(file = "./Source/TSSA.R")
source(file = "./Source/tsgen.R")
source(file = "./Source/utils.R")

S.true <- 2
A <- c(1, 2)
rates.true <- rep(0, S.true)
freqs.true <- c(0.2, 0.1)
poly.true <- 1 %o% rep(1, S.true)
phases.true <- rep(0, S.true)
complex.signal <- FALSE
# N <- 198
N <- 199

signal <- construct_ts(N, A, rates.true, freqs.true, phases.true, complex.signal, poly.true)

L <- 100
K <- N - L + 1
ss <- ssa(signal, L)
ssu <- ss$U[, 1:4]
ssv <- calc.v(ss, idx=1:4)
sssigm <- diag(ss$sigma[1:4])
sssigm

ssu_lag <- rbind(ssu[-1,], c(ssu[L - 1 / freqs.true[2] + 1, 1:2], ssu[L - 1 / freqs.true[1] + 1, 3:4]))
ssv_lag <- rbind(ssv[-1,], c(ssv[K - 1 / freqs.true[2] + 1, 1:2], ssv[K - 1 / freqs.true[1] + 1, 3:4]))

ssu_blag <- rbind(c(ssu[1 / freqs.true[2], 1:2], ssu[1 / freqs.true[1], 3:4]), ssu[-L,])
ssv_blag <- rbind(c(ssv[1 / freqs.true[2], 1:2], ssv[1 / freqs.true[1], 3:4]), ssv[-K,])
ssm <- as.matrix(Rssa:::.get.or.create.trajmat.1d.ssa(ss))

sssigm_lag <- crossprod(ssu_lag, ssm) %*% ssv_lag
sssigm_lag |> zapsmall()

sssigm_blag <- crossprod(ssu_blag, ssm) %*% ssv_blag
sssigm_blag |> zapsmall()

sst1 <- tens_ssa_decompose(signal, c(100, 50), delete.repeated = FALSE)
sst2 <- tens_ssa_decompose(signal, c(100, 50), delete.repeated = TRUE)
sst1$Z@data[1:4, 1:4, 1:4] |> round(2)
sst2$Z@data[1:4, 1:4, 1:4] |> round(2)
