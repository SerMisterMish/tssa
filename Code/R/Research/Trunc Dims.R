source("./Research/TSSA.R")

N <- 71
signal <- 30 * cos(2 * pi * 1:N / 12)
r <- 2
groups <- list(1:r)

# I <= J <= L
len.out <- 6

tens.dims <- matrix(ncol = 3, nrow = len.out)
colnames(tens.dims) <- c("I", "L", "J")
tens.dims[, "I"] <- c(12, 12, 7, 12, 19, 24)
tens.dims[, "L"] <- c(49, 37, 36, 31, 30, 25)
tens.dims[, "J"] <- N - tens.dims[, "I"] - tens.dims[, "L"] + 2

trunc_dims <- list(3, 2:3, 1:3)

sigma <- 5
R <- 500
rec.rmse <- list(ssa = numeric(length(Ls)), tssa = matrix(numeric(length(trunc_dims) * len.out), ncol = 3))
mse.tssa <- rep(list(matrix(numeric(R * len.out), ncol = len.out)), length(trunc_dims))

set.seed(5)
pb <- txtProgressBar(max = R * len.out, style = 3)
for (i in seq(R)) {
  noise <- rnorm(N, mean = 0, sd = sigma)
  x <- signal + noise
  for (j in seq(len.out)) {
    I <- tens.dims[j, 1]
    L <- tens.dims[j, 2]

    for (d in seq_along(trunc_dims)) {
      rec.tssa <- tens_ssa_reconstruct(
        x,
        I,
        L,
        groups = groups,
        decomp = "HOOI",
        trunc_dims = trunc_dims[[d]],
        status = FALSE
      )
      
      mse.tssa[[d]][i, j] <- mean(abs(Reduce(`+`, rec.tssa) - signal) ^ 2)
    }
    setTxtProgressBar(pb, (i - 1) * len.out + j)
  }
}

for (d in seq_along(trunc_dims)) {
  rec.rmse$tssa[, d] <- mse.tssa[[d]] |> colMeans() |> sqrt()
}

library(dplyr)
tssa.df <- cbind(as.data.frame(tens.dims), as.data.frame(rec.rmse$tssa))
names(tssa.df) <- c("I", "L", "J", "RMSE_3", "RMSE_2:3", "RMSE_1:3")
tssa.table <- tssa.df |> mutate(IxL = paste(I, L, sep = "x")) |> select(-I, -L, -J)
xtable(t(tssa.table))

# mse.ssa <- matrix(numeric(R * length(Ls)), ncol = length(Ls))
# 
# set.seed(5)
# pb <- txtProgressBar(max = R * length(Ls), style = 3)
# for (i in seq(R))
# {  
#   noise <- rnorm(N, mean = 0, sd = sigma)
#   x <- signal + noise
#   for (j in seq_along(Ls)) {
#     L <- Ls[j]
#     x.ssa <- ssa(x, L = L)
#     
#     rec.ssa <- reconstruct(x.ssa, groups = groups)
#     mse.ssa[i, j] <- mean(abs(Reduce(`+`, rec.ssa) - signal) ^ 2)
#     
#     setTxtProgressBar(pb, (i - 1) * len.out + j)
#   }
# }
# close(pb)
# rec.rmse$ssa <- mse.ssa |> colMeans() |> sqrt()
# ssa.df <- data.frame(L = Ls, RMSE = rec.rmse$ssa)
