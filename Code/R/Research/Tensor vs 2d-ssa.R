library(rTensor)
library(purrr)
library(Rssa)

set.seed(1)

N <- 71

signal1 <- 30 * cos(2 * pi * 1:N / 12)

signal2 <- 20 * cos(2 * pi * 1:N / 12)
signal3 <- 25 * cos(2 * pi * 1:N / 12)
signal4 <- 35 * cos(2 * pi * 1:N / 12)
r <- 2
r3 <- 1
rd <- 4
# signal2 <- 20 * cos(2 * pi * 1:N / 12 + pi / 4)
# signal3 <- 25 * cos(2 * pi * 1:N / 12)
# signal4 <- 35 * cos(2 * pi * 1:N / 12)
# r <- 2
# r3 <- 2
# rd <- 4
# signal2 <- 20 * cos(2 * pi * 1:N / 8 + pi / 4)
# signal3 <- 25 * cos(2 * pi * 1:N / 12)
# signal4 <- 35 * cos(2 * pi * 1:N / 12)
# r <- 4
# r3 <- 2
# rd <- 8
# signal2 <- 20 * cos(2 * pi * 1:N / 12 + pi / 4)
# signal3 <- 25 * cos(2 * pi * 1:N / 12 + pi / 2)
# signal4 <- 35 * cos(2 * pi * 1:N / 12 + 3 * pi / 4)
# r <- 2
# r3 <- 2
# rd <- 4
# signal2 <- 20 * cos(2 * pi * 1:N / 8 + pi / 4)
# signal3 <- 25 * cos(2 * pi * 1:N / 12 + pi / 2)
# signal4 <- 35 * cos(2 * pi * 1:N / 8 + 3 * pi / 4)
# r <- 4
# r3 <- 4
# rd <- 8
# signal2 <- 20 * cos(2 * pi * 1:N / 8 + pi / 4)
# signal3 <- 25 * cos(2 * pi * 1:N / 10 + pi / 2)
# signal4 <- 35 * cos(2 * pi * 1:N / 14 + 3 * pi / 4)
# r <- 8
# r3 <- 4
# rd <- 12

signal <- cbind(signal1, signal2, signal3, signal4)

series.number <- 4

R <- 100

sigma <- 5
# Ls <- c(12, 24, 36, 48, 60)
Ls <- seq(36, 50, 2)

mssa.errors <- function(Ls) {
  noise1 <- rnorm(N, sd = sigma)
  noise2 <- rnorm(N, sd = sigma)
  noise3 <- rnorm(N, sd = sigma)
  noise4 <- rnorm(N, sd = sigma)
  f1 <- signal1 + noise1
  f2 <- signal2 + noise2
  f3 <- signal3 + noise3
  f4 <- signal4 + noise4
  f <- cbind(f1, f2, f3, f4)
  f.list <- list(f1, f2, f3, f4)
  err.rec <- list(mat = numeric(length(Ls)), mat_2d = numeric(length(Ls)), tens = numeric(length(Ls)))
  names(err.rec$mat) <- Ls; names(err.rec$tens) <- Ls; names(err.rec$mat_2d) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s1 <- ssa(f.list, L = L, kind = "mssa")
    s2 <- ssa(f, L = c(L, series.number / 2), kind = "2d-ssa")
    rec1 <- reduce(reconstruct(s1, groups = list(1:r))[[1]], cbind)
    rec2 <- reconstruct(s2, groups = list(1:rd))[[1]]
    err.rec$mat[l] <- mean((rec1 - signal)^2)
    err.rec$mat_2d[l] <- mean((rec2 - signal)^2)

    f.hankel <- lapply(f.list, hankel, L)
    mat.f <- reduce(f.hankel, cbind)
    s.tens <- fold(mat.f, 1, 2:3, modes = c(L, N - L + 1, series.number))
    capture.output({ h <- rTensor::hosvd(s.tens) })
    tens <- ttl(h$Z[1:r, 1:r, 1:r3, drop = FALSE], list(h$U[[1]][, 1:r, drop = FALSE], h$U[[2]][, 1:r, drop = FALSE], h$U[[3]][, 1:r3, drop = FALSE]), 1:3)
    rec.tens <- reduce(apply(tens@data, 3, hankel, simplify = FALSE), cbind)
    err.rec$tens[l] <- mean((rec.tens - signal)^2)
  }
  err.rec
}

mres <- replicate(R, mssa.errors(Ls))

cat("RMSE for MSSA:\n", sqrt(rowMeans(simplify2array(mres[1,]))), "\n")
cat("RMSE for 2D-SSA:\n", sqrt(rowMeans(simplify2array(mres[2,]))), "\n")
cat("RMSE for HOSVD MSSA:\n", sqrt(rowMeans(simplify2array(mres[3,]))), "\n")

# Значимость
# t.test(simplify2array(mres[1,])[6,], simplify2array(mres[2,])[5,], paired = TRUE, alternative = "greater")