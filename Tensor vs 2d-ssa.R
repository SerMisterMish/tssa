library(rTensor)
library(purrr)
library(Rssa)

set.seed(1)

N <- 71

signal1 <- 30 * cos(2 * pi * 1:N / 12)
# signal2 <- 20 * cos(2 * pi * 1:N / 12) # case 1
# r <- 2
# r3 <- 1
# signal2 <- 20 * cos(2 * pi * 1:N / 12 + pi / 4) # case 2
# r <- 2
# r3 <- 2
signal2 <- 20 * cos(2 * pi * 1:N / 8 + pi / 4) # case 3
r <- 4
r3 <- 2

series.number <- 2

signal <- cbind(signal1, signal2)
R <- 100

sigma <- 5
Ls <- c(12, 24, 36, 48, 60)
# Ls <- seq(36, 50, 2)

mssa.errors <- function(Ls) {
  noise1 <- rnorm(N, sd = sigma)
  # noise1 <- arima.sim(n = N, model = list(ar = 0.5), sd = sqrt(5))
  # noise1 <- arima.sim(n = N, model = list(ar = 0.9), sd = sqrt(5))
  noise2 <- rnorm(N, sd = sigma)
  # noise2 <- arima.sim(n = N, model = list(ar = 0.5), sd = sqrt(5))
  # noise2 <- arima.sim(n = N, model = list(ar = 0.9), sd = sqrt(5))
  f1 <- signal1 + noise1
  f2 <- signal2 + noise2
  f <- cbind(f1, f2)
  f.list <- list(f1, f2)
  err.rec <- list(mat = numeric(length(Ls)), tens = numeric(length(Ls)))
  names(err.rec$mat) <- Ls; names(err.rec$tens) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    # s <- ssa(f, L = L, kind = "mssa")
    s <- ssa(f, L = c(L, series.number / 2), kind = "2d-ssa")
    # rec <- reduce(reconstruct(s, groups = list(1:r))[[1]], cbind)
    rec <- reconstruct(s, groups = list(1:r))[[1]]
    err.rec$mat[l] <- mean((rec - signal)^2)

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

cat("RMSE for 2D-SSA:\n", sqrt(rowMeans(simplify2array(mres[1,]))), "\n")
cat("RMSE for HOSVD MSSA:\n", sqrt(rowMeans(simplify2array(mres[2,]))), "\n")

# Значимость
# t.test(simplify2array(mres[1,])[6,], simplify2array(mres[2,])[5,], paired = TRUE, alternative = "greater")