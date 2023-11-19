library(rTensor)
library(purrr)
library(Rssa)

set.seed(1)

N <- 71

signal1 <- 30 * cos(2 * pi * 1:N / 12)
signal2 <- 20 * cos(2 * pi * 1:N / 12) # case 1
# signal2 <- 20 * cos(2 * pi * 1:N / 12 + pi / 4) # case 2
# signal2 <- 20 * cos(2 * pi * 1:N / 8 + pi / 4) # case 3

series.number <- 2

# series rank
# cases 1 & 2
r <- 2
# case 3
# r <- 4
# rank along the third dimension
# case 1
r3 <- 1
# cases 2 & 3
# r3 <- 2

signal <- cbind(signal1, signal2)
R <- 500

sigma <- 5
Ls <- c(12, 24, 36, 48, 60)

mssa.errors <- function(Ls) {
  f1 <- signal1 + rnorm(N, sd = sigma)
  f2 <- signal2 + rnorm(N, sd = sigma)
  f <- list(f1, f2)
  err.rec <- list(mat = numeric(length(Ls)), tens = numeric(length(Ls)))
  names(err.rec$mat) <- Ls; names(err.rec$tens) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s <- ssa(f, L = L, kind = "mssa")
    rec <- reduce(reconstruct(s, groups = list(1:r))[[1]], cbind)
    err.rec$mat[l] <- mean((rec - signal)^2)

    f.hankel <- lapply(f, hankel, L)
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
cat("RMSE for HOSVD MSSA:\n", sqrt(rowMeans(simplify2array(mres[2,]))), "\n")