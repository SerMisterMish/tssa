library(tensr)
library(rTensor)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
library(reshape2)
library(Rssa)
library(parallel)

cl <- makeCluster(detectCores()-1)

clusterEvalQ(cl, library(tensr))
clusterEvalQ(cl, library(rTensor))
clusterEvalQ(cl, library(purrr))
clusterEvalQ(cl, library(ggplot2))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tidyr))
clusterEvalQ(cl, library(xtable))
clusterEvalQ(cl, library(reshape2))
clusterEvalQ(cl, library(Rssa))

tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N - I - L + 2
  X <- outer(1:L, 1:J, function(l, j) l + j) |> outer(1:I, function(lj, i) s[lj + i - 2])
  return(as.tensor(X))
}

# diagonal averaging i + j + k = const, const = 3:(I+L+J)
reconstruct.group3 <- function(X.tens) {
  stopifnot(is(X.tens, "Tensor"))
  X <- X.tens@data
  I <- length(X[, 1, 1])
  L <- length(X[1, , 1])
  J <- length(X[1, 1,])
  s <- vector(mode = "numeric", length = I + L + J - 2)
  for (C in 3:(I + L + J)) {
    sum <- 0
    count <- 0
    for (i in 1:(C - 2)) {
      for (l in 1:(C - 1 - i)) {
        if (i <= I && l <= L && C - i - l <= J) {
          sum <- sum + X[i, l, C - i - l]
          count <- count + 1
        }
      }
    }
    s[C - 2] <- sum / count
  }
  return(s)
}

make.groupHOOI <- function(X, r) {
  hooi_x <- hooi(X@data, r = c(r, r, r))
  G <- hooi_x$G
  U <- hooi_x$U
  ## Reconstruct the hooi approximation.
  X_approx <- atrans(G, U)
  as.tensor(X_approx)
}

t3.reconstructHOOI <- function(X, r) {
  reconstruct.group3(make.groupHOOI(X, r))
}

# set.seed(1)
clusterSetRNGStream(cl, 1)

N <- 71
signal.l <- list()

signal.l[[1]] <- 30 * cos(2 * pi * 1:N / 12)

# Two signals
# series.number <- 2

# signal.l[[2]] <- 20 * cos(2 * pi * 1:N / 12)
# r <- 2
# r.ssa <- 2
# r3 <- 1
# rd <- 2
# plot_name <- "two-series-first"

# signal.l[[2]] <- 20 * cos(2 * pi * 1:N / 12 + pi / 4)
# r <- 2
# r.ssa <- 2
# r3 <- 2
# rd <- 2
# plot_name <- "two-series-second"
#
# signal.l[[2]] <- 20 * cos(2 * pi * 1:N / 8 + pi / 4)
# r <- 4
# r.ssa <- 2
# r3 <- 2
# rd <- 4
# plot_name <- "two-series-third"

# Five signals
# series.number <- 5
#
# signal.l[[2]] <- 30 * cos(2 * pi * 1:N / 12)
# signal.l[[3]] <- 30 * cos(2 * pi * 1:N / 12)
# signal.l[[4]] <- 30 * cos(2 * pi * 1:N / 12)
# signal.l[[5]] <- 30 * cos(2 * pi * 1:N / 12)
# r <- 2
# r.ssa <- 2
# r3 <- 1
# rd <- 2
# plot_name <- "five-series-first"
#
# signal.l[[2]] <- 20 * cos(2 * pi * 1:N / 12)
# signal.l[[3]] <- 25 * cos(2 * pi * 1:N / 12)
# signal.l[[4]] <- 30 * cos(2 * pi * 1:N / 12)
# signal.l[[5]] <- 20 * cos(2 * pi * 1:N / 12)
# r <- 2
# r.ssa <- 2
# r3 <- 1
# rd <- 6
# plot_name <- "five-series-second"
#
# signal.l[[2]] <- 30 * cos(2 * pi * 1:N / 12 + 2 * pi / 3)
# signal.l[[3]] <- 30 * cos(2 * pi * 1:N / 12 + 4 * pi / 3)
# signal.l[[4]] <- 30 * cos(2 * pi * 1:N / 12)
# signal.l[[5]] <- 30 * cos(2 * pi * 1:N / 12 + 2 * pi / 3)
# r <- 2
# r.ssa <- 2
# r3 <- 2
# rd <- 2
# plot_name <- "five-series-third"
#
# signal.l[[2]] <- 30 * cos(2 * pi * 1:N / 12 + 2 * pi / 3)
# signal.l[[3]] <- 30 * cos(2 * pi * 1:N / 12 + 4 * pi / 3)
# signal.l[[4]] <- 30 * cos(2 * pi * 1:N / 8)
# signal.l[[5]] <- 30 * cos(2 * pi * 1:N / 8 + 2 * pi / 3)
# r <- 4
# r.ssa <- 2
# r3 <- 4
# rd <- 10
# plot_name <- "five-series-fourth"
#
# Nine signals
series.number <- 9

# for (i in 2:series.number) {
#   signal.l[[i]] <- 30 * cos(2 * pi * 1:N / 12)
# }
# r <- 2
# r.ssa <- 2
# r3 <- 1
# rd <- 2
# plot_name <- "nine-series-first"
#
# ampl <- c(30, 28, 26, 24, 22)
# for (i in 2:series.number) {
#   signal.l[[i]] <- ampl[(i - 1) %% 5 + 1] * cos(2 * pi * 1:N / 12)
# }
# r <- 2
# r.ssa <- 2
# r3 <- 1
# rd <- 10
# plot_name <- "nine-series-second"
#
# for (i in 2:series.number) {
#   signal.l[[i]] <- 30 * cos(2 * pi * 1:N / 12 + 2 * (i - 1) * pi / 5)
# }
# r <- 2
# r.ssa <- 2
# r3 <- 2
# rd <- 2
# plot_name <- "nine-series-third"
#
for (i in 2:(series.number-2)) {
  signal.l[[i]] <- 30 * cos(2 * pi * 1:N / 12 + 2 * (i - 1) * pi / 5)
}
signal.l[[series.number - 1]] <- 30 * cos(2 * pi * 1:N / 8 + 2 * (series.number - 2) * pi / 5)
signal.l[[series.number]] <- 30 * cos(2 * pi * 1:N / 8 + 2 * (series.number - 1) * pi / 5)
r <- 4
r.ssa <- 2
r3 <- 4
rd <- 10
plot_name <- "nine-series-fourth"

signal <- reduce(signal.l, cbind)

R <- 500
# R <- 1

sigma <- 5
# Ls <- c(12, 24, 36, 48, 60)
Ls <- seq(36, 56, 2)
Is <- seq(12, 30, 6)

mssa.errors <- function(Ls, Is) {
  noise <- replicate(series.number, rnorm(N, sd = sigma), simplify = FALSE)
  # noise <- replicate(series.number,0, simplify = FALSE)
  f.list <- imap(signal.l, function(x, i) x + noise[[i]])
  f <- reduce(f.list, cbind)
  err.rec <- list(ssa = numeric(length(Ls)), hooi_ssa = matrix(numeric(length(Ls)*length(Is)), length(Ls)),
                  mssa = numeric(length(Ls)), hosvd_mssa = numeric(length(Ls)), `2d-ssa` = numeric(length(Ls)))
  err.rec$hooi_ssa[,] <- NA
  names(err.rec$ssa) <- Ls
  colnames(err.rec$hooi_ssa) <- Is
  rownames(err.rec$hooi_ssa) <- Ls
  names(err.rec$mssa) <- Ls
  names(err.rec$hosvd_mssa) <- Ls
  names(err.rec$`2d-ssa`) <- Ls
  for (l in seq_along(Ls)) {
    L <- Ls[l]
    s.ssa <- lapply(f.list, ssa, L = L, kind = "1d-ssa")
    rec.ssa <- reduce(lapply(lapply(s.ssa, reconstruct, groups = list(1:r.ssa)), unlist, use.names = FALSE), cbind)
    err.rec$ssa[l] <- mean((rec.ssa - signal)^2)

    s.mssa <- ssa(f.list, L = L, kind = "mssa")
    rec.mssa <- reduce(reconstruct(s.mssa, groups = list(1:r))[[1]], cbind)
    err.rec$mssa[l] <- mean((rec.mssa - signal)^2)

    s.2d_ssa <- ssa(f, L = c(L, (series.number+1) %/% 2), kind = "2d-ssa")
    rec.2d_ssa <- reconstruct(s.2d_ssa, groups = list(1:rd))[[1]]
    err.rec$`2d-ssa`[l] <- mean((rec.2d_ssa - signal)^2)

    f.hankel <- lapply(f.list, hankel, L)
    mat.f <- reduce(f.hankel, cbind)
    tens <- fold(mat.f, 1, 2:3, modes = c(L, N - L + 1, series.number))
    capture.output({ h <- rTensor::hosvd(tens) })
    s.hosvd.mssa <- ttl(h$Z[1:r, 1:r, 1:r3, drop = FALSE],
                        list(h$U[[1]][, 1:r, drop = FALSE], h$U[[2]][, 1:r, drop = FALSE],
                             h$U[[3]][, 1:r3, drop = FALSE]), 1:3)
    rec.hosvd.mssa <- reduce(apply(s.hosvd.mssa@data, 3, hankel, simplify = FALSE), cbind)
    err.rec$hosvd_mssa[l] <- mean((rec.hosvd.mssa - signal)^2)

    visited <- list()
    for (i in seq_along(Is)) {
      I <- Is[i]
      K <- N - L - I + 2
      if (K >= r.ssa && !(list(c(K, L)) %in% visited || list(c(K, I)) %in% visited ||
        c(L, K) %in% visited || c(I, K) %in% visited)) {
        s.hooi.ssa <- lapply(f.list, tens3, L = L, I = I)
        rec.hooi.ssa <- reduce(lapply(s.hooi.ssa, t3.reconstructHOOI, r = r.ssa), cbind)
        err.rec$hooi_ssa[l,i] <- mean((rec.hooi.ssa - signal)^2)
        visited <- c(visited, list(c(I, L)))
      }
    }
  }
  err.rec
}

clusterExport(cl, c("signal", "signal.l", "sigma", "series.number", "N", "r", "r.ssa", "rd", "r3",
                   "tens3", "reconstruct.group3", "make.groupHOOI", "t3.reconstructHOOI",
                    "Ls", "Is", "mssa.errors"))

system.time( {mres <- parSapply(cl, 1:R, function(i, ...) mssa.errors(Ls, Is)) })
stopCluster(cl)
# mres <- replicate(R, mssa.errors(Ls, Is))

save(mres, file = "nine-series-results-par4.rdata")

rmse.df <- list()
cat("\nRMSE for SSA:\n")
(rmse.df$ssa <- sqrt(rowMeans(simplify2array(mres[1,]))))
cat("\nRMSE for HOOI SSA:\n")
(rmse.df$hooi_ssa <- sqrt(apply(simplify2array(mres[2,]), c(1,2), mean)))
cat("\nRMSE for MSSA:\n")
(rmse.df$mssa <- sqrt(rowMeans(simplify2array(mres[3,]))))
cat("\nRMSE for HOSVD MSSA:\n")
(rmse.df$hosvd_mssa <- sqrt(rowMeans(simplify2array(mres[4,]))))
cat("\nRMSE for 2D-SSA:\n")
(rmse.df$`2d_ssa` <- sqrt(rowMeans(simplify2array(mres[5,]))))

rmse.df |> lapply(min, na.rm = TRUE) |> as.data.frame() |> xtable()

hooi.ssa.rmse <- rmse.df$hooi_ssa
rmse.df$hooi_ssa <- NULL
rmse.df <- rmse.df |> as.data.frame()
names(rmse.df)[4] <- "2d_ssa"

pdf(paste0("./src/img/", plot_name, ".pdf"), width = 8, height = 4)
rmse.df |>
  mutate(L = Ls) |>
  pivot_longer(cols = c(ssa, mssa, hosvd_mssa, `2d_ssa`), names_to = "Method") |>
  ggplot(aes(x = L, y = value)) + geom_line(aes(color = Method, linetype = Method)) +
  ylab("RMSE") + theme_minimal()
dev.off()

pdf(paste0("./src/img/", plot_name, "_hooi.pdf"), width = 4, height = 4)
hooi.ssa.rmse |>
  melt() |> filter(!is.na(value)) |>
  ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill = value)) +
  xlab("I") + ylab("L") + labs(fill = "RMSE") +
  theme_minimal()
dev.off()

# Значимость
# t.test(simplify2array(mres[1,])[6,], simplify2array(mres[2,])[5,], paired = TRUE, alternative = "greater")
