library(rTensor)
library(purrr)
library(Rssa)
library(Metrics)

tssa3 <- function(s, I = (length(s) + 2) %/% 3, L = I) {
  X <- tens3(s, I, L)
  result <- list()
  result$X <- X
  result$modes <- list(I = I, L = L, J = length(s) - I - L + 2)
  result$hosvd <- rTensor::hosvd(X)
  return(result)
}

tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N - I - L + 2
  X <- array(NA, c(I, L, J))
  for (i in 1:J) {
    X[, , i] <- outer(1:I, 1:L, function(x, y) s[i + x + y - 2])
  }
  return(as.tensor(X))
}

# diagonal averaging i + j + k = const, const = 3:(I+L+J)
reconstruct.group3 <- function(X.tens) {
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

make.group <- function(hosvd, group) {
  ttl(hosvd$Z[group, group, group, drop = FALSE], list(
    as.matrix(hosvd$U[[1]][, group]),
    as.matrix(hosvd$U[[2]][, group]),
    as.matrix(hosvd$U[[3]][, group])),
      1:3)
}

make.groupHOOI <- function(X, group) {
  r <- length(group)
  hooi_x <- hooi(X@data, r = c(r, r, r), itermax = 500, tol = 1.e-7)
  G <- hooi_x$G
  U <- hooi_x$U
  ## Reconstruct the hooi approximation.
  X_approx <- atrans(G, U)
  as.tensor(X_approx)
}

t3.reconstruct <- function(p, groups) {
  stopifnot(is.list(groups))
  lapply(lapply(groups, make.group, hosvd = p$hosvd), reconstruct.group3)
}

t3.reconstructHOOI <- function(p, groups) {
  require("tensr")
  stopifnot(is.list(groups))
  reconstruct.group3(make.groupHOOI(p$X, group = groups[[1]]))
}

# White gaussian noise

N <- 71

s <- 30 * cos(2 * pi * 1:N / 12)

reps <- 500

res.ssa <- list(numeric(reps), numeric(reps), numeric(reps), numeric(reps))
names(res.ssa) <- c(12, 24, 30, 36)
res.tssa.hosvd <- list(numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps))
res.tssa.hooi <- list(numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps))
names(res.tssa.hosvd) <- c("12x12", "12x24", "12x30", "24x24", "24x30", "30x36")
names(res.tssa.hooi) <- c("12x12", "12x24", "12x30", "24x24", "24x30", "30x36")
pairs <- list(c(12, 12), c(12, 24), c(12, 30), c(24, 24), c(24, 30), c(30, 36))

set.seed(5)
sigma <- 5

for (k in 1:reps) {
  s.noise <- rnorm(N, 0, sigma)
  for (L in c(12, 24, 30, 36)) {
    t.s <- ssa(s + s.noise, L = L)
    t.r <- reconstruct(t.s, list(1:2))
    res.ssa[[as.character(L)]][k] <- mse(t.r[[1]], s)
  }

  for (IL in pairs) {
    capture.output({ t.t <- tssa3((s + s.noise), IL[1], IL[2]) })
    t.tr.hosvd <- t3.reconstruct(t.t, list(1:2))
    t.tr.hooi <- t3.reconstructHOOI(t.t, list(1:2))
    res.tssa.hosvd[[paste(as.character(IL[1]), as.character(IL[2]), sep = "x")]][k] <- mse(t.tr.hosvd[[1]], s)
    res.tssa.hooi[[paste(as.character(IL[1]), as.character(IL[2]), sep = "x")]][k] <- mse(t.tr.hooi, s)
  }
}


me.ssa <- lapply(lapply(lapply(res.ssa, mean), sqrt), round, digits = 4)
me.tssa.hosvd <- lapply(lapply(lapply(res.tssa.hosvd, mean), sqrt), round, digits = 4)
me.tssa.hooi <- lapply(lapply(lapply(res.tssa.hooi, mean), sqrt), round, digits = 4)

cat("RMSE for SSA:\n", unlist(me.ssa), "\n")
cat("RMSE for HOSVD SSA:\n", unlist(me.tssa.hosvd), "\n")
cat("RMSE for HOOI:\n", unlist(me.tssa.hooi), "\n")

# Red gaussian noise

N <- 71

s <- 30 * cos(2 * pi * 1:N / 12)

reps <- 500

res.ssa <- list(numeric(reps), numeric(reps), numeric(reps), numeric(reps))
names(res.ssa) <- c(12, 24, 30, 36)
res.tssa.hosvd <- list(numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps))
res.tssa.hooi <- list(numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps), numeric(reps))
names(res.tssa.hosvd) <- c("12x12", "12x24", "12x30", "24x24", "24x30", "30x36")
names(res.tssa.hooi) <- c("12x12", "12x24", "12x30", "24x24", "24x30", "30x36")
pairs <- list(c(12, 12), c(12, 24), c(12, 30), c(24, 24), c(24, 30), c(30, 36))

set.seed(5)

delta <- sqrt(5)
phi <- 0.5
# phi <- 0.9
for (k in 1:reps) {
  s.noise <- arima.sim(n = N, model = list(ar = phi), sd = delta)
  for (L in c(12, 24, 30, 36)) {
    t.s <- ssa(s + s.noise, L = L)
    t.r <- reconstruct(t.s, list(1:2))
    res.ssa[[as.character(L)]][k] <- mse(t.r[[1]], s)
  }

  for (IL in pairs) {
    capture.output({ t.t <- tssa3((s + s.noise), IL[1], IL[2]) })
    t.tr.hosvd <- t3.reconstruct(t.t, list(1:2))
    t.tr.hooi <- t3.reconstructHOOI(t.t, list(1:2))
    res.tssa.hosvd[[paste(as.character(IL[1]), as.character(IL[2]), sep = "x")]][k] <- mse(t.tr.hosvd[[1]], s)
    res.tssa.hooi[[paste(as.character(IL[1]), as.character(IL[2]), sep = "x")]][k] <- mse(t.tr.hooi, s)
  }
}

me.ssa <- lapply(lapply(lapply(res.ssa, mean), sqrt), round, digits = 4)
me.tssa.hosvd <- lapply(lapply(lapply(res.tssa.hosvd, mean), sqrt), round, digits = 4)
me.tssa.hooi <- lapply(lapply(lapply(res.tssa.hooi, mean), sqrt), round, digits = 4)

cat("RMSE for SSA:\n", unlist(me.ssa), "\n")
cat("RMSE for HOSVD SSA:\n", unlist(me.tssa.hosvd), "\n")
cat("RMSE for HOOI:\n", unlist(me.tssa.hooi), "\n")