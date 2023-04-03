tssa <- function(s, rank, l = floor(sqrt(length(s))), l1 = (l + 1) %/% 2, o = 1, max_iter = 25, tol = 1e-05) {
  require("rTensor")
  X <- tens(s, l, l1, o)
  result <- list()
  result$X <- X
  result$l <- l
  result$l1 <- l1
  result$o <- o
  result$cp <- cp(X, num_components = rank, max_iter = max_iter, tol = tol)
  return(result)
}

tens <- function(s, l, l1, o) {
  v <- as.vector(s)
  N <- length(v)
  I <- N %/% l
  X.cap <- matrix(v, ncol = l, nrow = I, byrow = TRUE)
  J <- ((l - l1) %/% o + 1)
  X <- array(NA, c(J, l1, I))
  for (i in 1:I) {
    X[, , i] <- outer(0:(J - 1), 1:l1, function(x, y) X.cap[i, x * o + y])
  }
  return(as.tensor(X))
}

reconstruct.group <- function(X.tens, l, o) {
  X <- X.tens
  I <- length(X[1, 1,])
  L <- nrow(as.matrix(X[, , 1]))
  K <- ncol(as.matrix(X[, , 1]))
  X.cap <- matrix(NA, nrow = I, ncol = l)
  left <- c(1:(L * o), L * ((o + 1):K))
  right <- c(1 + L * seq.int(0, K - 1, o), rep(K * L, l - K %/% o))
  for (i in 1:I) {
    X.cap[i,] <- sapply(1:l, function(j) mean(X[, , i][seq.int(left[j], right[j], by = o * L - 1)]))
  }
  return(as.vector(t(X.cap)))
}

t.reconstruct <- function(p, groups) {
  stopifnot(is.list(groups))
  lapply(lapply(groups, make.group, cp = p$cp), reconstruct.group, l = p$l, o = p$o)
}


make.group <- function (cp, group) {
  Reduce("+", lapply(group, function (i) cp$lambdas[i] * (cp$U[[1]][,i] %o% cp$U[[2]][,i] %o% cp$U[[3]][,i])))
}

# Test 1: const
s <- rep(3, 20)
p <- tssa(s, 1)
rec <- t.reconstruct(p, list(1))
mse(rec[[1]], s)

# Test 2: sin
s <- sin(2 * pi * 0:125 / 3 + pi / 3)
p <- tssa(s, 5, 9, 6)
rec <- t.reconstruct(p, list(1:2))
mse(rec[[1]], s)

# Test 3: const + sin
s.const <- 3
s.sin <- sin(2 * pi * 0:120 / 3)
s <- s.const + s.sin
p <- tssa(s, 39, 39)
rec <- t.reconstruct(p, list(1, 2:3))
mse(rec[[1]], s.const)
mse(rec[[2]], s.sin)

# Test 4: tssa vs ssa with noise
library(Rssa)
library(Metrics)
M <- 1000
tssa.deviation <- numeric(0)
ssa.deviation <- numeric(0)
N <- 9
s.sin <- sin(2 * pi * 1:N / 3 + pi / 2)# + sin(2 * pi * 1:N / 20 + pi / 3)
r <- 2
pb <- txtProgressBar(max = M, style = 3)
for (i in 1:M) {
  # with this noise tssa performs not worse than ssa
  s.noise <- 0.1 * arima.sim(n=N, model = list(ar = (0.9)))
  s <- s.sin + s.noise
  s.ssa <- ssa(s, L = 4, kind = "1d-ssa")
  capture.output({ s.tssa <- tssa(s, 4, 2) })
  s.trec <- t.reconstruct(s.tssa, groups = list(1:r))
  s.rec <- reconstruct(s.ssa, groups = list(1:r))
  tssa.deviation <- c(tssa.deviation, mse(s.sin, s.trec[[1]]))
  ssa.deviation <- c(ssa.deviation, mse(s.sin, s.rec$F1))
  setTxtProgressBar(pb,i)
}
close(pb)

print(c(" ssa  ", round(sqrt(mean(ssa.deviation)), digits = 4)), quote = FALSE)
print(c(" tssa ", round(sqrt(mean(tssa.deviation)), digits = 4)), quote = FALSE)

