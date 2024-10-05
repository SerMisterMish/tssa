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

tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N - I - L + 2
  X <- outer(1:L, 1:J, function(l, j) l + j) |> outer(1:I, function(lj, i) s[lj + i - 2])
  return(as.tensor(X))
}

reconstruct.group <- function(X.tens, l, o) {
  stopifnot(is(X.tens, "Tensor"))
  X <- X.tens@data
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

tssa3 <- function(s, rank, I = (length(s) + 2) %/% 3, L = I, max_iter = 100, tol = 1e-05) {
  X <- tens3(s, I, L)
  result <- list()
  result$X <- X
  result$modes <- list(I = I, L = L, J = length(s)-I-L+2)
  result$cp <- cp(X, num_components = rank, max_iter = max_iter, tol=tol)
  return(result)
}

tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N-I-L + 2
  X <- array(NA, c(I, L, J))
  # print(c(I, L, J))
  for (i in 1:J) {
    X[, , i] <- outer(1:I, 1:L, function(x, y) s[i+ x + y - 2])
  }
  return(as.tensor(X))
}

t3.reconstruct <- function(p, groups) {
  stopifnot(is.list(groups))
  lapply(lapply(groups, make.group, cp = p$cp), reconstruct.group3)
}

reconstruct.group3 <- function (X.tens) {
  X <- X.tens
  I <- length(X[,1,1])
  L <- length(X[1,,1])
  J <- length(X[1,1,])
  s <- vector(mode = "numeric", length = I + L + J - 2)
  for (C in 3:(I+L+J)){
    sum <- 0
    count <- 0
    for (i in 1:(C - 2)){
      for (l in 1:(C - 1 - i)){
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

require("Metrics")
# Test 1: const
s <- rep(3, 16)
p <- tssa(s, 1, 6, 6)
p <- tssa3(s, 1, 6, 6)
rec <- t.reconstruct(p, list(1))
mse(rec[[1]], s)

# Test 2: sin
s <- sin(2 * pi * 0:15 / 3)
p <- tssa(s, 2, 9, 6)
while (p$cp$conv != TRUE || p$cp$norm_percent <=95) {
  p <- tssa3(s, 3, 6, 6, max_iter = 10000, tol=1e-8)
}
rec <- t.reconstruct(p, list(1:2))
mse(rec[[1]], s)

# Test 3: const + sin
s.const <- 3
s.sin <- sin(2 * pi * 0:15/ 3)
s <- s.const + s.sin
p <- tssa(s, 3, 11, 6, max_iter = 10000, tol=1e-12)
while (p2$cp$norm_percent <= 95)
  p2 <- tssa3(s, 4, 6, 6, max_iter = 1000, tol=1e-6)
for (i in 1:4) lines(t.reconstruct(p, list(i))[[1]], col=(rgb(20*i/255, 10*i/255, 1-10*i/255)))
rec <- t.reconstruct(p, list(1, 2:4))
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

tnsr <- rand_tensor(c(6,7,8))
hosvdD <- hosvd(tnsr)
plot(hosvdD$fnorm_resid)
hosvdD2 <- hosvd(tnsr,ranks=c(3,3,4))
plot(hosvdD2$fnorm_resid)
