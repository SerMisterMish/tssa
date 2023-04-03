tssa3 <- function(s, I = (length(s) + 2) %/% 3, L = I) {
  X <- tens3(s, I, L)
  result <- list()
  result$X <- X
  result$modes <- list(I = I, L = L, J = length(s)-I-L+2)
  result$hosvd <- hosvd(X)
  return(result)
}

tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N-I-L + 2
  X <- array(NA, c(I, L, J))
  print(c(I, L, J))
  for (i in 1:J) {
    X[, , i] <- outer(1:I, 1:L, function(x, y) s[i+ x + y - 2])
  }
  return(as.tensor(X))
}

# diagonal averaging i + j + k = const, const = 3:(I+L+J)
reconstruct.group3 <- function (X.tens) {
  X <- X.tens@data
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

make.group <- function (hosvd, group) {
  ttl(hosvd$Z[group,group,group,drop=FALSE], list(
    as.matrix(hosvd$U[[1]][,group]),
    as.matrix(hosvd$U[[2]][,group]),
    as.matrix(hosvd$U[[3]][,group])),
    1:3)
}

t3.reconstruct <- function(p, groups) {
  stopifnot(is.list(groups))
  lapply(lapply(groups, make.group, hosvd = p$hosvd), reconstruct.group3)
}


# Test 1: const
s <- rep(3, 200)
p <- tssa3(s)
rec <- t3.reconstruct(p, list(1))
mse(rec[[1]], s)

# Test 2: sin
s <- sin(2 * pi * 0:120 / 3 + pi / 3)
p <- tssa3(s, 39, 39)
rec <- t3.reconstruct(p, list(1:2))
mse(rec[[1]], s)

# Test 3: const + sin
s.const <- 3
s.sin <- sin(2 * pi * 0:120 / 3)
s <- s.const + s.sin
p <- tssa3(s, 39, 39)
rec <- t3.reconstruct(p, list(1, 2:3))
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
  capture.output({ s.tssa <- tssa3(s, 4, 2) })
  s.trec <- t3.reconstruct(s.tssa, groups = list(1:r))
  s.rec <- reconstruct(s.ssa, groups = list(1:r))
  tssa.deviation <- c(tssa.deviation, mse(s.sin, s.trec[[1]]))
  ssa.deviation <- c(ssa.deviation, mse(s.sin, s.rec$F1))
  setTxtProgressBar(pb,i)
}
close(pb)

print(c(" ssa  ", round(sqrt(mean(ssa.deviation)), digits = 4)), quote = FALSE)
print(c(" tssa ", round(sqrt(mean(tssa.deviation)), digits = 4)), quote = FALSE)

