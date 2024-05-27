### HO-SSA implementation

## Main function, implements the stages of embedding and decomposition
tssa3 <- function(s, I = (length(s) + 2) %/% 3, L = I, approx.method = c("HOSVD", "HOOI")) {
  require("rTensor")
  require("tensr")
  X <- tens3(s, I, L)
  if (length(approx.method) == 2 || approx.method == "HOSVD")
    result <- list(X = X,
                 modes = c(I = I, L = L, J = length(s) - I - L + 2),
                 hosvd = rTensor::hosvd(X),
                 method = "HOSVD")
  else if (approx.method == "HOOI")
    result <- list(X = X,
                   modes = c(I = I, L = L, J = length(s) - I - L + 2),
                   hooi = NA,
                   method = "HOOI")
  result
}

## Embedding operator
tens3 <- function(s, I, L) {
  require("rTensor")
  v <- as.vector(s)
  N <- length(v)
  J <- N - I - L + 2
  X <- array(NA, c(I, L, J))
  for (i in 1:J) {
    X[, , i] <- outer(1:I, 1:L, function(x, y) s[i+ x + y - 2])
  }
  as.tensor(X)
}

## Diagonal averaging of a tensor along i + j + k = const, where const = 3:(I+L+J)
reconstruct.group3 <- function (X.tens) {
  if (length(dim(X.tens)) != 3)
    stop("length(dim(X.tens)) != 3")
  if (class(X.tens) == "Tensor") {
    X <- X.tens@data
  }
  else if (class(X.tens) == "array") {
    X <- X.tens
  }
  
  I <- length(X[,1,1])
  L <- length(X[1,,1])
  J <- length(X[1,1,])
  s <- numeric(I + L + J - 2)
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

## Creates a tensor by adding the components of the HOSVD of the trajectory tensor
## over all dimensions with indices specified in group
make.group <- function (hosvd, group) {
  ttl(hosvd$Z[group,group,group,drop=FALSE], list(
    as.matrix(hosvd$U[[1]][,group]),
    as.matrix(hosvd$U[[2]][,group]),
    as.matrix(hosvd$U[[3]][,group])),
    1:3)
}

## Creates a tensor by adding the components of the HOSVD of the trajectory tensor
## with indices corresponding to certain dimension specified in corresponding
## element of a group list
make.group.arb <- function (hosvd, group) {
  ttl(hosvd$Z[group[[1]], group[[2]], group[[3]], drop=FALSE], 
      list(
        hosvd$U[[1]][, group[[1]], drop = FALSE], 
        hosvd$U[[2]][, group[[2]], drop = FALSE], 
        hosvd$U[[3]][, group[[3]], drop = FALSE]), 1:3)
}

## Calculates the (r, r, r)-rank approximation of the trajectory tensor 
## using HOOI algorithm
make.group.HOOI <- function(t, r) {
  if (t$method != "HOOI")
    warning("Used HOOI related function for a list without a hooi element")
  t$hooi <- hooi(t$X@data, r = c(r, r, r))
  t$approx <- atrans(t$hooi$G, t$hooi$U)
  t$approx
}

## Reconstructs the time series by implementing the steps of grouping and reconstruction
## on the list t, which is obtainable from the tssa3 function
t3.reconstruct <- function(t, groups) {
  if (t$method == "HOOI")
  {
    stopifnot(is.numeric(groups))
    reconstruct.group3(make.group.HOOI(t, groups))
  }
  else
  {
    stopifnot(is.list(groups))
    lapply(lapply(groups, make.group, hosvd = t$hosvd), reconstruct.group3)
  }
}

