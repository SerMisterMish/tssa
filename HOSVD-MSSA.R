### HOSVD-MSSA implementation

library(rTensor)
library(Rssa)
library(purrr)

## Main function, implements the stages of injection and decomposition.
## s must be a matrix with one-dimensional time series as its columns.
## progress = TRUE provides a progress bar for HOSVD step.
## n-ranks of the desired tensor estimate can be provided with ranks.
tmssa3 <- function(s, L = (N + 1) %% 2, progress = FALSE, ranks = hankel.list$hankel.tensor@modes) {
  require("rTensor")
  require("Rssa")
  require("purrr")
  
  hankel.list <- t.hankel(s, L)
  if (progress) {
    hosvd <- rTensor::hosvd(hankel.list$hankel.tensor, ranks = ranks) 
  } else {
    capture.output({ hosvd <- rTensor::hosvd(hankel.list$hankel.tensor, ranks = ranks) }, file = nullfile())
  }
  
  result <- append(hankel.list, list(
                 hosvd = hosvd,
                 modes = hankel.list$hankel.tensor@modes,
                 ranks = ranks))
  result
}

## Injection operator. Returns a list with the trajectory tensor and the
## trajectory matrix in terms of the MSSA method.
t.hankel <- function(s, L) {
  hankel.matrix <- purrr::reduce(apply(s, 2, Rssa::hankel, L = L, simplify = "list"), cbind)
  hankel.tensor <- rTensor::fold(hankel.matrix, 1, 2:3, modes = c(L, nrow(s) + 1 - L, ncol(s)))
  list(hankel.matrix = hankel.matrix, hankel.tensor = hankel.tensor)
}

## Creates a tensor by adding the components of the HOSVD of the trajectory tensor
## with the indices, contained in r1, r2 and r3 parameters, which 
## correspond to the 1-st, 2-nd and 3-rd dimensions respectfully.
tmssa.reconstruct.group <- function(t, r1, r2, r3) {
  apply(
    rTensor::ttl(t$hosvd$Z[r1, r2, r3, drop = FALSE],
        list(t$hosvd$U[[1]][, r1, drop = FALSE],
             t$hosvd$U[[2]][, r2, drop = FALSE],
             t$hosvd$U[[3]][, r3, drop = FALSE]),
        m = 1:3)@data,
    3, Rssa::hankel)
}

## Reconstructs the time series by implementing the steps of grouping and reconstruction
## on the list t, which is obtainable from the tmssa3 function
tmssa.reconstruct <- function(t, groups) {
  lapply(groups, function(r) tmssa.reconstruct.group(t, r[1], r[2], r[3]))
}




