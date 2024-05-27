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
                 s = s,
                 L = L,
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
  lapply(groups, function(r) tmssa.reconstruct.group(t, r[[1]], r[[2]], r[[3]]))
}

## Plots k-mode singular vectors of HOSVD of the trajectory tensor with
## indecies specified in idx[[k]]
plot.tmssa3.vectors <- function(t, idx) {
  
  old.par <- par(mfrow = c(3, max(length(idx[[1]]), 
                                  length(idx[[2]]),
                                  length(idx[[3]]))))
  
  for (k in 1:3) {
    for (i in idx[[k]]) {
      plot(t$hosvd$U[[k]][,i], type = "l", ylab = paste0(k, "-mode singular vector #", i))
    }
  }
  
  par(old.par)
}

## Plots i-th k-mode singular vector of HOSVD of the trajectory tensor
## against i+1-st, where i runs through idx[[k]]
plot.tmssa3.paired <- function(t, idx) {
  old.par <- par(mfrow = c(3, max(length(idx[[1]]), 
                                  length(idx[[2]]),
                                  length(idx[[3]]))))
  
  for (k in 1:3) {
    for (i in idx[[k]]) {
      plot(t$hosvd$U[[k]][,i], t$hosvd$U[[k]][,i + 1], type = "l", 
           xlab = paste0(k, "-mode singular vector #", i),
           ylab = paste0(k, "-mode singular vector #", i+1), 
           main = paste0(i, "-", i+1))
    }
  }
  
  par(old.par)
}

### Examples
library(purrr)
library(rTensor)
library(Rssa)
library(ggplot2)

## Separability of cosine and const
P <- 2
N <- 11
signals <- list(s1 = rep(1, N) %o% c(3, 2), 
                s2 = cos(2 * pi / 3 * 1:N) %o% c(-1, 1.5))
s <- reduce(signals, `+`)

# Parameters
L <- 6
r <- list(1, 2:3)
r3 <- list(1:2, 1:2)
max.r <- 3
max.r3 <- 2

# Creating a list, containing trajectory tensor and its HOSVD
t <- tmssa3(s, L)
# k-mode singular vectors plots
plot.tmssa3.vectors(t, list(1:max.r, 1:max.r, 1:max.r3))
# k-mode singular vectors paired plots
plot.tmssa3.paired(t, list(1:(max.r - 1), 1:(max.r - 1), 1:(max.r3 - 1)))
# Reconstructing const and cosine
rec <- tmssa.reconstruct(t, groups = list(const = list(r[[1]], r[[1]], r3[[1]]), 
                                          cos = list(r[[2]], r[[2]], r3[[2]])))
print(paste0("Const RMSE: ", sqrt(mean((rec$const - signals[[1]])^2))))
print(paste0("Cos RMSE: ", sqrt(mean((rec$cos - signals[[2]])^2))))

# Creating a list, containing trajectory tensor and its HOSVD, using
# a priori information about n-ranks
t <- tmssa3(s, L, ranks = c(max.r, max.r, max.r3))
plot.tmssa3.vectors(t, list(1:max.r, 1:max.r, 1:max.r3))
plot.tmssa3.paired(t, list(1:(max.r - 1), 1:(max.r - 1), 1:(max.r3 - 1)))
rec <- tmssa.reconstruct(t, groups = list(const = list(r[[1]], r[[1]], r3[[1]]), 
                                          cos = list(r[[2]], r[[2]], r3[[2]])))
print(paste0("Const RMSE: ", sqrt(mean((rec$const - signals[[1]])^2))))
print(paste0("Cos RMSE: ", sqrt(mean((rec$cos - signals[[2]])^2))))

## Separability of two cosines with orthogonal 3-rd dimensions
P <- 12
N <- 29
signals <- list(s1 = 2 * cos(2 * pi / 5 * 1:N) %o% cos(2 * pi / 3 * 1:P),
                s2 = cos(2 * pi / 3 * 1:N) %o% (0.5 * cos(2 * pi / 6 * 1:P)))
s <- reduce(signals, `+`)

# Parameters
L <- 15
r <- list(1:2, 3:4)
r3 <- list(1, 2)
max.r <- 4
max.r3 <- 2

t <- tmssa3(s, L)
plot.tmssa3.vectors(t, list(1:max.r, 1:max.r, 1:max.r3))
plot.tmssa3.paired(t, list(1:(max.r - 1), 1:(max.r - 1), 1:(max.r3 - 1)))
rec <- tmssa.reconstruct(t, groups = list(cos5 = list(r[[1]], r[[1]], r3[[1]]), 
                                          cos3 = list(r[[2]], r[[2]], r3[[2]])))
print(paste0("5-period cos RMSE: ", sqrt(mean((rec$cos5 - signals[[1]])^2))))
print(paste0("3-period cos RMSE: ", sqrt(mean((rec$cos3 - signals[[2]])^2))))

## Cosine signal extraction
P <- 12
N <- 44
set.seed(1)
c1 <- rnorm(P)
c2 <- rnorm(P)
signals <- list(s1 = 2 * cos(2 * pi / 5 * 1:N) %o% c1,
                s2 = cos(2 * pi / 3 * 1:N) %o% c2)
s <- reduce(signals, `+`)

# Adding noise to the signal
sd <- 0.5
set.seed(1)
s.noised <- s + rnorm(P * N, sd = sd)

# Parameters
L <- 15
r <- list(1:2, 3:4)
r3 <- list(1:2, 1:2)
max.r <- 4
max.r3 <- 2

# Extracting the signal
t <- tmssa3(s.noised, L)
plot.tmssa3.vectors(t, list(1:max.r, 1:max.r, 1:max.r3))
plot.tmssa3.paired(t, list(1:(max.r - 1), 1:(max.r - 1), 1:(max.r3 - 1)))
rec <- tmssa.reconstruct(t, groups = list(cos5 = list(r[[1]], r[[1]], r3[[1]]), 
                                          cos3 = list(r[[2]], r[[2]], r3[[2]])))
print(paste0("5-period cos RMSE: ", sqrt(mean((rec$cos5 - signals[[1]])^2))))
print(paste0("3-period cos RMSE: ", sqrt(mean((rec$cos3 - signals[[2]])^2))))
