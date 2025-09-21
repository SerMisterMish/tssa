source("HOSVD-MSSA.R")

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

plot(rec[[1]][,1], type = "l")
lines(signals[[1]][,1], col = "red")

plot(rec[[2]][,1], type = "l")
lines(signals[[2]][,1], col = "red")
