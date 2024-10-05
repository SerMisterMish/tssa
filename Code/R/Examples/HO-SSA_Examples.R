source("HO-SSA.R")

## Separability of sine and constant
N <- 16
x.const <- rep(3, N)
x.sin <- sin(2 * pi * (0:(N-1)) / 3 + pi / 3)
# Creating a list, containing trajectory tensor and its HOSVD
t <- tssa3(x.const + x.sin, I = 6, L = 6)
# 1-mode singular vectors
print(t$hosvd$U[[1]][,1:4] |> zapsmall())
# Core tensor
print(t$hosvd$Z@data[1:4, 1:4, 1:4] |> zapsmall())
# Calculates an estimate of the signal components
t.rec <- t3.reconstruct(t, groups = list(1, 2:3))
# Estimate of the first signal component
plot(t.rec[[1]], type = "l")
# Estimate of the second signal component
plot(t.rec[[2]], type = "l")

## Mixing of cosines
N <- 34
x.cos.3 <- cos(2 * pi * 0:(N - 1) / 3)
x.cos.4 <- cos(2 * pi * 0:(N - 1) / 4)
t <- tssa3(x.cos.3 + x.cos.4, I = 12, L = 12)
# 1-mode singular vectors
print(t$hosvd$U[[1]][,1:5] |> zapsmall())
# Core tensor
print(t$hosvd$Z@data[1:5, 1:5, 1:5] |> zapsmall())

## Sine signal extraction with HOOI
N <- 71
s <- 30 * cos(2 * pi * 1:N / 12)
sigma <- 5
set.seed(1)
x <- s + rnorm(N, sd = sigma)
t <- tssa3(x, I = 7, L = 36, approx.method = "HOOI")
t.rec <- t3.reconstruct(t, groups = 2)

plot(t.rec, type = "l")
lines(s, col = "red")
