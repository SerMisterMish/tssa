library(rTensor)
library(tensr)
library(purrr)
library(Rssa)
library(Metrics)
library(ggplot2)


Q <- 2
N <- 40
signals <- list(s1 = rep(1, N) %o% c(3, -1.5),
                s2 = cos(2 * pi / 30 * 1:N) %o% c(1, 2))

s <- reduce(signals, `+`)

r <- list(1, 2:3)
r3 <- rep(list(1:2), length(r))

max.r <- 3
Ls <- max.r:(N - max.r + 1)

mult.signal.separability <- function(L, MSSA = TRUE, Tens = TRUE, sigma = 1) {
  if (sigma)
    s.n <- s + rnorm(N * Q, sd = sigma)
  else
    s.n <- s
  
  estimates <- array(dim = c(nrow(s), Q, 2 * length(r)))
  
  if (MSSA) {
    s.mssa <- ssa(s.n, L = L["MSSA"], kind = "mssa")
    for (i in seq_along(r)) {
      estimates[, , i] <- reconstruct(s.mssa, groups = list(r[[i]]))[[1]]
    }
  }
  
  if (Tens) {
    hm <- reduce(apply(s.n, 2, hankel, L = L["Tens"], simplify = "list"), cbind)
    ht <- fold(hm, 1, 2:3, modes = c(L["Tens"], N + 1 - L["Tens"], Q))
    capture.output({ ht.hosvd <- rTensor::hosvd(ht) }, file = nullfile())
    
    for (i in seq_along(r)) {
      estimates[, , i + length(r)] <- apply(
        ttl(ht.hosvd$Z[r[[i]], r[[i]], r3[[i]], drop = FALSE],
            list(ht.hosvd$U[[1]][, r[[i]], drop = FALSE],
                 ht.hosvd$U[[2]][, r[[i]], drop = FALSE],
                 ht.hosvd$U[[3]][, r3[[i]], drop = FALSE]),
            m = 1:3)@data,
        3, hankel)
    }
  }
  estimates
}


## No noise
sigma <- 0
pb <- txtProgressBar(min = max.r, max = N - max.r + 1, initial = 4, style = 3)
res.mssa <- matrix(nrow = Q, ncol = length(Ls))
res.tens <- matrix(nrow = Q, ncol = length(Ls))
for (L.ind in seq_along(Ls)) {
  L <- c(MSSA = Ls[L.ind], Tens = Ls[L.ind])
  mult.sep.res <- mult.signal.separability(L = L, sigma = sigma)
  res.list <- list(MSSA = mult.sep.res[, , seq_along(r)],
                                      Tens = mult.sep.res[, , length(r) + seq_along(r)]) |>
                     lapply(function(x)
                              sapply(seq_along(signals), function(i) mean(sqrt(apply(abs(x[, , i] - signals[[i]])^2, 1:2, mean))))
                     )
  res.mssa[, L.ind] <- res.list$MSSA
  res.tens[, L.ind] <- res.list$Tens
  setTxtProgressBar(pb, L['MSSA'])
}
close(pb)


df.plot <- data.frame(L = Ls, MSSA.const = res.mssa[1,], MSSA.cos = res.mssa[2,],
                      Tens.const = res.tens[1,], Tens.cos = res.tens[2,])
df.plot |> ggplot(aes(x = L, colour = "Method")) + 
  geom_line(aes(y = MSSA.const, colour = "red", linetype = "solid")) + 
  geom_line(aes(y = MSSA.cos, colour = "red", linetype = "dashed")) +
  geom_line(aes(y = Tens.const, colour = "purple", linetype = "solid")) + 
  geom_line(aes(y = Tens.cos, colour = "purple", linetype = "dashed")) +
  xlab("L") + ylab("RMSE") + 
  scale_colour_manual(name = "Method", values = c(`red` = "red", `purple` = "purple"), labels = c("HOSVD-MSSA", "MSSA")) +
  scale_linetype_manual(name = "Signal", values = c(`solid` = "solid", `dashed` = "dashed"), labels = c("Cos", "Const")) 


## Small noise
sigma <- 1
R <- 500
signal.arr <- lapply(signals, function(s) s %o% rep(1, R))
res.mssa.noise1 <- matrix(nrow = Q, ncol = length(Ls))
res.tens.noise1 <- matrix(nrow = Q, ncol = length(Ls))
pb <- txtProgressBar(min = max.r, max = N - max.r + 1, initial = max.r, style = 3)
for (L.ind in seq_along(Ls)) {
  set.seed(1)
  L <- c(MSSA = Ls[L.ind], Tens = Ls[L.ind])
  mult.sep.res <- replicate(R, mult.signal.separability(L = L, sigma = sigma), simplify = "array")
  res.list <- list(MSSA = mult.sep.res[, , seq_along(r),],
                   Tens = mult.sep.res[, , length(r) + seq_along(r),]) |>
  lapply(function(x)
           sapply(seq_along(signals), function(i) mean(sqrt(apply(abs(x[, , i,] - signal.arr[[i]])^2, 1:2, mean))))
  )
  res.mssa.noise1[, L.ind] <- res.list$MSSA
  res.tens.noise1[, L.ind] <- res.list$Tens
  setTxtProgressBar(pb, Ls[L.ind])
}
close(pb)

df.plot2 <- data.frame(L = Ls, MSSA.const = res.mssa.noise1[1,], 
                      MSSA.cos = res.mssa.noise1[2,],
                      Tens.const = res.tens.noise1[1,], 
                      Tens.cos = res.tens.noise1[2,])

df.plot2 |> ggplot(aes(x = L, colour = "Method")) + 
  geom_line(aes(y = MSSA.const, colour = "red", linetype = "solid")) + 
  geom_line(aes(y = MSSA.cos, colour = "red", linetype = "dashed")) +
  geom_line(aes(y = Tens.const, colour = "purple", linetype = "solid")) + 
  geom_line(aes(y = Tens.cos, colour = "purple", linetype = "dashed")) +
  xlab("L") + ylab("RMSE") +
  scale_colour_manual(name = "Method", values = c(`red` = "red", `purple` = "purple"), labels = c("HOSVD-MSSA", "MSSA")) +
  scale_linetype_manual(name = "Signal", values = c(`solid` = "solid", `dashed` = "dashed"), labels = c("Cos", "Const"))


## Large noise
sigma <- 3
R <- 500
signal.arr <- lapply(signals, function(s) s %o% rep(1, R))
res.mssa.noise2 <- matrix(nrow = Q, ncol = length(Ls))
res.tens.noise2 <- matrix(nrow = Q, ncol = length(Ls))
pb <- txtProgressBar(min = max.r, max = N - max.r + 1, initial = max.r, style = 3)
for (L.ind in seq_along(Ls)) {
  set.seed(1)
  L <- c(MSSA = Ls[L.ind], Tens = Ls[L.ind])
  mult.sep.res <- replicate(R, mult.signal.separability(L = L, sigma = sigma), simplify = "array")
  res.list <- list(MSSA = mult.sep.res[, , seq_along(r),],
                   Tens = mult.sep.res[, , length(r) + seq_along(r),]) |>
  lapply(function(x)
      sapply(seq_along(signals), function(i) mean(sqrt(apply(abs(x[, , i,] - signal.arr[[i]])^2, 1:2, mean))))
  )
  res.mssa.noise2[, L.ind] <- res.list$MSSA
  res.tens.noise2[, L.ind] <- res.list$Tens
  setTxtProgressBar(pb, Ls[L.ind])
}
close(pb)

df.plot3 <- data.frame(L = Ls, MSSA.const = res.mssa.noise2[1,], 
                       MSSA.cos = res.mssa.noise2[2,],
                       Tens.const = res.tens.noise2[1,], 
                       Tens.cos = res.tens.noise2[2,])

df.plot3 |> ggplot(aes(x = L, colour = "Method")) + 
  geom_line(aes(y = MSSA.const, colour = "red", linetype = "solid")) + 
  geom_line(aes(y = MSSA.cos, colour = "red", linetype = "dashed")) +
  geom_line(aes(y = Tens.const, colour = "purple", linetype = "solid")) + 
  geom_line(aes(y = Tens.cos, colour = "purple", linetype = "dashed")) +
  xlab("L") + ylab("RMSE") +
  scale_colour_manual(name = "Method", values = c(`red` = "red", `purple` = "purple"), labels = c("HOSVD-MSSA", "MSSA")) +
  scale_linetype_manual(name = "Signal", values = c(`solid` = "solid", `dashed` = "dashed"), labels = c("Cos", "Const")) 


## For component identification
# ht.hosvd$Z@data[1:4, 1:4, 1:2] |> zapsmall()
# par(mfrow = c(3, 2))
# plot(ht.hosvd$U[[1]][,1], type = "l")
# plot(ht.hosvd$U[[1]][,2], type = "l")
# plot(ht.hosvd$U[[1]][,3], type = "l")
# plot(ht.hosvd$U[[2]][,1], type = "l")
# plot(ht.hosvd$U[[2]][,2], type = "l")
# plot(ht.hosvd$U[[2]][,3], type = "l")
# par(mfrow = c(2, 1))
# plot(ht.hosvd$U[[3]][,1], type = "l")
# plot(ht.hosvd$U[[3]][,2], type = "l")
# par(mfrow = c(2, 2))
# plot(ht.hosvd$U[[1]][,1], ht.hosvd$U[[1]][,2], type = "l")
# plot(ht.hosvd$U[[1]][,2], ht.hosvd$U[[1]][,3], type = "l")
# plot(ht.hosvd$U[[2]][,1], ht.hosvd$U[[2]][,2], type = "l")
# plot(ht.hosvd$U[[2]][,2], ht.hosvd$U[[2]][,3], type = "l")
