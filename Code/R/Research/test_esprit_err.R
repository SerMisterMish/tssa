library(Rssa)
N <- 25
alpha <- -0.02
omega <- 0.2
sig <- exp(alpha*(1:N))*sin(2*pi*(1:N)*omega)

err.alpha <- numeric(0)
err.omega <- numeric(0)
M <- 100
set.seed(42)
for(i in 1:M){
  ts <- sig + rnorm(N, sd = 0.02)
  s <- ssa(ts, L=12)
  par <- parestimate(s, groups = list(1:2))
  err.alpha <- c(err.alpha,(par$rates[1]-alpha)^2)
  err.omega <- c(err.omega,(par$frequencies[2]-omega)^2)
}

#absolute errors are similar
print(sqrt(mean(err.alpha)))
print(sqrt(mean(err.omega)))
print(sqrt(mean(err.omega))*2*pi)

print(100 * sqrt(mean(err.alpha))/alpha)
print(100 * sqrt(mean(err.omega))/omega)
print(sqrt(mean(err.alpha))/alpha / 
      (sqrt(mean(err.omega))/omega))
# print(100 * sqrt(mean(err.omega))/omega*2*pi)
