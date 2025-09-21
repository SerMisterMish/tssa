library(xtable)

rrmse <- function(true, estimate) {
  100 * sqrt(mean(abs(true - estimate)^2)) / abs(true)
}

plot.save <- function(plt, name, format, ..., save = TRUE) {
  print(plt)
  if (save){
    format(name, ...)
    print(plt)
    dev.off()
  }
}

rmse <- function(true, pred) {
  sqrt(mean((pred - true)^2))
}

rmse_mat <- function(true, pred) {
  apply(cbind(true, pred), 1, function(v) rmse(v[1], v[-1]))
}

rmse_ts <- function(true, pred) {
  sqrt(mean(apply(pred, 2, function(p) mean((true - p)^2))))
}

rrmse_freq <- function(true, pred) {
  rmse_mat(true, pred) * 100 / (0.5 * abs(true[1] - true[2]))
}

sd_to_snr <- function(s, sd) {
  10 * log10(sum(abs(s)^2) / sd^2)
}

snr_to_sd <- function(s, snr) {
  sqrt(sum(abs(s)^2) / 10^(snr / 10))
}