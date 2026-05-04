if (!require("Rssa"))
  install.packages("Rssa")
if (!require("rTensor"))
  install.packages("rTensor")
if (!require("Rcpp"))
  install.packages("Rcpp")
if (!require("RcppArmadillo"))
  install.packages("RcppArmadillo")

library(Rssa)
library(rTensor)
library(Rcpp)
library(RcppArmadillo)

if (!(exists("reconstruct_group_t3", mode = "function") &&
      exists("reconstruct_group_tn", mode = "function"))) {
  if (!require("TssaCppDev")) {
    if (!require("devtools")) {
      install.packages("devtools")
      library(devtools)
    }
    oldwd <- getwd()
    setwd("TssaCppDev")
    Rcpp::compileAttributes(".")
    devtools::install()
    setwd(oldwd)
    library(TssaCppDev)
  }
}

get_unfolds_t3 <- function(x, L, rm.repeated = FALSE) {
  N <- length(x)
  I <- L[1]
  L <- L[2]
  K <- N - I - L + 2
  unfolds_r <- list()
  if (!rm.repeated) {
    h12 <- hankel(x, I + L - 1)
    h3 <- hankel(x, I + K - 1)
    sR12 <- Re(h12)
    sR3 <- Re(h3)
    
    unfolds_r[[1]] <- Rssa::new.hbhmat(sR12, c(I, 1))
    unfolds_r[[2]] <- Rssa::new.hbhmat(sR12, c(L, 1))
    unfolds_r[[3]] <- Rssa::new.hbhmat(sR3, c(K, 1))
  } else {
    unfolds_r[[1]] <- Rssa::new.hmat(Re(x), I)
    unfolds_r[[2]] <- Rssa::new.hmat(Re(x), L)
    unfolds_r[[3]] <- Rssa::new.hmat(Re(x), K)
  }
  
  if (is.complex(x)) {
    unfolds_i <- list()
    
    if (!rm.repeated) {
      sI12 <- Im(h12)
      sI3 <- Im(h3)
      
      unfolds_i[[1]] <- Rssa::new.hbhmat(sI12, c(I, 1))
      unfolds_i[[2]] <- Rssa::new.hbhmat(sI12, c(L, 1))
      unfolds_i[[3]] <- Rssa::new.hbhmat(sI3, c(K, 1))
    } else {
      unfolds_i[[1]] <- Rssa::new.hmat(Im(x), I)
      unfolds_i[[2]] <- Rssa::new.hmat(Im(x), L)
      unfolds_i[[3]] <- Rssa::new.hmat(Im(x), K)
    }
  } else unfolds_i <- NULL
  invisible(list(Re = unfolds_r, Im = unfolds_i))
}

setMethod("get_unfolds_t3", "Tensor", function(x, L = NULL, rm.repeated = NULL) {
  unfolds <- lapply(1:3, rs_unfold, tnsr = x)
  unfolds_r <- lapply(unfolds, \(u) Re(u@data))
  if (is.complex(x@data)) unfolds_i <- lapply(unfolds, \(u) Im(u@data)) else unfolds_i <- NULL
  invisible(list(Re = unfolds_r, Im = unfolds_i))
})

get_unfolds <- function(x, L, rm.repeated = TRUE) {
  N <- length(x)
  n <- length(L)
  stopifnot(is.numeric(L) && n > 1)
  Lsum <- sum(L)
  L <- c(L, N - Lsum + n)
  n <- n + 1
  
  unfolds_r <- list()
  if (!rm.repeated) {
    hn <- hankel(x, Lsum - 1)
    hlast <- hankel(x, sum(L[-(n - 1)]) - 1)

    unfolds_r <- lapply(L[-n], \(l) Rssa::new.hbhmat(Re(hn), c(l, 1)))
    unfolds_r[[n]] <- Rssa::new.hbhmat(Re(hlast), c(L[n], 1))
  } else {
    unfolds_r <- lapply(L, FUN = Rssa::new.hmat, F = Re(x))
  }
  
  if (is.complex(x)) {
    unfolds_i <- list()
    
    if (!rm.repeated) {
      unfolds_i <- lapply(L[-n], \(l) Rssa::new.hbhmat(Im(hn), c(l, 1)))
      unfolds_i[[n]] <- Rssa::new.hbhmat(Im(hlast), c(L[n], 1))
    } else {
      unfolds_i <- lapply(L, FUN = Rssa::new.hmat, F = Im(x))
    }
  } else unfolds_i <- NULL
  invisible(list(Re = unfolds_r, Im = unfolds_i))
}

# Embedding series into tensor
tens3 <- function(x, L, kind = c("SSA", "MSSA", "CP"), rm.repeated = FALSE) {
  kind <- toupper(kind)
  kind <- match.arg(kind)
  
  if (kind == "SSA") {
    stopifnot(length(L) == 2)
    N <- length(x)
    I <- L[1]
    L <- L[2]
    K <- N - I - L + 2
    X <- outer(1:I, 1:L, `+`) |>
      outer(1:K, function(il, j)
        x[il + j - 2]) |> as.tensor()
    
    unf <- get_unfolds_t3(x, c(I, L), rm.repeated)
    attr(X, "unfolds_r") <- unf$Re
    
    if (is.complex(x)) {
      attr(X, "unfolds_i") <- unf$Im
    }
  }
  else if (kind == "CP") {
    stopifnot(length(L) == 2)
    N <- length(x)
    l <- L[1]
    J <- L[2]
    
    I <- N %/% l
    K <- l - J + 1
    X_mat <- matrix(x[1:(I * l)], nrow = I, byrow = TRUE) |>
      apply(1, function(v)
        outer(1:J, 1:K, function(k, j)
          v[k + j - 1])) |>
      as.matrix()
    dim(X_mat) <- c(J * K, I)
    X <- rTensor::fold(X_mat, 2:3, 1, modes = c(I, J, K))
  }
  else if (kind == "MSSA") {
    N <- nrow(x)
    K <- N - L + 1
    Q <- ncol(x)
    X <- apply(x, 2, Rssa::hankel, L = L, simplify = FALSE) |>
      Reduce(cbind, x = _) |>
      rTensor::fold(1, 2:3, modes = c(L, K, Q))
    
    sR <- Re(x)
    unfolds_r <- list()
    
    unfolds_r[[1]] <- Rssa::new.hbhmat(sR, c(L, 1))
    unfolds_r[[2]] <- Rssa::new.hbhmat(sR, c(K, 1))
    
    mulR <- function(v)
      as.numeric(t(sR) %*% v)
    tmulR <- function(v)
      as.numeric(sR %*% v)
    unfolds_r[[3]] <- extmat(mulR, tmulR, nrow = Q, ncol = N)
    
    attr(X, "unfolds_r") <- unfolds_r
    
    if (is.complex(x)) {
      sI <- Im(x)
      unfolds_i <- list()
      
      unfolds_i[[1]] <- Rssa::new.hbhmat(sI, c(L, 1))
      unfolds_i[[2]] <- Rssa::new.hbhmat(sI, c(K, 1))
      
      mulI <- function(v)
        as.numeric(t(sI) %*% v)
      tmulI <- function(v)
        as.numeric(sI %*% v)
      unfolds_i[[3]] <- extmat(mulI, tmulI, nrow = Q, ncol = N)
      
      attr(X, "unfolds_i") <- unfolds_i
    }
  }
  return(X)
}

tens <- function(x, L, kind = c("SSA", "MSSA", "CP"), rm.repeated = FALSE) {
  kind <- toupper(kind)
  kind <- match.arg(kind)
  
  if (kind == "SSA") {
    n <- length(L)
    stopifnot(is.numeric(L) && n > 1)
    if (n == 2) return(tens3(x, L, kind, rm.repeated))
    N <- length(x)
    Lsum <- sum(L)
    L <- c(L, N - Lsum + n)
    n <- n + 1
    
    idx <- Reduce(\(l, k) outer(l, 1:L[k], `+`), 3:n, outer(1:L[1], 1:L[2], `+`)) - n + 1
    X <- x[idx]
    dim(X) <- dim(idx)
    X <- as.tensor(X)
    
    unf <- get_unfolds(x, L[-n], rm.repeated)
    attr(X, "unfolds_r") <- unf$Re
    
    if (is.complex(x)) {
      attr(X, "unfolds_i") <- unf$Im
    }
  }
  else if (kind == "CP") {
    stopifnot(length(L) == 2)
    N <- length(x)
    l <- L[1]
    J <- L[2]
    
    I <- N %/% l
    K <- l - J + 1
    X_mat <- matrix(x[1:(I * l)], nrow = I, byrow = TRUE) |>
      apply(1, function(v)
        outer(1:J, 1:K, function(k, j)
          v[k + j - 1])) |>
      as.matrix()
    dim(X_mat) <- c(J * K, I)
    X <- rTensor::fold(X_mat, 2:3, 1, modes = c(I, J, K))
  }
  else if (kind == "MSSA") {
    N <- nrow(x)
    K <- N - L + 1
    Q <- ncol(x)
    X <- apply(x, 2, Rssa::hankel, L = L, simplify = FALSE) |>
      Reduce(cbind, x = _) |>
      rTensor::fold(1, 2:3, modes = c(L, K, Q))
    
    sR <- Re(x)
    unfolds_r <- list()
    
    unfolds_r[[1]] <- Rssa::new.hbhmat(sR, c(L, 1))
    unfolds_r[[2]] <- Rssa::new.hbhmat(sR, c(K, 1))
    
    mulR <- function(v)
      as.numeric(t(sR) %*% v)
    tmulR <- function(v)
      as.numeric(sR %*% v)
    unfolds_r[[3]] <- extmat(mulR, tmulR, nrow = Q, ncol = N)
    
    attr(X, "unfolds_r") <- unfolds_r
    
    if (is.complex(x)) {
      sI <- Im(x)
      unfolds_i <- list()
      
      unfolds_i[[1]] <- Rssa::new.hbhmat(sI, c(L, 1))
      unfolds_i[[2]] <- Rssa::new.hbhmat(sI, c(K, 1))
      
      mulI <- function(v)
        as.numeric(t(sI) %*% v)
      tmulI <- function(v)
        as.numeric(sI %*% v)
      unfolds_i[[3]] <- extmat(mulI, tmulI, nrow = Q, ncol = N)
      
      attr(X, "unfolds_i") <- unfolds_i
    }
  }
  return(X)
}

reconstruct.group3 <- function(X.tens, kind = c("SSA", "MSSA", "CP")) {
  stopifnot(is(X.tens, "Tensor"))
  X <- X.tens@data
  kind <- toupper(kind)
  kind <- match.arg(kind)
  if (kind == "SSA") {
    s <- reconstruct_group_t3(X)
  } else if (kind == "CP") {
    s <- Reduce(c, lapply(seq(dim(X)[1]), function(i)
      hankel(as.matrix(X[i, , ]))))
  } else if (kind == "MSSA") {
    s <- Reduce(cbind, apply(X, 3, Rssa::hankel, simplify = FALSE))
  }
  if (is.complex(X)) return(s) else return(Re(s))
}

reconstruct.group <- function(X.tens, kind = c("SSA", "MSSA", "CP")) {
  stopifnot(is(X.tens, "Tensor"))
  X <- X.tens@data
  modes <- X.tens@modes
  kind <- toupper(kind)
  kind <- match.arg(kind)
  if (kind == "SSA") {
    s <- reconstruct_group_tn(X)
  } else if (kind == "CP") {
    s <- Reduce(c, lapply(seq(dim(X)[1]), function(i)
      hankel(as.matrix(X[i, , ]))))
  } else if (kind == "MSSA") {
    s <- Reduce(cbind, apply(X, 3, Rssa::hankel, simplify = FALSE))
  }
  if (is.complex(X)) return(s) else return(Re(s))
}

.pad0 <- function(x, N) c(x, rep(0, N - length(x)))
.mvpad0 <- function(x, N) rbind(x, matrix(0, nrow = N - nrow(x), ncol = ncol(x)))

reconstruct.rank1.from.decomp <- function(dec, idx) {
  modes <- attr(dec, "modes")
  if (is.null(modes)) modes <- sapply(dec$U, dim)[1,]
  nm <- length(modes)
  stopifnot(nm == length(idx) && all(modes >= idx))
  
  N <- attr(dec, "N")
  if (is.null(N)) N <- sum(modes) - nm + 1
  
  us_fft <- attr(dec, "us_fft")
  if (is.null(us_fft)) {
    u_fft <- lapply(seq(nm), \(j) fft(.pad0(dec$U[[j]][, idx[j]], N)))
  } else {
    u_fft <- lapply(seq(nm), \(j) us_fft[[j]][, idx[j]])
  }
  
  gp <- fft(Reduce(`*`, u_fft), inverse = TRUE)
  
  # TODO: find formula for weights
  w <- attr(dec, "w")
  if (is.null(w)) w <- 1 / as.vector(table(tens(seq(N), modes[-nm])@data)) / N
  
  sig <- do.call(`[`, c(list(dec$Z@data), as.list(idx)))
  
  if (is.complex(dec$Z@data)) return(sig * gp * w) else return(Re(sig * gp * w))
}

reconstruct.group.from.decomp <- function(dec, group) {
  w <- attr(dec, "w")
  if (is.null(w)) {
    modes <- sapply(dec$U, dim)[1,]
    nm <- length(dec$U)
    N <- sum(modes) - nm + 1
    w <- 1 / as.vector(table(tens(seq(N), modes[-nm])@data)) / N
  }
  
  group_grid <- expand.grid(group)
  res <- 0
  for (i in seq(nrow(group_grid))) {
    res <- res + reconstruct.rank1.from.decomp(dec, unlist(group_grid[i,]))
  }
  res
}

# Complex norm
fnorm_complex <- function(x) {
  sqrt(sum(Re(x * Conj(x))))
}

setMethod("fnorm_complex", "Tensor", function(x) {
  fnorm_complex(x@data)
})

# Calculate tensor SINGULAR values (sorry for naming)
get_tensor_eigenvalues <- function(X, dims = seq(num_dims)) {
  num_dims <- length(dim(X))
  
  lapply(dims, function(d) {
    mask <- as.list(rep(TRUE, num_dims))
    sapply(seq(dim(X)[d]), function(i) {
      mask[[d]] <- i
      fnorm_complex(do.call(`[`, c(list(X), mask, list(drop = FALSE))))
    })
  })
}

setMethod("get_tensor_eigenvalues", "Tensor", function(X, dims = seq(X@num_modes)) {
  get_tensor_eigenvalues(X@data, dims)
})

setMethod("get_tensor_eigenvalues", "list", function(X, dims = seq(X$Z@num_modes)) {
  get_tensor_eigenvalues(X$Z, dims)
})

# HOSVD and HOOI modifications for complex cases
hosvd_mod <- function(tnsr,
                      ranks = NULL,
                      svd.method = c("svd", "primme"),
                      status = FALSE,
                      calc.est = FALSE)
{
  stopifnot(is(tnsr, "Tensor"))
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  if (all(tnsr@data == 0))
    stop("Zero tensor detected")
  num_modes <- tnsr@num_modes
  if (is.null(ranks)) {
    ranks <- tnsr@modes
  }
  else {
    if (sum(ranks > tnsr@modes) != 0)
      stop("ranks must be smaller than the corresponding mode")
  }
  svd.method <- match.arg(svd.method)
  if (identical(svd.method, "primme") &&
      !requireNamespace("PRIMME", quietly = TRUE))
    stop("PRIMME package is required for SVD method `primme'")
  
  if (status)
    pb <- txtProgressBar(min = 0, max = num_modes, style = 3)
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    R <- attr(tnsr, "unfolds_r")[[m]]
    Ilist <- attr(tnsr, "unfolds_i")
    if (!is.null(Ilist))
      I <- Ilist[[m]]
    else
      I <- NULL
    
    if (identical(svd.method, "svd")) {
      if (is.null(I)) {
        unfld <- as.matrix(R)
      } else {
        unfld <- as.matrix(R) + 1i * as.matrix(I)
      }
      
      U_list[[m]] <- svd(unfld, nu = ranks[m])$u
    } else if (identical(svd.method, "primme")) {
      matmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = FALSE)
        else
          ematmul(x, y, transposed = FALSE)
      }
      
      tmatmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = TRUE)
        else
          ematmul(x, y, transposed = TRUE)
      }
      
      if (is.null(I)) {
        A <- function(x, trans) {
          if (identical(trans, "c"))
            tmatmul(R, x)
          else
            matmul(R, x)
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = min(ranks[m], ncol(R)),
          m = nrow(R),
          n = ncol(R),
          isreal = TRUE
        )$u
      } else {
        A <- function(x, trans) {
          rX <- Re(x)
          iX <- Im(x)
          if (identical(trans, "c")) {
            (tmatmul(R, rX) + tmatmul(I, iX)) +
              1i * (tmatmul(R, iX) - tmatmul(I, rX))
          } else {
            (matmul(R, rX) - matmul(I, iX)) +
              1i * (matmul(R, iX) + matmul(I, rX))
          }
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = min(ranks[m], ncol(R)),
          m = nrow(R),
          n = ncol(R),
          isreal = FALSE
        )$u
      }
    } else
      stop("Unsupported svd method")
    
    if (status)
      setTxtProgressBar(pb, m)
  }
  if (status)
    close(pb)
  Z <- ttl(tnsr, lapply(U_list, (\(.) Conj(t(
    .
  )))), ms = 1:num_modes)
  if (calc.est) {
    est <- ttl(Z, U_list, ms = 1:num_modes)
    resid <- fnorm_complex(est - tnsr)
    ret <- list(
      Z = Z,
      U = U_list,
      est = est,
      fnorm_resid = resid
    )
  } else {
    ret <- list(
      Z = Z,
      U = U_list
    ) 
  }
  ret
}

tucker_mod <- function(tnsr,
                       ranks = NULL,
                       svd.method = c("svd", "primme"),
                       max_iter = 25,
                       tol = 1e-05,
                       status = FALSE)
{
  stopifnot(is(tnsr, "Tensor"))
  if (is.null(ranks))
    stop("ranks must be specified")
  if (sum(ranks > tnsr@modes) != 0)
    stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  if (all(tnsr@data == 0))
    stop("Zero tensor detected")
  svd.method <- match.arg(svd.method)
  if (identical(svd.method, "primme") &&
      !requireNamespace("PRIMME", quietly = TRUE))
    stop("PRIMME package is required for SVD method `primme'")
  
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    R <- attr(tnsr, "unfolds_r")[[m]]
    Ilist <- attr(tnsr, "unfolds_i")
    if (!is.null(Ilist))
      I <- Ilist[[m]]
    else
      I <- NULL
    
    if (identical(svd.method, "svd")) {
      if (is.null(I)) {
        unfld <- as.matrix(R)
      } else {
        unfld <- as.matrix(R) + 1i * as.matrix(I)
      }
      
      U_list[[m]] <- svd(unfld, nu = ranks[m])$u
    } else if (identical(svd.method, "primme")) {
      
      matmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = FALSE)
        else
          ematmul(x, y, transposed = FALSE)
      }
      
      tmatmul <- function(x, y) {
        if (is.matrix(y))
          apply(y, 2, ematmul, emat = x, transposed = TRUE)
        else
          ematmul(x, y, transposed = TRUE)
      }
      
      if (is.null(I)) {
        A <- function(x, trans) {
          if (identical(trans, "c"))
            tmatmul(R, x)
          else
            matmul(R, x)
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = ranks[m],
          m = nrow(R),
          n = ncol(R),
          isreal = TRUE
        )$u
      } else {
        A <- function(x, trans) {
          rX <- Re(x)
          iX <- Im(x)
          if (identical(trans, "c")) {
            (tmatmul(R, rX) + tmatmul(I, iX)) +
              1i * (tmatmul(R, iX) - tmatmul(I, rX))
          } else {
            (matmul(R, rX) - matmul(I, iX)) +
              1i * (matmul(R, iX) + matmul(I, rX))
          }
        }
        U_list[[m]] <- PRIMME::svds(
          A,
          NSvals = ranks[m],
          m = nrow(R),
          n = ncol(R),
          isreal = FALSE
        )$u
      }
    } else
      stop("Unsupported svd method")
  }
  tnsr_norm <- fnorm_complex(tnsr)
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(Z, U_list) {
    est <- ttl(Z, U_list, ms = 1:num_modes)
    curr_resid <- fnorm_complex(tnsr - est)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter == 1)
      return(FALSE)
    if (abs(curr_resid - fnorm_resid[curr_iter - 1]) / tnsr_norm <
        tol)
      return(TRUE)
    else {
      return(FALSE)
    }
  }
  
  if (status)
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  while ((curr_iter < max_iter) && (!converged)) {
    if (status)
      setTxtProgressBar(pb, curr_iter)
    modes <- tnsr@modes
    modes_seq <- 1:num_modes
    for (m in modes_seq) {
      X <- ttl(tnsr, lapply(U_list[-m], (\(.) Conj(t(.)))), ms = modes_seq[-m])
      
      # TODO: Change this part so it supports primme and pseudo-hosvd
      U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
    }
    Z <- ttm(X, mat = Conj(t(U_list[[num_modes]])), m = num_modes)
    if (CHECK_CONV(Z, U_list)) {
      converged <- TRUE
      if (status)
        setTxtProgressBar(pb, max_iter)
    }
    else {
      curr_iter <- curr_iter + 1
    }
  }
  if (status)
    close(pb)
  
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) *
    100
  est <- ttl(Z, U_list, ms = 1:num_modes)
  invisible(
    list(
      Z = Z,
      U = U_list,
      conv = converged,
      iter = curr_iter,
      est = est,
      norm_percent = norm_percent,
      fnorm_resid = tail(fnorm_resid, 1),
      all_resids = fnorm_resid
    )
  )
}

# CPD modification for complex cases

cp_mod <- function(tnsr,
                   num_components = NULL,
                   max_iter = 25,
                   tol = 1e-05,
                   start = c("rand", "svd"),
                   status = FALSE)
{
  if (is.null(num_components))
    stop("num_components must be specified")
  stopifnot(is(tnsr, "Tensor"))
  if (all(tnsr@data == 0))
    stop("Zero tensor detected")
  start <- match.arg(start)
  num_modes <- tnsr@num_modes
  modes <- tnsr@modes
  U_list <- vector("list", num_modes)
  unfolded_mat <- vector("list", num_modes)
  tnsr_norm <- fnorm_complex(tnsr)
  for (m in 1:num_modes) {
    unfolded_mat[[m]] <- rs_unfold(tnsr, m = m)@data
    if (identical(start, "rand")) {
      U_list[[m]] <- matrix(rnorm(modes[m] * num_components),
                            nrow = modes[m],
                            ncol = num_components)
    } else if (identical(start, "svd")) {
      U_list[[m]] <- svd(unfolded_mat[[m]], nu = num_components)$u
    } else
      stop("Unsupported `start` parameter value")
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est) {
    curr_resid <- fnorm_complex(est - tnsr)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter == 1)
      return(FALSE)
    if (abs(curr_resid - fnorm_resid[curr_iter - 1]) / tnsr_norm < tol)
      return(TRUE)
    else {
      return(FALSE)
    }
  }
  
  if (status)
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  
  norm_vec <- function(vec) {
    norm(as.matrix(vec))
  }
  while ((curr_iter < max_iter) && (!converged)) {
    if (status)
      setTxtProgressBar(pb, curr_iter)
    
    for (m in 1:num_modes) {
      V <- hadamard_list(lapply(U_list[-m], function(x) {
        t(x) %*% x
      }))
      V_inv <- solve(V)
      tmp <- unfolded_mat[[m]] %*% khatri_rao_list(U_list[-m], reverse = TRUE) %*% V_inv
      lambdas <- apply(tmp, 2, norm_vec)
      U_list[[m]] <- sweep(tmp, 2, lambdas, "/")
      Z <- rTensor:::.superdiagonal_tensor(num_modes = num_modes,
                                           len = num_components,
                                           elements = lambdas)
      est <- ttl(Z, U_list, ms = 1:num_modes)
    }
    if (CHECK_CONV(est)) {
      converged <- TRUE
      if (status)
        setTxtProgressBar(pb, max_iter)
    }
    else {
      curr_iter <- curr_iter + 1
    }
  }
  if (status && !converged) {
    setTxtProgressBar(pb, max_iter)
  }
  if (status)
    close(pb)
  
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) *
    100
  invisible(
    list(
      lambdas = lambdas,
      U = U_list,
      conv = converged,
      est = est,
      norm_percent = norm_percent,
      fnorm_resid = tail(fnorm_resid, 1),
      all_resids = fnorm_resid
    )
  )
}

# Partially reconstruct CPD
cp_reconstruct_part <- function(cp, ind) {
  num_modes <- length(cp$U)
  ttl(
    rTensor:::.superdiagonal_tensor(
      num_modes = num_modes,
      len = length(ind),
      elements = cp$lambdas[ind]
    ),
    list(cp$U[[1]][, ind, drop = FALSE], cp$U[[2]][, ind, drop = FALSE], cp$U[[3]][, ind, drop = FALSE]),
    seq(num_modes)
  )
}

# TSVD modification
tsvd_mod <- function(tnsr, status = FALSE)
{
  if (tnsr@num_modes != 3)
    stop("T-SVD only implemented for 3d so far")
  if (all(tnsr@data == 0))
    stop("Zero tensor detected")
  modes <- tnsr@modes
  n1 <- modes[1]
  n2 <- modes[2]
  n3 <- modes[3]
  if (status)
    pb <- txtProgressBar(min = 0, max = n3, style = 3)
  fftz <- aperm(apply(tnsr@data, MARGIN = 1:2, fft), c(2, 3, 1))
  U_arr <- array(0, dim = c(n1, n1, n3))
  V_arr <- array(0, dim = c(n2, n2, n3))
  m <- min(n1, n2)
  S_arr <- array(0, dim = c(n1, n2, n3))
  for (j in 1:n3) {
    if (status)
      setTxtProgressBar(pb, j)
    decomp <- svd(fftz[, , j], nu = n1, nv = n2)
    U_arr[, , j] <- decomp$u
    V_arr[, , j] <- decomp$v
    S_arr[, , j] <- diag(decomp$d, nrow = n1, ncol = n2)
  }
  if (status)
    close(pb)
  
  IFFT <- function(x) {
    as.numeric(fft(x, inverse = TRUE)) / length(x)
  }
  U <- as.tensor(aperm(apply(U_arr, MARGIN = 1:2, rTensor:::.ifft), c(2, 3, 1)))
  V <- as.tensor(aperm(apply(V_arr, MARGIN = 1:2, rTensor:::.ifft), c(2, 3, 1)))
  S <- as.tensor(aperm(apply(S_arr, MARGIN = 1:2, rTensor:::.ifft), c(2, 3, 1)))
  
  invisible(list(U = U, V = V, S = S))
}

# Partially reconstruct TSVD
tsvd_reconstruct_part <- function(tsvd, group) {
  t_mult(
    t_mult(
      tsvd$U[, group, , drop = FALSE], 
      tsvd$S[group, group, , drop = FALSE]),
    t(tsvd$V[, group, , drop = FALSE]))
}

# HO-ESPRIT
tens_esprit <- function(s,
                        L,
                        groups,
                        kind = c("TSSA", "TMSSA"),
                        decomp = c("HOOI", "HOSVD"),
                        svd.method = c("svd", "primme"),
                        rm.repeated = FALSE,
                        est_dim,
                        r3 = NULL,
                        status = FALSE,
                        qrtol = 1e-07)
{
  max_rank <- max(sapply(groups, max))
  
  kind <- toupper(kind)
  kind <- match.arg(kind)
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  is.decomposition <- attr(s, "is.decomposition")
  
  if (is.null(is.decomposition) || !is.decomposition) {
    if (identical(kind[1], "TSSA"))
      H <- tens(s, L, kind = "SSA", rm.repeated = rm.repeated)
    else
    {
      if (is.null(r3)) {
        simpleWarning("r3 argument was not provided, setting
                      r3 as maximum across groups")
        r3 <- max_rank
      }
      H <- tens(s, L, kind = "MSSA")
    }
    
    max_ranks <- pmin(max_rank, H@modes)
    
    if (identical(kind[1], "TMSSA")) {
      max_ranks[3] <- min(max_ranks[3], r3)
      if (identical(decomp[1], "HOOI")) {
        H.decomp <- tucker_mod(H, max_ranks, status = status)
      } else {
        H.decomp <- hosvd_mod(
          H,
          max_ranks,
          svd.method = svd.method,
          status = status
        )
      }
    } else if (identical(decomp[1], "HOOI")) {
      H.decomp <- tucker_mod(H, max_ranks, status = status)
    } else {
      H.decomp <- hosvd_mod(H,
                            max_ranks,
                            svd.method = svd.method,
                            status = status)
    }
  } else {
    H.decomp <- s
  }
  
  estimates <- list()
  if (is.list(est_dim)) est_dim <- unlist(est_dim)

  for (i in seq(groups)) {
    poles <- sapply(est_dim, \(ed) {
      U <- H.decomp$U[[ed]][, groups[[i]], drop = FALSE]
      Z <- qr.solve(U[-nrow(U), ], U[-1, ], qrtol)
      as.complex(eigen(Z, only.values = TRUE)$values)
    })
    if (is.matrix(poles)) {
      poles <- rowMeans(poles)
    }
    estimates[[i]] <- list(
      poles = poles,
      rates = Re(log(poles)),
      frequencies = Im(log(poles)) / 2 / pi
    )
  }
  estimates
}

# TSSA
tens_ssa_decompose <- function(s,
                               L,
                               neig = NULL,
                               decomp = c("HOSVD", "HOOI", "CP"),
                               svd.method = c("svd", "primme"),
                               rm.repeated = FALSE,
                               status = FALSE,
                               cp_span = c("mean", "last"),
                               ...) {
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  if (decomp == "CP") {
    if (is.null(neig))
      stop("For CP decomposition `neig` argument must be provided")
    H <- tens(s, L, kind = "CP")
  }
  else
    H <- tens(s, L, kind = "SSA", rm.repeated = rm.repeated)
  
  
  if (is.null(neig))
    neig <- 50
  
  if (decomp != "CP" && length(neig) == 1)
    neig <- pmin(neig, c(L, length(s) - sum(L) + length(L)))
  
  H.dec <- switch(
    decomp,
    HOSVD = hosvd_mod(
      H,
      ranks = neig,
      svd.method = svd.method,
      status = status,
      ...
    ),
    HOOI = tucker_mod(H, ranks = neig, status = status, ...),
    CP = cp_mod(H, num_components = neig, status = status, ...)
  )
  attr(H.dec, "is.decomposition") <- TRUE
  H.dec
}

tens_ssa_reconstruct <- function(s,
                                 L,
                                 groups,
                                 decomp = c("HOSVD", "HOOI", "CP"),
                                 svd.method = c("svd", "primme"),
                                 rm.repeated = FALSE,
                                 rec.from.decomp = TRUE,
                                 trunc_dims = NULL,
                                 status = FALSE,
                                 cp_span = c("mean", "last"),
                                 ...) {
  is.decomposition <- attr(s, "is.decomposition")
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  if (is.null(is.decomposition) || !is.decomposition) {
    if (!is(s, "Tensor")) {
      if (decomp == "CP")
        H <- tens(s, L, kind = "CP")
      else
        H <- tens(s, L, kind = "SSA", rm.repeated = rm.repeated)
    } else {
      H <- s
      L <- H@modes[-(H@num_modes)]
      if (is.null(attr(H, "unfolds_r"))) {
        unfolds <- get_unfolds(H)
        attr(H, "unfolds_r") <- unfolds$Re
        if (is.complex(H@data))
          attr(H, "unfolds_i") <- unfolds$Im
      }
    }
    
    H.dec <- switch(
      decomp,
      HOSVD = hosvd_mod(
        H,
        ranks = trunc_ranks,
        svd.method = svd.method,
        status = status,
        ...
      ),
      HOOI = tucker_mod(
        H,
        ranks = trunc_ranks,
        svd.method = svd.method,
        status = status,
        ...
      ),
      CP = cp_mod(H, num_components = max_rank, status = status, ...)
    )
  } else {
    H.dec <- s  
  }
  
  modes <- sapply(H.dec$U, dim)[1,]
  nm <- length(H.dec$U)

  if (missing(groups)) {
    groups = as.list(seq(min(modes)))
  }
  max_rank <- max(sapply(groups, max))

  if (decomp != "CP") {
    if (is.null(trunc_dims)) {
      trunc_dims <- seq(nm)
      trunc_ranks <- rep(max_rank, nm)
    } else {
      trunc_ranks <- modes
      trunc_ranks[trunc_dims] <- max_rank
    }
  }
    
  
  rec <- list()
  if (is.null(names(groups)))
    group.names <- paste0("F", seq_along(groups))
  else
    group.names <- names(groups)
  
  if (rec.from.decomp) {
    N <- sum(modes) - nm + 1
    attr(H.dec, "w") <- 1 / as.vector(table(tens(seq(N), modes[-nm])@data)) / N
    attr(H.dec, "N") <- N
    attr(H.dec, "modes") <- modes
    
    attr(H.dec, "us_fft") <- lapply(seq(nm), \(j) {
      mvfft(.mvpad0(H.dec$U[[j]][, seq(trunc_ranks[j]), drop = FALSE], N))
    })
  }
  
  for (i in seq_along(groups)) {
    if (decomp == "CP") {
      H.rec <- cp_reconstruct_part(H.dec, groups[[i]])
      rec[[i]] <- reconstruct.group(H.rec, kind = "CP")
      if (!missing(cp_span)) {
        if (is.character(cp_span)) {
          cp_span <- tolower(cp_span)
          cp_span <- match.arg(cp_span)
          rec[[i]] <- c(rec[[i]], rep(switch(
            cp_span,
            mean = mean(rec[[i]]),
            last = tail(rec[[i]], 1)
          ), length(s) - length(rec[[i]])))
        } else if (is.numeric(cp_span)) {
          rec[[i]] <- c(rec[[i]], rep(cp_span, length(s) - length(rec[[i]])))
        } else {
          cp_span <- match.fun(cp_span)
          rec[[i]] <- c(rec[[i]], rep(cp_span(rec[[i]]), length(s) - length(rec[[i]])))
        }
      }
    } else {
      group <- rep(list(groups[[i]]), nm)
      group[-trunc_dims] <- lapply(modes[-trunc_dims], seq)
      
      if (rec.from.decomp) {
        rec[[i]] <- reconstruct.group.from.decomp(H.dec, group = group)
      } else {
        H.rec <- ttl(
          do.call(`[`, c(list(H.dec$Z), group, list(drop = FALSE))), 
          lapply(seq(nm), \(i) as.matrix(H.dec$U[[i]][, group[[i]]])), 
          seq(nm))
        
        rec[[i]] <- reconstruct.group(H.rec, kind = "SSA")
      }
    }
  }
  
  names(rec) <- group.names
  rec
}

# Complex multi-channel ESPRIT
cmesprit <- function(s, L, groups, qrtol = 1e-07) {
  H <- apply(s, 2, Rssa::hankel, L = L, simplify = FALSE) |>
    Reduce(cbind, x = _)
  
  max_rank <- max(sapply(groups, max))
  H.svd <- svd(H, nu = max_rank, nv = 0)
  estimates <- list()
  for (i in seq(groups)) {
    U <- H.svd$u[, groups[[i]], drop = FALSE]
    Z <- qr.solve(U[-nrow(U), ], U[-1, ], qrtol)
    poles <- as.complex(eigen(Z, only.values = TRUE)$values)
    estimates[[i]] <- list(
      poles = as.complex(poles),
      rates = Re(log(poles)),
      frequencies = Im(log(poles)) / 2 / pi
    )
  }
  estimates
}

# Complex multi-channel SSA (decomposition and reconstruction)
cmssa_reconstruct <- function(s, L, groups) {
  stopifnot(is.matrix(s) && is.list(groups))
  Q <- ncol(s)
  K <- nrow(s) - L + 1
  H <- apply(s, 2, Rssa::hankel, L = L, simplify = FALSE) |>
    Reduce(cbind, x = _)
  
  max_rank <- max(sapply(groups, max))
  H.dec <- svd(H, nu = max_rank, nv = max_rank)
  rec <- list()
  
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    rec.mat <- H.dec$u[, group] %*% diag(H.dec$d[group], nrow = length(group)) %*% Conj(t(H.dec$v[, group]))
    rec.tens <- rTensor::fold(rec.mat, 1, 2:3, modes = c(L, K, Q))
    rec[[i]] <- reconstruct.group3(rec.tens, kind = "MSSA")
  }
  rec
}

# Complex Circular White Gaussian Noise Generator
CCSWGN <- function(n, mean = 0, sd = 1) {
  rnorm(n, mean = Re(mean), sd = sd / sqrt(2)) + 1i * rnorm(n, mean = Im(mean), sd = sd / sqrt(2))
}

# HOSVD-MSSA
tens_mssa_decompose <- function(s,
                                L,
                                decomp = c("HOSVD", "HOOI", "TSVD", "CP"),
                                svd.method = c("svd", "primme"),
                                neig = NULL,
                                status = FALSE,
                                max_iter = 25,
                                tol = 1e-05,
                                cp_start = c("rand", "svd")) {
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  
  if (decomp == "CP" && is.null(neig))
    stop("For CP decomposition `neig` argument must be provided")
  else if (is.null(neig))
    neig <- 50
  
  if (decomp != "CP" && length(neig) == 1)
    neig <- c(min(neig, L), min(neig, nrow(s) - L + 1), min(neig, ncol(s)))
  
  H <- tens3(s, L, kind = "MSSA")
  H.dec <- switch(
    decomp,
    HOSVD = hosvd_mod(
      H,
      ranks = neig,
      svd.method = svd.method,
      status = status
    ),
    HOOI = tucker_mod(
      H,
      ranks = neig,
      status = status,
      max_iter = max_iter,
      tol = tol
    ),
    TSVD = tsvd_mod(H, status = status),
    CP = cp_mod(
      H,
      num_components = neig,
      status = status,
      max_iter = max_iter,
      tol = tol,
      start = cp_start
    )
  )
  
  H.dec
}

tens_mssa_reconstruct <- function(s,
                                  L,
                                  groups,
                                  groups3,
                                  decomp = c("HOSVD", "HOOI", "TSVD", "CP"),
                                  svd.method = c("svd", "primme"),
                                  status = FALSE,
                                  max_iter = 25,
                                  tol = 1e-05,
                                  cp_start = c("rand", "svd")) {
  if (!is.list(groups))
    groups <- as.list(groups)
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  
  H <- tens3(s, L, kind = "MSSA")
  max_rank <- max(sapply(groups, max))
  if (decomp == "HOSVD" || decomp == "HOOI") {
    if (!is.list(groups3))
      groups3 <- as.list(groups3)
    if (length(groups) != length(groups3))
      simpleError(paste0(
        "Lengths of groups and groups3 are not equal: ",
        length(groups),
        " != ",
        length(groups3)
      ))
    max_rank3 <- max(sapply(groups3, max))
  }
  
  H.dec <- switch(
    decomp,
    HOSVD = hosvd_mod(
      H,
      ranks = c(max_rank, max_rank, max_rank3),
      svd.method = svd.method,
      status = status
    ),
    HOOI = tucker_mod(
      H,
      ranks = c(max_rank, max_rank, max_rank3),
      status = status,
      max_iter = max_iter,
      tol = tol
    ),
    TSVD = tsvd_mod(H, status = status),
    CP = cp_mod(
      H,
      num_components = max_rank,
      status = status,
      max_iter = max_iter,
      tol = tol,
      start = cp_start
    )
  )
  
  rec <- list()
  group.names <- names(groups)
  
  if (is.null(group.names))
    group.names <- paste0("F", seq_along(groups))
  
  if (decomp == "CP") {
    l_ord <- order(H.dec$lambdas, decreasing = TRUE)
  }
  
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    if (decomp == "CP") {
      H.rec <- cp_reconstruct_part(H.dec, l_ord[group])
    } else if (decomp == "TSVD") {
      H.rec <- tsvd_reconstruct_part(H.dec, group)
    } else {
      group3 <- groups3[[i]]
      H.rec <- ttl(H.dec$Z[group, group, group3, drop = FALSE], list(
        as.matrix(H.dec$U[[1]][, group]),
        as.matrix(H.dec$U[[2]][, group]),
        as.matrix(H.dec$U[[3]][, group3])
      ), 1:3)
    }
    rec[[i]] <- reconstruct.group3(H.rec, kind = "MSSA")
  }
  
  names(rec) <- group.names
  rec
}

Dstack <- function(x,
                   D = N %/% out_len,
                   out_len = NULL) {
  N <- length(x)
  M <- N %/% D
  if (N %% D != 0)
    simpleWarning("N %% D != 0, last elements of x will not be considered")
  
  outer(1:M, 1:D, function(m, d)
    x[(m - 1) * D + d])
}

HTLSDstack <- function(x,
                       D = N %/% out_len,
                       L,
                       groups,
                       method = c("base", "hosvd", "hooi"),
                       r3,
                       est_dim,
                       out_len,
                       ...) {
  method <- match.arg(method)
  
  xmat <- Dstack(x, D)
  if (method == "base") {
    estimates <- cmesprit(xmat, L, groups, ...)
  }
  else {
    estimates <- tens_esprit(
      xmat,
      L = L,
      kind = "TMSSA",
      decomp = method,
      groups = groups,
      r3 = r3,
      est_dim = est_dim,
      ...
    )
  }
  
  for (i in seq(estimates)) {
    estimates[[i]]$freq <- estimates[[i]]$freq / D
  }
  estimates
}

unstack <- function(X) {
  as.vector(t(X))
}

SSADstack <- function(x,
                      D = N %/% out_len,
                      L,
                      groups,
                      method = c("base", "hosvd", "hooi"),
                      r3,
                      out_len,
                      ...) {
  method <- match.arg(method)
  
  xmat <- Dstack(x, D)
  if (method == "base") {
    s <- reconstruct(ssa(xmat, L, kind = "mssa", ...), groups = groups)
  }
  else {
    s <- tens_mssa_reconstruct(xmat, L, groups, list(1:r3), decomp = method, ...)
  }
  
  lapply(s, unstack)
}