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

if (!exists("reconstruct_group_t3", mode = "function")) {
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

# Embedding series into tensor
tens3 <- function(s, L, kind = c("SSA", "MSSA", "CP")) {
  kind <- toupper(kind)
  kind <- match.arg(kind)
  
  if (kind == "SSA") {
    stopifnot(length(L) == 2)
    N <- length(s)
    I <- L[1]
    L <- L[2]
    K <- N - I - L + 2
    X <- outer(1:I, 1:L, `+`) |>
      outer(1:K, function(il, j)
        s[il + j - 2]) |> as.tensor()
    
    h12 <- hankel(s, I + L - 1)
    h3 <- hankel(s, I + K - 1)
    sR12 <- Re(h12)
    sR3 <- Re(h3)
    unfolds_r <- list()
    
    unfolds_r[[1]] <- Rssa::new.hbhmat(sR12, c(I, 1))
    unfolds_r[[2]] <- Rssa::new.hbhmat(sR12, c(L, 1))
    unfolds_r[[3]] <- Rssa::new.hbhmat(sR3, c(K, 1))
    
    attr(X, "unfolds_r") <- unfolds_r
    
    if (is.complex(s)) {
      sI12 <- Im(h12)
      sI3 <- Im(h3)
      unfolds_i <- list()
      
      unfolds_i[[1]] <- Rssa::new.hbhmat(sI12, c(I, 1))
      unfolds_i[[2]] <- Rssa::new.hbhmat(sI12, c(L, 1))
      unfolds_i[[3]] <- Rssa::new.hbhmat(sI3, c(K, 1))
      
      attr(X, "unfolds_i") <- unfolds_i
    }
  }
  else if (kind == "CP") {
    stopifnot(length(L) == 2)
    N <- length(s)
    l <- L[1]
    J <- L[2]
    
    I <- N %/% l
    K <- l - J + 1
    X_mat <- matrix(s[1:(I * l)], nrow = I, byrow = TRUE) |>
      apply(1, function(v)
        outer(1:J, 1:K, function(k, j)
          v[k + j - 1])) |>
      as.matrix()
    dim(X_mat) <- c(J * K, I)
    X <- rTensor::fold(X_mat, 2:3, 1, modes = c(I, J, K))
  }
  else if (kind == "MSSA") {
    N <- nrow(s)
    K <- N - L + 1
    Q <- ncol(s)
    X <- apply(s, 2, Rssa::hankel, L = L, simplify = FALSE) |>
      Reduce(cbind, x = _) |>
      rTensor::fold(1, 2:3, modes = c(L, K, Q))
    
    sR <- Re(s)
    unfolds_r <- list()
    
    unfolds_r[[1]] <- Rssa::new.hbhmat(sR, c(L, 1))
    unfolds_r[[2]] <- Rssa::new.hbhmat(sR, c(K, 1))
    
    mulR <- function(v) as.numeric(t(Re(s)) %*% v)
    tmulR <- function(v) as.numeric(Re(s) %*% v)
    unfolds_r[[3]] <- extmat(mulR, tmulR, nrow = Q, ncol = N)
    
    attr(X, "unfolds_r") <- unfolds_r
    
    if (is.complex(s)) {
      sI <- Im(s)
      unfolds_i <- list()
      
      unfolds_i[[1]] <- Rssa::new.hbhmat(sI, c(L, 1))
      unfolds_i[[2]] <- Rssa::new.hbhmat(sI, c(K, 1))
      
      mulI <- function(v) as.numeric(t(Im(s)) %*% v)
      tmulI <- function(v) as.numeric(Im(s) %*% v)
      unfolds_i[[3]] <- extmat(mulI, tmulI, nrow = Q, ncol = N)
      
      attr(X, "unfolds_i") <- unfolds_i
    }
  }
  return(X)
}

.USE.R.RECONSTRUCT <- FALSE
reconstruct.group3 <- function(X.tens, kind = c("SSA", "MSSA", "CP")) {
  stopifnot(is(X.tens, "Tensor"))
  X <- X.tens@data
  kind <- toupper(kind)
  kind <- match.arg(kind)
  if (kind == "SSA") {
    if (!.USE.R.RECONSTRUCT)
      s <- reconstruct_group_t3(X)
    else {
      I <- length(X[, 1, 1])
      L <- length(X[1, , 1])
      K <- length(X[1, 1, ])
      s <- vector(mode = "numeric", length = I + L + K - 2)
      for (C in 3:(I + L + K)) {
        sum <- 0
        count <- 0
        for (i in 1:(C - 2)) {
          for (l in 1:(C - 1 - i)) {
            if (i <= I && l <= L && C - i - l <= K) {
              sum <- sum + X[i, l, C - i - l]
              count <- count + 1
            }
          }
        }
        s[C - 2] <- sum / count
      }
    }
  } else if (kind == "CP") {
    s <- Reduce(c, lapply(seq(dim(X)[1]), function(i)
      hankel(as.matrix(X[i, , ]))))
  } else if (kind == "MSSA") {
    s <- Reduce(cbind, apply(X, 3, Rssa::hankel, simplify = FALSE))
  }
  return(s)
}

# Complex norm
fnorm_complex <- function(x) {
  sqrt(sum(Re(x * Conj(x))))
}

setMethod("fnorm_complex", "Tensor", function(x) {
  fnorm_complex(x@data)
})

# Calculate tensor eigenvalues
get_tensor_eigenvalues <- function(X, dims = seq(num_dims)) {
  num_dims = length(dim(X))
  
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
                      status = TRUE)
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
  if (identical(svd.method, "primme") && !requireNamespace("PRIMME", quietly = TRUE))
      stop("PRIMME package is required for SVD method `primme'")
    
  if (status)
    pb <- txtProgressBar(min = 0, max = num_modes, style = 3)
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    temp_mat <- rs_unfold(tnsr, m = m)@data
    
    if (identical(svd.method, "svd")) {
      U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
    } else if (identical(svd.method, "primme")) {
      R <- attr(tnsr, "unfolds_r")[[m]]
      
      Ilist <- attr(tnsr, "unfolds_i")
      if (!is.null(Ilist))
        I <- Ilist[[m]]
      else
        I <- NULL
      
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
      } else {
        A <- function(x, trans) {
          rX <- Re(x); iX <- Im(x)
          if (identical(trans, "c")) {
            (tmatmul(R, rX) + tmatmul(I, iX)) +
              1i*(tmatmul(R, iX) - tmatmul(I, rX))
          } else {
            (matmul(R, rX) - matmul(I, iX)) +
              1i*(matmul(R, iX) + matmul(I, rX))
          }
        }
      }
      
      U_list[[m]] <- PRIMME::svds(A, NSvals = ranks[m], m = nrow(R), n = ncol(R), isreal = FALSE)$u
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
  est <- ttl(Z, U_list, ms = 1:num_modes)
  resid <- fnorm_complex(est - tnsr)
  list(
    Z = Z,
    U = U_list,
    est = est,
    fnorm_resid = resid
  )
}

tucker_mod <- function(tnsr,
                       ranks = NULL,
                       svd.method = c("svd", "primme"),
                       max_iter = 25,
                       tol = 1e-05,
                       status = TRUE)
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
  if (identical(svd.method, "primme") && !requireNamespace("PRIMME", quietly = TRUE))
    stop("PRIMME package is required for SVD method `primme'")
  
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    temp_mat <- rs_unfold(tnsr, m = m)@data
    if (identical(svd.method, "svd")) {
      U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
    } else if (identical(svd.method, "primme")) {
      R <- attr(tnsr, "unfolds_r")[[m]]
      I <- attr(tnsr, "unfolds_i")[[m]]
      
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
      
      U_list[[m]] <- PRIMME::svds(A, NSvals = ranks[m], m = nrow(R), n = ncol(R), isreal = FALSE)$u
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
      X <- ttl(tnsr, lapply(U_list[-m], (\(.) Conj(t(
        .
      )))), ms = modes_seq[-m])
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
                   status = TRUE)
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
tsvd_mod <- function(tnsr, status = TRUE) 
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
  fftz <- aperm(apply(tnsr@data, MARGIN = 1:2, fft), c(2, 
                                                       3, 1))
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
    as.numeric(fft(x, inverse = TRUE))/length(x)
  }
  U <- as.tensor(aperm(apply(U_arr, MARGIN = 1:2, rTensor:::.ifft),
                       c(2, 3, 1)))
  V <- as.tensor(aperm(apply(V_arr, MARGIN = 1:2, rTensor:::.ifft),
                       c(2, 3, 1)))
  S <- as.tensor(aperm(apply(S_arr, MARGIN = 1:2, rTensor:::.ifft),
                       c(2, 3, 1)))

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
                        kind = c("HO-SSA", "HO-MSSA"),
                        decomp = c("HOOI", "HOSVD"),
                        svd.method = c("svd", "primme"),
                        est_dim,
                        r3 = NULL,
                        status = TRUE,
                        qrtol = 1e-07)
{
  max_rank <- max(sapply(groups, max))
  
  kind <- toupper(kind)
  kind <- match.arg(kind)
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  
  if (identical(kind[1], "HO-SSA"))
    H <- tens3(s, L)
  else
  {
    if (is.null(r3)) {
      simpleWarning("r3 argument was not provided, setting
                    r3 as maximum across groups")
      r3 <- max_rank
    }
      H <- tens3(s, L, kind = "MSSA")
  }
  
  max_rank1 <- min(max_rank, H@modes[1])
  max_rank2 <- min(max_rank, H@modes[2])
  max_rank3 <- min(max_rank, H@modes[3])
  
  if (identical(kind[1], "HO-MSSA")) {
    max_rank3 <- min(max_rank3, r3)
    if (identical(decomp[1], "HOOI")) {
      H.decomp <- tucker_mod(H, c(max_rank1, max_rank2, max_rank3), status = status)
    } else {
      H.decomp <- hosvd_mod(H, c(max_rank1, max_rank2, max_rank3), svd.method = svd.method, status = status)
    }
  } else if (identical(decomp[1], "HOOI")) {
    H.decomp <- tucker_mod(H, c(max_rank1, max_rank2, max_rank3), status = status)
  } else {
    H.decomp <- hosvd_mod(H, c(max_rank1, max_rank2, max_rank3), svd.method = svd.method, status = status)
  }
  estimates <- list()
  for (i in seq(groups)) {
    U <- H.decomp$U[[est_dim]][, groups[[i]], drop = FALSE]
    Z <- qr.solve(U[-nrow(U), ], U[-1, ], qrtol)
    poles <- as.complex(eigen(Z, only.values = TRUE)$values)
    estimates[[i]] <- list(
      poles = poles,
      rates = Re(log(poles)),
      frequencies = Im(log(poles)) / 2 / pi
    )
  }
  estimates
}

# TSSA

tens_ssa_reconstruct <- function(s,
                                 L,
                                 groups,
                                 decomp = c("HOSVD", "HOOI", "CP"),
                                 svd.method = c("svd", "primme"),
                                 trunc_dims = NULL,
                                 status = TRUE,
                                 cp_span = c("mean", "last"),
                                 ...) {
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  if (decomp == "CP")
    H <- tens3(s, L, kind = "CP")
  else
    H <- tens3(s, L, kind = "SSA")
  if (missing(groups)) {
    groups = as.list(seq(min(dim(H))))
  }
  max_rank <- max(sapply(groups, max))
  
  if (decomp != "CP") {
    if (is.null(trunc_dims))
      trunc_ranks <- rep(max_rank, 3)
    else {
      trunc_ranks <- dim(H)
      trunc_ranks[trunc_dims] <- max_rank
    }
  }
  
  H.dec <- switch(
    decomp,
    HOSVD = hosvd_mod(H, ranks = trunc_ranks, svd.method = svd.method, status = status),
    HOOI = tucker_mod(H, ranks = trunc_ranks, status = status, ...),
    CP = cp_mod(H, num_components = max_rank, status = status, ...)
  )
  
  rec <- list()
  if (is.null(names(groups)))
    group.names <- paste0("F", seq_along(groups))
  else
    group.names <- names(groups)
  
  for (i in seq_along(groups)) {
    if (decomp == "CP") {
      H.rec <- cp_reconstruct_part(H.dec, groups[[i]])
      rec[[i]] <- reconstruct.group3(H.rec, kind = "CP")
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
      group <- rep(list(groups[[i]]), 3)
      group[-trunc_dims] <- lapply(dim(H)[-trunc_dims], seq)
      H.rec <- ttl(H.dec$Z[group[[1]], group[[2]], group[[3]], drop = FALSE], list(
        as.matrix(H.dec$U[[1]][, group[[1]]]),
        as.matrix(H.dec$U[[2]][, group[[2]]]),
        as.matrix(H.dec$U[[3]][, group[[3]]])
      ), 1:3)
      rec[[i]] <- reconstruct.group3(H.rec, kind = "SSA")
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
                                neig3 = NULL,
                                status = TRUE,
                                max_iter = 25,
                                tol = 1e-05,
                                cp_start = c("rand", "svd")) {
  decomp <- toupper(decomp)
  decomp <- match.arg(decomp)
  svd.method <- match.arg(svd.method)
  
  if (decomp == "CP" && is.null(neig))
    stop("For CP decomposition `neig` argument must be provided")
  else if (is.null(neig))
    neig = min(50, L, nrow(s) - L + 1)
  
  if (decomp != "CP" && is.null(neig3))
    neig3 = min(50, ncol(s))
  
  H <- tens3(s, L, kind = "MSSA")
  H.dec <- switch(
    decomp,
    HOSVD = hosvd_mod(H, ranks = c(neig, neig, neig3), svd.method = svd.method, status = status),
    HOOI = tucker_mod(
      H,
      ranks = c(neig, neig, neig3),
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
                                  status = TRUE,
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

Dstack <- function(x, D = N %/% out_len, out_len = NULL) {
  N <- length(x)
  M <- N %/% D
  if (N %% D != 0)
    simpleWarning("N %% D != 0, last elements of x will not be considered")
  
  outer(1:M, 1:D, function(m, d) x[(m - 1) * D + d])
}

HTLSDstack <- function(x, D = N %/% out_len, L, groups, method = c("base", "hosvd", "hooi"), r3, est_dim, out_len, ...) {  
  method <- match.arg(method)
  
  xmat <- Dstack(x, D)
  if (method == "base") {
    estimates <- cmesprit(xmat, L, groups, ...)
  }
  else {
    estimates <- tens_esprit(
      xmat,
      L = L,
      kind = "HO-MSSA",
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

SSADstack <- function(x, D = N %/% out_len, L, groups, method = c("base", "hosvd", "hooi"), r3, out_len, ...) {  
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