if (!require(Rssa))
  install.packages("Rssa")
if (!require(rTensor))
  install.packages("rTensor")

# Single-channel series tensorisation

tens3 <- function(s, I, L, kind = c("HO-SSA", "HO-MSSA")) {
  require("rTensor")
  if (identical(kind[1], "HO-SSA")) {
    N <- length(s)
    J <- N - I - L + 2
    X <- outer(1:I, 1:L, `+`) |>
      outer(1:J, function(il, j)
        s[il + j - 2]) |> as.tensor()
  }
  else if (identical(kind[1], "HO-MSSA")) {
    if (missing(L))
      L <- I
    N <- nrow(s)
    K <- N - L + 1
    Q <- ncol(s)
    X <- apply(s, 2, Rssa::hankel, L = L, simplify = FALSE) |>
      Reduce(cbind, x = _) |>
      rTensor::fold(1, 2:3, modes = c(L, K, Q))
  }
  return(X)
}

reconstruct.group3 <- function(X.tens, kind = c("HO-SSA", "HO-MSSA")) {
  stopifnot(is(X.tens, "Tensor"))
  X <- X.tens@data
  
  if (identical(kind[1], "HO-SSA")) {
    I <- length(X[, 1, 1])
    L <- length(X[1, , 1])
    J <- length(X[1, 1, ])
    s <- vector(mode = "numeric", length = I + L + J - 2)
    for (C in 3:(I + L + J)) {
      sum <- 0
      count <- 0
      for (i in 1:(C - 2)) {
        for (l in 1:(C - 1 - i)) {
          if (i <= I && l <= L && C - i - l <= J) {
            sum <- sum + X[i, l, C - i - l]
            count <- count + 1
          }
        }
      }
      s[C - 2] <- sum / count
    }
  } else if (identical(kind[1], "HO-MSSA")) {
    s <- Reduce(cbind, apply(X, 3, Rssa::hankel, simplify = FALSE))
    # s <- reduce(apply(X, 3, Rssa::hankel, simplify = FALSE), cbind)
  }
  else {
    simpleError(paste("Unknown kind", kind))
  }
  return(s)
}

# Complex norm

fnorm_complex <- function(x) {
  sqrt(sum(abs(x) ^ 2))
}

setMethod("fnorm_complex", "Tensor", function(x) {
  sqrt(sum(abs(x@data) ^ 2))
})

# HOSVD and HOOI modifications for complex cases

hosvd_mod <- function(tnsr,
                      ranks = NULL,
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
  if (status) {
    pb <- txtProgressBar(min = 0,
                         max = num_modes,
                         style = 3)
    U_list <- vector("list", num_modes)
    for (m in 1:num_modes) {
      temp_mat <- rs_unfold(tnsr, m = m)@data
      U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
      setTxtProgressBar(pb, m)
    }
    close(pb)
  }
  else {
    U_list <- vector("list", num_modes)
    for (m in 1:num_modes) {
      temp_mat <- rs_unfold(tnsr, m = m)@data
      U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
    }
  }
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
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    temp_mat <- rs_unfold(tnsr, m = m)@data
    U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
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
  if (status) {
    pb <- txtProgressBar(min = 0,
                         max = max_iter,
                         style = 3)
    while ((curr_iter < max_iter) && (!converged)) {
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
        setTxtProgressBar(pb, max_iter)
      }
      else {
        curr_iter <- curr_iter + 1
      }
    }
    close(pb)
  }
  else {
    while ((curr_iter < max_iter) && (!converged)) {
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
      }
      else {
        curr_iter <- curr_iter + 1
      }
    }
  }
  
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

# HO-ESPRIT

tens_esprit <- function(s,
                        I,
                        L,
                        groups,
                        kind = c("HO-SSA", "HO-MSSA"),
                        decomp = c("HOOI", "HOSVD"),
                        est_dim,
                        r3 = NULL,
                        status = TRUE,
                        qrtol = 1e-07)
{
  max_rank <- max(sapply(groups, max))
  
  if (identical(kind[1], "HO-SSA"))
    H <- tens3(s, I, L)
  else
  {
    if (is.null(r3)) {
      simpleWarning("r3 argument was not provided, setting
                    r3 as maximum across groups")
      r3 <- max_rank
    }
    if (missing(L))
      L <- I
    H <- tens3(s, L, kind = "HO-MSSA")
  }
  
  if (identical(kind[1], "HO-MSSA")) {
    if (identical(decomp[1], "HOOI")) {
      H.decomp <- tucker_mod(H, c(max_rank, max_rank, r3), status = status)
    } else {
      H.decomp <- hosvd_mod(H, c(max_rank, max_rank, r3), status = status)
    }
  } else if (identical(decomp[1], "HOOI")) {
    H.decomp <- tucker_mod(H, rep(max_rank, 3), status = status)
  } else {
    H.decomp <- hosvd_mod(H, rep(max_rank, 3), status = status)
  }
  estimates <- list()
  for (i in seq(groups)) {
    U <- H.decomp$U[[est_dim]][, groups[[i]], drop = FALSE]
    Z <- qr.solve(U[-nrow(U), ], U[-1, ], qrtol)
    poles <- eigen(Z, only.values = TRUE)$values
    estimates[[i]] <- list(
      poles = poles,
      rates = Re(log(poles)),
      frequencies = Im(log(poles)) / 2 / pi
    )
  }
  estimates
}

# HO-SSA

tens_ssa_reconstruct <- function(s,
                                 I,
                                 L,
                                 groups,
                                 decomp = c("HOSVD", "HOOI"),
                                 trunc_dims = NULL,
                                 status = TRUE) {
  H <- tens3(s, I, L, kind = "HO-SSA")
  max_rank <- max(sapply(groups, max))
  if (is.null(trunc_dims))
    trunc_ranks <- rep(max_rank, 3)
  else {
    trunc_ranks <- dim(H)
    trunc_ranks[trunc_dims] <- max_rank
  }
  
  if (identical(decomp[1], "HOSVD"))
    H.dec <- hosvd_mod(H, ranks = trunc_ranks, status = status)
  else
    H.dec <- tucker_mod(H, ranks = trunc_ranks, status = status)
  
  rec <- list()
  if (is.null(names(groups)))
    group.names <- paste0("F", seq_along(groups))
  else
    group.names <- names(groups)
  
  for (i in seq_along(groups)) {
    group <- rep(list(groups[[i]]), 3)
    group[-trunc_dims] <- lapply(dim(H)[-trunc_dims], seq)
    H.rec <- ttl(H.dec$Z[group[[1]], group[[2]], group[[3]], drop = FALSE], list(
      as.matrix(H.dec$U[[1]][, group[[1]]]),
      as.matrix(H.dec$U[[2]][, group[[2]]]),
      as.matrix(H.dec$U[[3]][, group[[3]]])
    ), 1:3)
    rec[[i]] <- reconstruct.group3(H.rec)
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
    poles <- eigen(Z, only.values = TRUE)$values
    estimates[[i]] <- list(
      poles = poles,
      rates = Re(log(poles)),
      frequencies = Im(log(poles)) / 2 / pi
    )
  }
  estimates
}

# Complex Circular White Gaussian Noise Generator

CCSWGN <- function(n, mean = 0, sd = 1) {
  rnorm(n, mean = Re(mean), sd = sd / sqrt(2)) + 1i * rnorm(n, mean = Im(mean), sd = sd / sqrt(2))
}

# HOSVD-MSSA

tens_mssa_reconstruct <- function(s,
                                  L,
                                  groups,
                                  groups3,
                                  decomp = c("HOSVD", "HOOI"),
                                  status = TRUE,
                                  max_iter = 25,
                                  tol = 1e-05) {
  if (!is.list(groups))
    groups <- as.list(groups)
  if (!is.list(groups3))
    groups3 <- as.list(groups3)
  if (length(groups) != length(groups3))
    simpleError(paste0(
      "Lengths of groups and groups3 are not equal: ",
      length(groups),
      " != ",
      length(groups3)
    ))
  
  H <- tens3(s, L, kind = "HO-MSSA")
  max_rank <- max(sapply(groups, max))
  max_rank3 <- max(sapply(groups3, max))
  
  if (identical(decomp[1], "HOSVD"))
    H.dec <- hosvd_mod(H,
                       ranks = c(max_rank, max_rank, max_rank3),
                       status = status)
  else
    H.dec <- tucker_mod(
      H,
      ranks = c(max_rank, max_rank, max_rank3),
      status = status,
      max_iter = max_iter,
      tol = tol
    )
  
  rec <- list()
  if (is.null(names(groups)))
    group.names <- paste0("F", seq_along(groups))
  else
    group.names <- names(groups)
  
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    group3 <- groups3[[i]]
    H.rec <- ttl(H.dec$Z[group, group, group3, drop = FALSE], list(
      as.matrix(H.dec$U[[1]][, group]),
      as.matrix(H.dec$U[[2]][, group]),
      as.matrix(H.dec$U[[3]][, group3])
    ), 1:3)
    rec[[i]] <- reconstruct.group3(H.rec, kind = "HO-MSSA")
  }
  
  names(rec) <- group.names
  rec
}
