pred_ssa_rec <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  rec <- reconstruct(ssa(x, L = L_ssa), groups = groups_ssa)
  Reduce("+", rec)
}

pred_ssa_sep <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  rec <- reconstruct(ssa(x, L = L_ssa), groups = groups_ssa)
  Reduce(cbind, rec)
}

pred_ssa_est <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  est <- parestimate(ssa(x, L = L_ssa), groups = groups_ssa)
  c(sort(est$freq), sort(est$rates))
}

pred_ssa_several <- function(rec = FALSE, sep = FALSE, est = FALSE) {
  stopifnot(any(rec, sep, est))
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  x.ssa <- ssa(x, L = L_ssa)
  result <- list()
  if (rec || sep) recon <- reconstruct(x.ssa, groups = unlist(groups_ssa, recursive = FALSE))
  if (rec) result$rec <- Reduce("+", recon)
  if (sep) result$sep <- Reduce(cbind, recon)
  if (est) {
    estim <- parestimate(x.ssa, groups = unlist(groups_ssa, recursive = FALSE))
    estim_all <- Reduce(\(l, r) list(frequencies = c(l$frequencies, r$frequencies), 
                                     rates = c(l$rates, r$rates)), estim)
    result$est <- c(sort(estim_all$freq), sort(estim_all$rates))
  }
  result
}


pred_hossa_rec <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  
  rec <- tens_ssa_reconstruct(
    x,
    unlist(L_hossa),
    groups = unlist(groups_hossa, recursive = FALSE),
    svd.method = "primme",
    rm.repeated = use.hutd,
    trunc_dims = unlist(td_hossa),
    status = FALSE
  )
  
  Reduce("+", rec)
}

pred_hossa_sep <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  
  rec <- tens_ssa_reconstruct(
    x,
    unlist(L_hossa),
    groups = groups_hossa,
    svd.method = "primme",
    rm.repeated = use.hutd,
    trunc_dims = unlist(td_hossa),
    status = FALSE
  )
  
  Reduce(cbind, rec)
}

pred_nest_hossa_sep <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  neig <- max(sapply(groups_nested_hossa, max))
  
  dec <- tens_ssa_decompose(
    x,
    unlist(L_hossa),
    neig = rep(neig, 3),
    svd.method = "primme",
    rm.repeated = TRUE,
    status = FALSE
  )
  
  rec <- tens_ssa_reconstruct(dec$est,
                              groups = groups_nested_hossa,
                              svd.method = "svd",
                              status = FALSE)
  
  Reduce(cbind, rec)
}

pred_hossa_est <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  
  est <- tens_esprit(
    x,
    unlist(L_hossa),
    groups = groups_hossa,
    decomp = "HOSVD",
    svd.method = "primme",
    rm.repeated = use.hutd,
    est_dim = est_dim_hossa,
    status = FALSE
  )[[1]]
  
  c(sort(est$freq), sort(est$rates))
}

pred_hossa_several <- function(rec = FALSE,
                               sep = FALSE,
                               sep_nest = FALSE,
                               est = FALSE) {
  stopifnot(any(rec, sep, sep_nest, est))
  if (is.complex(signal)) {
    noise <- CCSWGN(N, 0, sd)
  } else {
    noise <- rnorm(N, 0, sd)
  }
  x <- signal + noise
  Lvec <- unlist(L_hossa)
  neig <- rep(max(sapply(unlist(groups_hossa, recursive = FALSE), max)), length(Lvec) + 1)
  result <- list()
  
  td <- unlist(td_hossa)
  neig[-td] <- c(Lvec, length(x) - sum(Lvec) + 2)[-td]
  
  if (max(neig) >= 55) .svd.method <- "svd"
  else .svd.method <- "primme"
  
  x.dec <- tens_ssa_decompose(
    x,
    L = Lvec,
    neig = neig,
    svd.method = .svd.method,
    rm.repeated = use.hutd
  )
  
  if (rec)
    result$rec <- Reduce('+',
                         tens_ssa_reconstruct(
                           x.dec,
                           groups = unlist(groups_hossa, recursive = FALSE),
                           trunc_dims = td
                         ))
  if (sep)
    result$sep <- Reduce(cbind,
                         tens_ssa_reconstruct(
                           x.dec,
                           groups = unlist(groups_hossa_sep, recursive = FALSE),
                           trunc_dims = td
                         ))
  if (est) {
    estim <- tens_esprit(x.dec,
                         groups = unlist(groups_hossa, recursive = FALSE),
                         est_dim = est_dim_hossa)[[1]]
    result$est <- c(sort(estim$freq), sort(estim$rates))
  }
  if (sep_nest) {
    x.rec <- ttl(x.dec$Z, x.dec$U, ms = seq(x.dec$Z@num_modes))
    result$sep_nest <- Reduce(
      cbind,
      tens_ssa_reconstruct(x.rec, groups = groups_nested_hossa, svd.method = "primme")
    )
  }
  result
}

pred_mssa_rec <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N * P, 0, sd)
  } else {
    noise <- rnorm(N * P, 0, sd)
  }
  x <- signal + noise
  rec <- cmssa_reconstruct(x, L = L_mssa, groups = groups_mssa)
  Reduce("+", rec)
}

pred_mssa_sep <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N * P, 0, sd)
  } else {
    noise <- rnorm(N * P, 0, sd)
  }
  x <- signal + noise
  rec <- cmssa_reconstruct(x, L = L_mssa, groups = groups_mssa)
  Reduce(\(m1, m2) abind(m1, m2, along = 3), rec)
}

pred_mssa_est <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N * P, 0, sd)
  } else {
    noise <- rnorm(N * P, 0, sd)
  }
  x <- signal + noise
  est <- cmesprit(x, L = L_mssa, groups = groups_mssa)[[1]]
  c(sort(est$freq), sort(est$rates))
}

pred_homssa_rec <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N * P, 0, sd)
  } else {
    noise <- rnorm(N * P, 0, sd)
  }
  x <- signal + noise
  
  rec <- tens_mssa_reconstruct(
    x,
    L_homssa,
    groups = groups_homssa,
    groups3 = groups3_homssa,
    decomp = "HOSVD",
    svd.method = "primme",
    status = FALSE
  )
  
  Reduce("+", rec)
}

pred_homssa_sep <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N * P, 0, sd)
  } else {
    noise <- rnorm(N * P, 0, sd)
  }
  x <- signal + noise
  
  rec <- tens_mssa_reconstruct(
    x,
    L_homssa,
    groups = groups_homssa,
    groups3 = groups3_homssa,
    decomp = "HOSVD",
    svd.method = "primme",
    status = FALSE
  )
  
  Reduce(\(m1, m2) abind(m1, m2, along = 3), rec)
}

pred_homssa_est <- function() {
  if (is.complex(signal)) {
    noise <- CCSWGN(N * P, 0, sd)
  } else {
    noise <- rnorm(N * P, 0, sd)
  }
  x <- signal + noise
  
  est <- tens_esprit(
    x,
    L_homssa,
    groups = groups_homssa,
    kind = "HO-MSSA",
    decomp = "HOSVD",
    svd.method = "primme",
    rm.repeated = use.hutd,
    est_dim = est_dim_homssa,
    r3 = max.r3,
    status = FALSE
  )[[1]]
  
  c(sort(est$freq), sort(est$rates))
}

rmse_several <- function(true, pred) {
  ret <- list()
  if (!is.null(pred$rec)) ret$rec <- rmse_ts_nd(true$rec, pred$rec)
  if (!is.null(pred$sep)) ret$sep <- rmse_ts_nd(true$sep, pred$sep)
  if (!is.null(pred$est)) ret$est <- rmmse_params_mat(true$est, pred$est)
  if (!is.null(pred$sep_nest)) ret$sep_nest <- rmse_ts_nd(true$sep, pred$sep_nest)
  ret
}

calc_errors <- function(R, pred_func, err_func, true, progbar, seed = 5) {
  if (Sys.getenv("RSTUDIO") == "1" ||
      .Platform$OS.type == "windows")
    prl_plan <- multisession
  else
    prl_plan <- multicore
  with(plan(prl_plan), {
    batch_recs <- function(R) {
      future_replicate(R, {
        progbar()
        pred_func()
      }, simplify = "array", future.seed = seed)
    }
    if (is.numeric(sd) && sd == 0)
      recs <- batch_recs(2) # don't need R because no noise present. Set to 2, because with 1 will be simplified
    else
      recs <- batch_recs(R)
  })
  res <- list()
  if (!is.null(names(recs))) {
    nm <- names(recs)
    if ("rec" %in% nm)
      res$rec <- Reduce(cbind, recs[which(nm == "rec")])
    if ("sep" %in% nm)
      res$sep <- Reduce(\(m1, m2) abind(m1, m2, along = 3), recs[which(nm == "sep")])
    if ("sep_nest" %in% nm)
      res$sep_nest <- Reduce(\(m1, m2) abind(m1, m2, along = 3), recs[which(nm == "sep_nest")])
    if ("est" %in% nm)
      res$est <- Reduce(cbind, recs[which(nm == "est")])
  } else if (!is.null(dimnames(recs))){
    dmn <- dimnames(recs)[[1]]
    if ("rec" %in% dmn)
      res$rec <- Reduce(cbind, recs["rec",])
    if ("sep" %in% dmn)
      res$sep <- Reduce(\(m1, m2) abind(m1, m2, along = 3), recs["sep",])
    if ("sep_nest" %in% dmn)
      res$sep_nest <- Reduce(\(m1, m2) abind(m1, m2, along = 3), recs["sep_nest",])
    if ("est" %in% dmn)
      res$est <- Reduce(cbind, recs["est",])
  }
  list(error = err_func(true, res),
       results = res)
}

# params_grid_list: list with the parameter names as names and vectors of the parameters values as values
calc_errors_by_params <- function(params_grid_list, R, comp_func, err_func, true, seed = 5, save_snap_to = NULL) {
  params_grid <- expand.grid(params_grid_list)
  results <- list()
  q <- progressor(nrow(params_grid) * R)
  for (i in seq(nrow(params_grid))) {
    print(params_grid[i, ])
    list2env(params_grid[i, ], envir = .GlobalEnv)
    res <- calc_errors(R, comp_func, err_func, true, q, seed)
    results[i] <- list(list(params = params_grid[i, ], errors = res$error, results = res$results))
  }
  if (!is.null(save_snap_to)) save(list = ls(all.names = TRUE), file = save_snap_to)
  return(results)
}