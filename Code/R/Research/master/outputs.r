setwd(rprojroot::find_rstudio_root_file())

library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)
library(grDevices)

entry_to_row <- function(entry) {
  params <- entry$params
  errs <- entry$errors
  types <- names(errs)
  vals <- as.numeric(errs)
  params_df <- if (length(params) > 0) {
    as.data.frame(params, stringsAsFactors = FALSE)[rep(1, length(vals)), , drop = FALSE]
  } else {
    data.frame(matrix(nrow = length(vals), ncol = 0))
  }
  out <- data.frame(params_df, error = vals, type = types, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(out) <- NULL
  out
}

results_to_df <- function(results) {
  rows <- lapply(results, entry_to_row)
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

results_to_dflist <- function(results, method_name = NULL) {
  rows <- lapply(results, entry_to_row)
  if (length(rows) == 0) return(data.frame())
  df_all <- do.call(rbind, rows)
  if (!is.null(method_name)) df_all$method = method_name
  types <- unique(df_all$type)
  ret <- list()
  for (tp in types) {
    ret[[tp]] <- df_all |> filter(type == tp) |> select(-type)
  }
  ret
}

print_errors <- function(err_df) {
  err_df$L <- err_df$L
  err_df |> 
    mutate(zapped_error = zapsmall(error)) |> 
    group_by(num, method, type) |> 
    arrange(zapped_error, L, .by_group = TRUE)
}

summarise_errors <- function(err_df) {
  err_df$L <- err_df$L
  err_df |>
    mutate(
      zapped_error = zapsmall(error),
      ed_hossa = if_else(type == "est", est_dim_hossa, NA)
    ) |>
    group_by(num, method, type) |> summarise(
      best = min(zapped_error),
      best_L = L[which.min(zapped_error)],
      avg = mean(zapped_error),
      worst = max(zapped_error),
      hossa_ed = ed_hossa[which.min(zapped_error)],
      .groups = "drop"
    ) |> arrange(num, method, type)
}

min_errors <- function(err_df) {
  err_df |> 
    mutate(zapped_error = zapsmall(error)) |> 
    group_by(num, method, type) |> 
    filter(zapped_error == min(zapped_error)) |>
    arrange(num, method, type, .by_group = TRUE)
}

make_table <- function(df) {
  df |> dplyr::select(num, method, type, best, best_L, hossa_ed) |>
    pivot_wider(names_from = type,
                values_from = c(best, best_L, hossa_ed)) |>
    select(!hossa_ed_rec & !hossa_ed_sep) |>
    mutate(
      Метод = factor(Map(\(m) switch(
        m,
        ssa = "SSA",
        hossa = "HOSSA",
        utssa = "3D UTSSA",
        utssa_4d = "4D UTSSA"
      ), method), levels = c("SSA", "HOSSA", "3D UTSSA", "4D UTSSA")),
      Выделение = sprintf(
        "\\begin{tabular}{c}%.4f\\\\ $L = (%s)$\\end{tabular}",
        round(best_rec, 4),
        Map(\(L) paste0(L, collapse = ", "), best_L_rec)
      ),
      Разделение = sprintf(
        "\\begin{tabular}{c}%.4f\\\\ $L = (%s)$\\end{tabular}",
        round(best_sep, 4),
        Map(\(L) paste0(L, collapse = ", "), best_L_sep)
      ),
      Оценка = sprintf(
        "\\begin{tabular}{c}%.3e\\\\ $L = (%s)%s$\\end{tabular}",
        round(best_est, 8),
        Map(\(L) paste0(L, collapse = ", "), best_L_est),
        if_else(method == "hossa" & !is.na(best_est), sprintf(", d = %s", hossa_ed_est), "")
      )
    ) |> 
    select(num, Метод, Выделение, Разделение, Оценка) |>
    pivot_longer(cols = c(Выделение, Разделение, Оценка), names_to = "Задача", values_to = "RMSE") |>
    pivot_wider(names_from = Метод, values_from = RMSE) |>
    arrange(num, Задача) |>
    group_by(num) |>
    group_map(
      ~ print(
        xtable(
          .x,
          caption = sprintf("Оптимальные значения RMSE и соответствующие параметры для ряда~\\ref{enum:compcase-1d-%s}", .y[[1]]),
          align = "lccccc",
          label = sprintf("tab:compcase-1d-%s", unlist(.y[[1]]))
        ),
        file = "Research/master/tables.txt",
        append = TRUE,
        table.placement = "!ht",
        caption.placement = "top",
        include.rownames = FALSE,
        sanitize.text.function = \(c) ifelse(c == "\\begin{tabular}{c}NA\\\\ $L = ()$\\end{tabular}", "---", c),
        comment = FALSE
      )
    )
}

make_plot <- function(df) {
  df$L <- sapply(df$L, \(L) if (length(L) == 1)
    L
    else
      paste0("(", paste0(L, collapse = ', '), ")"))
  df$L <- factor(df$L, levels = unique(df$L))
  df |> 
    mutate(is_est = type == "est") |>
    mutate(type = ifelse(is_est & method == "hossa", paste0(type, ", ed = ", est_dim_hossa), type)) |>
    group_by(num, method, is_est) |>
    group_map(\(.x, .y) {
      num <- unlist(.y[, 1])
      method <- unlist(.y[, 2])
      is_est <- unlist(.y[, 3])
      if (is_est) est_str <- "_est" else est_str <- ""
      fname <- sprintf("Research/master/img/compcase-1d-%s_%s%s.pdf", num, method, est_str)
      plt <- ggplot(.x, aes(x = L, y = error)) +
        geom_line(aes(
          color = type,
          linetype = type,
          group = type
        )) +
        theme_minimal() +
        theme(axis.text.x = element_text(
          angle = ifelse(method == "ssa", 0, 90),
          vjust = 0.1,
          hjust = 1,
        ), text = element_text(size=15)) +
        guides(color = guide_legend("Задача"),
               linetype = guide_legend("Задача")) +
        ylab("RMSE")
      cairo_pdf(fname, width = 14, height = 3.9)
      print(plt)
      dev.off()
    })
}

BASE_DIR <- "Research/master/snapshots/"
result_fnames <- dir(BASE_DIR)
nums <- unique(substr(result_fnames, 1, 5))

if (!exists("df_all")) {
  df_catlist <- list()
  summary_tables <- list()
  plot_tables <- list()
  
  for (num in nums) {
    fnames <- paste0(BASE_DIR, sort(grep(num, result_fnames, value = TRUE)))
    tmpenv <- new.env()
    for (fname in fnames) {
      load(fname, envir = tmpenv)
      method <- sub(".*\\/\\d+_\\d+_\\d+_(.*)\\.RData", "\\1", fname)
      with(tmpenv, {
        df <- results_to_df(results)
        df$method <- method
        df$num <- num
        df$L <- switch(EXPR = method,
                       ssa = as.list(df$L_ssa),
                       df$L_hossa)
        df_catlist[[paste0(num, method)]] <<- df
      })
    }
  }
  df_all <- do.call(bind_rows, df_catlist)
}
# if (file.exists("Research/master/tables.txt")) file.remove("Research/master/tables.txt")
# summarise_errors(df_all) |> make_table()
df_all |> make_plot()
# plts <- df_all |> make_plot()
# print(plts[[1]])
