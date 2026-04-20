# ============================================================
# Simulation Part 2b: OOS Validation — Vary T (fixed gamma)
#
# Goal: Show that more periods (longer panel) increases the
#       benefit of temporal smoothing (gamma_opt vs gamma=0).
#
# Run Part 2a first to determine gamma_opt, then set it below.
# ============================================================
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(devtools)

save_dir <- "xxx"
load_all(".")

set.seed(123)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
# !! Set gamma_opt from Part 2a results !!
gamma_opt  <- 0.3

n_rep_T    <- 1000
N          <- 100
M          <- 120
K          <- 2
row_ids    <- paste0("r", seq_len(N))
col_ids    <- paste0("c", seq_len(M))

mask_frac  <- 0.3
T_grid     <- c(2, 3, 5, 7, 10)

# ---- DGP parameters (same as Part 2a) ----
sigma_eps      <- 0.35
phi_alpha      <- 0.75
phi_beta       <- 0.75
phi_U          <- 0.95
phi_V          <- 0.95
sd_alpha0      <- 0.60
sd_beta0       <- 0.60
sd_U0          <- 2.5
sd_V0          <- 2.5
sd_alpha_innov <- 0.18
sd_beta_innov  <- 0.18
sd_U_innov     <- 0.4
sd_V_innov     <- 0.4
rho_row        <- 0.25
rho_col        <- 0.25
coef_row_true  <- c(0.8, -0.6)
coef_col_true  <- c(0.5,  0.4)
coef_dyad_true <- c(0.7)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
long_to_matrix <- function(df, value_col, row_ids, col_ids) {
  out   <- matrix(NA_real_, nrow = length(row_ids), ncol = length(col_ids),
                  dimnames = list(row_ids, col_ids))
  r_idx <- match(df$node_row, row_ids)
  c_idx <- match(df$node_col, col_ids)
  ok    <- !is.na(r_idx) & !is.na(c_idx)
  out[cbind(r_idx[ok], c_idx[ok])] <- df[[value_col]][ok]
  out
}

compute_cross_dyad_means <- function(base_mat) {
  N       <- nrow(base_mat)
  M       <- ncol(base_mat)
  row_sum <- rowSums(base_mat)
  col_sum <- colSums(base_mat)
  same_row_mean <- matrix(0, N, M)
  same_col_mean <- matrix(0, N, M)
  for (j in seq_len(M)) same_row_mean[, j] <- (row_sum - base_mat[, j]) / (M - 1)
  for (i in seq_len(N)) same_col_mean[i, ]  <- (col_sum - base_mat[i, ]) / (N - 1)
  list(same_row_mean = same_row_mean, same_col_mean = same_col_mean)
}

make_covariate_dfs_from_long <- function(dyad_df) {
  list(
    row_cov_df  = dyad_df %>% distinct(year, node_row, x_row_1, x_row_2),
    col_cov_df  = dyad_df %>% distinct(year, node_col, x_col_1, x_col_2),
    dyad_cov_df = dyad_df %>% distinct(year, node_row, node_col, x_dyad_1)
  )
}

predict_last_period <- function(res_t, coef_df_t, analysis_df, yr) {
  pkg_row_ids <- names(res_t$alpha)
  pkg_col_ids <- names(res_t$beta)
  N_t <- length(pkg_row_ids)
  M_t <- length(pkg_col_ids)

  fitted_mat <- matrix(res_t$alpha, N_t, M_t, byrow = FALSE) +
    matrix(res_t$beta,  N_t, M_t, byrow = TRUE)  +
    res_t$U %*% t(res_t$V)
  rownames(fitted_mat) <- pkg_row_ids
  colnames(fitted_mat) <- pkg_col_ids

  orig_rows <- sub("^row_", "", pkg_row_ids)
  orig_cols <- sub("^col_", "", pkg_col_ids)
  yr_df <- analysis_df %>% filter(year == yr)

  cr <- setNames(coef_df_t$coefficient[coef_df_t$type == "row"],
                 coef_df_t$covariate[coef_df_t$type == "row"])
  cc <- setNames(coef_df_t$coefficient[coef_df_t$type == "col"],
                 coef_df_t$covariate[coef_df_t$type == "col"])
  cd <- setNames(coef_df_t$coefficient[coef_df_t$type == "dyad"],
                 coef_df_t$covariate[coef_df_t$type == "dyad"])

  if (length(cr) > 0) {
    row_cov_df  <- yr_df %>% distinct(node_row, x_row_1, x_row_2)
    row_cov_mat <- as.matrix(
      row_cov_df[match(orig_rows, row_cov_df$node_row), names(cr), drop = FALSE]
    )
    fitted_mat <- fitted_mat +
      matrix(as.vector(row_cov_mat %*% cr), N_t, M_t, byrow = FALSE)
  }
  if (length(cc) > 0) {
    col_cov_df  <- yr_df %>% distinct(node_col, x_col_1, x_col_2)
    col_cov_mat <- as.matrix(
      col_cov_df[match(orig_cols, col_cov_df$node_col), names(cc), drop = FALSE]
    )
    fitted_mat <- fitted_mat +
      matrix(as.vector(col_cov_mat %*% cc), N_t, M_t, byrow = TRUE)
  }
  if (length(cd) > 0) {
    for (p in seq_along(cd)) {
      cov_name  <- names(cd)[p]
      dyad_wide <- yr_df %>%
        dplyr::select(node_row, node_col, all_of(cov_name)) %>%
        tidyr::pivot_wider(names_from = node_col, values_from = all_of(cov_name)) %>%
        tibble::column_to_rownames("node_row") %>%
        as.matrix()
      fitted_mat <- fitted_mat +
        cd[p] * dyad_wide[orig_rows, orig_cols, drop = FALSE]
    }
  }
  fitted_mat
}

# ------------------------------------------------------------
# 3. Storage
# ------------------------------------------------------------
# 2D matrix: [rep, T_val]  — only gamma_opt
oos_rmse_arr <- matrix(NA_real_, n_rep_T, length(T_grid),
                       dimnames = list(NULL, as.character(T_grid)))
oos_cor_arr  <- matrix(NA_real_, n_rep_T, length(T_grid),
                       dimnames = list(NULL, as.character(T_grid)))
success_arr  <- matrix(FALSE,    n_rep_T, length(T_grid),
                       dimnames = list(NULL, as.character(T_grid)))

# ------------------------------------------------------------
# 4. Monte Carlo loop
# ------------------------------------------------------------
for (rr in seq_len(n_rep_T)) {
  cat(sprintf("\nReplication %d / %d\n", rr, n_rep_T))

  for (tt_idx in seq_along(T_grid)) {
    TT_val    <- T_grid[tt_idx]
    years_v   <- 2001:(2001 + TT_val - 1)
    last_yr_v <- years_v[TT_val]
    # interpolate delta across TT_val periods (same range as Part 2a)
    delta_v   <- setNames(seq(-0.20, 0.20, length.out = TT_val),
                          as.character(years_v))

    # ---- Generate DGP for TT_val periods ----
    alpha_list <- vector("list", TT_val)
    beta_list  <- vector("list", TT_val)
    U_list     <- vector("list", TT_val)
    V_list     <- vector("list", TT_val)

    for (tt in seq_len(TT_val)) {
      if (tt == 1) {
        alpha_t <- rnorm(N, sd = sd_alpha0)
        beta_t  <- rnorm(M, sd = sd_beta0)
        U_t     <- matrix(rnorm(N * K, sd = sd_U0), N, K)
        V_t     <- matrix(rnorm(M * K, sd = sd_V0), M, K)
      } else {
        alpha_t <- phi_alpha * alpha_list[[tt-1]] + rnorm(N, sd = sd_alpha_innov)
        beta_t  <- phi_beta  * beta_list[[tt-1]]  + rnorm(M, sd = sd_beta_innov)
        U_t     <- phi_U * U_list[[tt-1]] + matrix(rnorm(N * K, sd = sd_U_innov), N, K)
        V_t     <- phi_V * V_list[[tt-1]] + matrix(rnorm(M * K, sd = sd_V_innov), M, K)
      }
      names(alpha_t) <- row_ids;  names(beta_t)  <- col_ids
      rownames(U_t)  <- row_ids;  rownames(V_t)  <- col_ids
      alpha_list[[tt]] <- alpha_t;  beta_list[[tt]] <- beta_t
      U_list[[tt]]     <- U_t;      V_list[[tt]]    <- V_t
    }

    dyad_panel <- vector("list", TT_val)
    for (tt in seq_len(TT_val)) {
      yr        <- years_v[tt]
      base_grid <- expand.grid(node_row = row_ids, node_col = col_ids,
                               stringsAsFactors = FALSE) %>% mutate(year = yr)
      if (tt == 1) {
        row_cov_t  <- data.frame(node_row = row_ids,
                                 x_row_1 = rnorm(N), x_row_2 = rnorm(N),
                                 stringsAsFactors = FALSE)
        col_cov_t  <- data.frame(node_col = col_ids,
                                 x_col_1 = rnorm(M), x_col_2 = rnorm(M),
                                 stringsAsFactors = FALSE)
        dyad_cov_t <- base_grid %>% transmute(node_row, node_col, x_dyad_1 = rnorm(n()))
      } else {
        prev_rows  <- dyad_panel[[tt-1]] %>% distinct(node_row, x_row_1, x_row_2)
        prev_cols  <- dyad_panel[[tt-1]] %>% distinct(node_col, x_col_1, x_col_2)
        prev_dyad  <- dyad_panel[[tt-1]] %>% distinct(node_row, node_col, x_dyad_1)
        row_cov_t  <- prev_rows %>%
          mutate(x_row_1 = 0.50 * x_row_1 + rnorm(n(), sd = 0.8),
                 x_row_2 = 0.50 * x_row_2 + rnorm(n(), sd = 0.8))
        col_cov_t  <- prev_cols %>%
          mutate(x_col_1 = 0.50 * x_col_1 + rnorm(n(), sd = 0.8),
                 x_col_2 = 0.50 * x_col_2 + rnorm(n(), sd = 0.8))
        dyad_cov_t <- prev_dyad %>%
          mutate(x_dyad_1 = 0.40 * x_dyad_1 + rnorm(n(), sd = 0.9))
      }

      dyad_df_t <- base_grid %>%
        left_join(row_cov_t,  by = "node_row") %>%
        left_join(col_cov_t,  by = "node_col") %>%
        left_join(dyad_cov_t, by = c("node_row", "node_col")) %>%
        mutate(
          alpha_t_val = alpha_list[[tt]][node_row],
          beta_t_val  = beta_list[[tt]][node_col],
          u1 = U_list[[tt]][node_row, 1],  u2 = U_list[[tt]][node_row, 2],
          v1 = V_list[[tt]][node_col, 1],  v2 = V_list[[tt]][node_col, 2]
        ) %>%
        mutate(
          row_part    = coef_row_true[1] * x_row_1 + coef_row_true[2] * x_row_2,
          col_part    = coef_col_true[1] * x_col_1 + coef_col_true[2] * x_col_2,
          dyad_part   = coef_dyad_true[1] * x_dyad_1,
          year_part   = delta_v[as.character(yr)],
          latent_part = u1 * v1 + u2 * v2,
          base_signal = row_part + col_part + dyad_part + year_part +
            alpha_t_val + beta_t_val + latent_part
        )

      base_mat  <- long_to_matrix(dyad_df_t, "base_signal", row_ids, col_ids)
      dep_means <- compute_cross_dyad_means(base_mat)
      same_df   <- expand.grid(node_row = row_ids, node_col = col_ids,
                               stringsAsFactors = FALSE) %>%
        mutate(same_row_mean = as.vector(dep_means$same_row_mean),
               same_col_mean = as.vector(dep_means$same_col_mean))
      dyad_df_t <- dyad_df_t %>%
        left_join(same_df, by = c("node_row", "node_col")) %>%
        mutate(signal = base_signal + rho_row * same_row_mean + rho_col * same_col_mean,
               eps    = rnorm(n(), sd = sigma_eps),
               value  = signal + eps)
      dyad_panel[[tt]] <- dyad_df_t
    }
    analysis_df_v <- bind_rows(dyad_panel)

    # ---- Test mask for last period ----
    last_df_v  <- analysis_df_v %>% filter(year == last_yr_v)
    test_idx_v <- sample(seq_len(nrow(last_df_v)),
                         size = floor(mask_frac * nrow(last_df_v)))
    test_dyads_v <- last_df_v[test_idx_v, ]

    train_df_v <- analysis_df_v %>%
      anti_join(test_dyads_v %>% dplyr::select(year, node_row, node_col),
                by = c("year", "node_row", "node_col"))

    edge_panel_train_v <- train_df_v %>% dplyr::select(year, node_row, node_col, value)
    covars_v           <- make_covariate_dfs_from_long(analysis_df_v)

    # ---- Fit: gamma_opt only ----
    pkg_try <- try(
      fit_temporal_bipartite_als(
        edge_panel            = edge_panel_train_v,
        row_cov_df            = covars_v$row_cov_df,
        col_cov_df            = covars_v$col_cov_df,
        dyad_cov_df           = covars_v$dyad_cov_df,
        row_covar_names       = c("x_row_1", "x_row_2"),
        col_covar_names       = c("x_col_1", "x_col_2"),
        dyad_covar_names      = "x_dyad_1",
        time_col              = "year",
        row_col               = "node_row",
        col_col               = "node_col",
        value_col             = "value",
        transform             = identity,
        demean_row_covariates = TRUE,
        demean_col_covariates = TRUE,
        K      = K,
        lambda = 0.30,
        gamma  = gamma_opt,
        verbose = FALSE
      ),
      silent = TRUE
    )

    if (!inherits(pkg_try, "try-error")) {
      last_key_v  <- as.character(last_yr_v)
      res_last_v  <- pkg_try$results[[last_key_v]]
      coef_last_v <- pkg_try$coef_df %>% filter(period == last_yr_v)
      fitted_mat_v <- predict_last_period(res_last_v, coef_last_v,
                                          analysis_df_v, last_yr_v)

      test_row_int_v <- paste0("row_", test_dyads_v$node_row)
      test_col_int_v <- paste0("col_", test_dyads_v$node_col)
      r_idx_v <- match(test_row_int_v, rownames(fitted_mat_v))
      c_idx_v <- match(test_col_int_v, colnames(fitted_mat_v))
      ok_v    <- !is.na(r_idx_v) & !is.na(c_idx_v)

      if (sum(ok_v) > 10) {
        pred_vals_v <- fitted_mat_v[cbind(r_idx_v[ok_v], c_idx_v[ok_v])]
        true_vals_v <- test_dyads_v$value[ok_v]
        oos_rmse_arr[rr, tt_idx] <- sqrt(mean((pred_vals_v - true_vals_v)^2, na.rm = TRUE))
        oos_cor_arr[rr,  tt_idx] <- cor(pred_vals_v, true_vals_v, use = "complete.obs")
        success_arr[rr,  tt_idx] <- TRUE
      }
    }
    cat(sprintf("  T=%2d | RMSE[γ=%.1f]=%.4f\n",
                TT_val, gamma_opt,
                mean(oos_rmse_arr[, tt_idx], na.rm = TRUE)))
  }
}

# ------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------
n_success_vec <- colSums(success_arr)

T_summary <- data.frame(
  T_periods = T_grid,
  mean_rmse = colMeans(oos_rmse_arr, na.rm = TRUE),
  sd_rmse   = apply(oos_rmse_arr, 2, sd, na.rm = TRUE),
  mean_cor  = colMeans(oos_cor_arr,  na.rm = TRUE),
  sd_cor    = apply(oos_cor_arr,  2, sd, na.rm = TRUE),
  n_success = n_success_vec
) %>%
  mutate(
    se_rmse = sd_rmse / sqrt(n_success),
    ci_lo   = mean_rmse - 1.96 * se_rmse,
    ci_hi   = mean_rmse + 1.96 * se_rmse
  )

cat("\n==============================\n")
cat(sprintf("OOS Validation Summary: Vary T  (gamma = %.1f)\n", gamma_opt))
cat("==============================\n\n")
print(T_summary %>% dplyr::select(T_periods, mean_rmse, sd_rmse, mean_cor, n_success),
      digits = 4)

# ------------------------------------------------------------
# 6. Visualization
# ------------------------------------------------------------
ggplot(T_summary, aes(x = T_periods, y = mean_rmse)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, fill = "steelblue") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 3.5, color = "steelblue") +
  scale_x_continuous(breaks = T_grid) +
  theme_minimal(base_size = 16) +
  labs(x = "Number of time periods (T)",
       y = sprintf("Out-of-sample RMSE  (γ = %.1f)", gamma_opt))


# ------------------------------------------------------------
# 7. Save results to disk
# ------------------------------------------------------------
save_dir <- "xxx"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 7a. Save all arrays + summary table as .RData
save_path <- file.path(save_dir, "part2b_varyT.RData")
save(oos_rmse_arr, oos_cor_arr, success_arr,
     T_summary, T_grid, gamma_opt, n_rep_T, mask_frac,
     file = save_path)
cat(sprintf("\nResults saved to:\n  %s\n", save_path))


p_out <- ggplot(T_summary, aes(x = T_periods, y = mean_rmse)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, fill = "steelblue") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 3.5, color = "steelblue") +
  scale_x_continuous(breaks = T_grid) +
  theme_minimal(base_size = 16) +
  labs(x = "Number of time periods (T)",
       y = sprintf("Out-of-sample RMSE  (\u03b3 = %.1f)", gamma_opt))
p_out


plot_path <- file.path(save_dir, "part2b_varyT.png")
ggsave(plot_path, plot = p_out, width = 7, height = 5)
cat(sprintf("Plot saved to:\n  %s\n", plot_path))

