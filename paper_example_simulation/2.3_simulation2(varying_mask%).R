# ============================================================
# Simulation Part 2c: OOS Validation — Vary mask_frac
#
# Goal: Show that OOS performance is stable across different
#       fractions of held-out dyads in the last period.
#
# Fixed:  T = 5,  gamma = 1.0
# Varied: mask_frac ∈ {0.1, 0.3, 0.5, 0.8, 0.9}
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
gamma_opt  <- 0.3
TT_fixed   <- 5
n_rep      <- 1000
N          <- 100
M          <- 120
K          <- 2
row_ids    <- paste0("r", seq_len(N))
col_ids    <- paste0("c", seq_len(M))

mask_frac_grid <- c(0.1, 0.3, 0.5, 0.8, 0.9)

years_v   <- 2001:(2001 + TT_fixed - 1)
last_yr_v <- years_v[TT_fixed]
delta_v   <- setNames(seq(-0.20, 0.20, length.out = TT_fixed),
                      as.character(years_v))

# ---- DGP parameters (same as Part 2b) ----
sigma_eps      <- 0.35
phi_alpha      <- 0.75;  phi_beta  <- 0.75
phi_U          <- 0.95;  phi_V     <- 0.95
sd_alpha0      <- 0.60;  sd_beta0  <- 0.60
sd_U0          <- 2.5;   sd_V0     <- 2.5
sd_alpha_innov <- 0.18;  sd_beta_innov  <- 0.18
sd_U_innov     <- 0.4;   sd_V_innov     <- 0.4
rho_row        <- 0.25;  rho_col        <- 0.25
coef_row_true  <- c(0.8, -0.6)
coef_col_true  <- c(0.5,  0.4)
coef_dyad_true <- c(0.7)

# ------------------------------------------------------------
# 2. Helper functions  (identical to Part 2b)
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
  N <- nrow(base_mat);  M <- ncol(base_mat)
  row_sum <- rowSums(base_mat);  col_sum <- colSums(base_mat)
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
  N_t <- length(pkg_row_ids);  M_t <- length(pkg_col_ids)
  fitted_mat  <- matrix(res_t$alpha, N_t, M_t, byrow = FALSE) +
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
      row_cov_df[match(orig_rows, row_cov_df$node_row), names(cr), drop = FALSE])
    fitted_mat <- fitted_mat +
      matrix(as.vector(row_cov_mat %*% cr), N_t, M_t, byrow = FALSE)
  }
  if (length(cc) > 0) {
    col_cov_df  <- yr_df %>% distinct(node_col, x_col_1, x_col_2)
    col_cov_mat <- as.matrix(
      col_cov_df[match(orig_cols, col_cov_df$node_col), names(cc), drop = FALSE])
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
oos_rmse_arr <- matrix(NA_real_, n_rep, length(mask_frac_grid),
                       dimnames = list(NULL, as.character(mask_frac_grid)))
oos_cor_arr  <- matrix(NA_real_, n_rep, length(mask_frac_grid),
                       dimnames = list(NULL, as.character(mask_frac_grid)))
success_arr  <- matrix(FALSE,    n_rep, length(mask_frac_grid),
                       dimnames = list(NULL, as.character(mask_frac_grid)))

# ------------------------------------------------------------
# 4. Monte Carlo loop
# Generate data ONCE per rep, then apply each mask_frac
# ------------------------------------------------------------
for (rr in seq_len(n_rep)) {
  cat(sprintf("\nReplication %d / %d\n", rr, n_rep))

  # ---- Generate DGP (T = TT_fixed periods) ----
  alpha_list <- beta_list <- U_list <- V_list <- vector("list", TT_fixed)
  for (tt in seq_len(TT_fixed)) {
    if (tt == 1) {
      alpha_t <- rnorm(N, sd = sd_alpha0);  beta_t <- rnorm(M, sd = sd_beta0)
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

  dyad_panel <- vector("list", TT_fixed)
  for (tt in seq_len(TT_fixed)) {
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
        v1 = V_list[[tt]][node_col, 1],  v2 = V_list[[tt]][node_col, 2],
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

  analysis_df <- bind_rows(dyad_panel)
  covars      <- make_covariate_dfs_from_long(analysis_df)
  last_df     <- analysis_df %>% filter(year == last_yr_v)

  # ---- Loop over mask_frac values (reuse same DGP draw) ----
  for (mf_idx in seq_along(mask_frac_grid)) {
    mf <- mask_frac_grid[mf_idx]

    test_idx   <- sample(seq_len(nrow(last_df)),
                         size = floor(mf * nrow(last_df)))
    test_dyads <- last_df[test_idx, ]
    train_df   <- analysis_df %>%
      anti_join(test_dyads %>% dplyr::select(year, node_row, node_col),
                by = c("year", "node_row", "node_col"))

    edge_panel_train <- train_df %>% dplyr::select(year, node_row, node_col, value)

    pkg_try <- try(
      fit_temporal_bipartite_als(
        edge_panel            = edge_panel_train,
        row_cov_df            = covars$row_cov_df,
        col_cov_df            = covars$col_cov_df,
        dyad_cov_df           = covars$dyad_cov_df,
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
      res_last   <- pkg_try$results[[as.character(last_yr_v)]]
      coef_last  <- pkg_try$coef_df %>% filter(period == last_yr_v)
      fitted_mat <- predict_last_period(res_last, coef_last, analysis_df, last_yr_v)

      test_row_int <- paste0("row_", test_dyads$node_row)
      test_col_int <- paste0("col_", test_dyads$node_col)
      r_idx <- match(test_row_int, rownames(fitted_mat))
      c_idx <- match(test_col_int, colnames(fitted_mat))
      ok    <- !is.na(r_idx) & !is.na(c_idx)

      if (sum(ok) > 10) {
        pred_vals <- fitted_mat[cbind(r_idx[ok], c_idx[ok])]
        true_vals <- test_dyads$value[ok]
        oos_rmse_arr[rr, mf_idx] <- sqrt(mean((pred_vals - true_vals)^2, na.rm = TRUE))
        oos_cor_arr[rr,  mf_idx] <- cor(pred_vals, true_vals, use = "complete.obs")
        success_arr[rr,  mf_idx] <- TRUE
      }
    }

    cat(sprintf("  mask=%.1f | RMSE=%.4f\n", mf,
                mean(oos_rmse_arr[, mf_idx], na.rm = TRUE)))
  }
}

# ------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------
n_success_vec <- colSums(success_arr)

mask_summary <- data.frame(
  mask_frac = mask_frac_grid,
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
cat(sprintf("OOS Validation Summary: Vary mask_frac  (T=%d, gamma=%.1f)\n",
            TT_fixed, gamma_opt))
cat("==============================\n\n")
print(mask_summary %>%
        dplyr::select(mask_frac, mean_rmse, sd_rmse, mean_cor, n_success),
      digits = 4)

# ------------------------------------------------------------
# 6. Visualization
# ------------------------------------------------------------
p_out <- ggplot(mask_summary, aes(x = mask_frac, y = mean_rmse)) +
  ## geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, fill = "steelblue") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 3.5, color = "steelblue") +
  scale_x_continuous(breaks = mask_frac_grid,
                     labels = scales::percent(mask_frac_grid)) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Fraction of last-period dyads held out",
    y = sprintf("Out-of-sample RMSE  (\u03b3 = %.1f)", gamma_opt),
    title = sprintf("OOS RMSE vs held-out fraction  (T=%d, N=%d)", TT_fixed, N)
  )
print(p_out)

# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------
save_dir  <- "xxxx"
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save_path <- file.path(save_dir, "part2c_varyMask.RData")
save(oos_rmse_arr, oos_cor_arr, success_arr,
     mask_summary, mask_frac_grid,
     gamma_opt, TT_fixed, n_rep,
     file = save_path)
cat(sprintf("\nResults saved to:\n  %s\n", save_path))

plot_path <- file.path(save_dir, "part2c_varyMask.png")
ggsave(plot_path, plot = p_out, width = 7, height = 5)
cat(sprintf("Plot saved to:\n  %s\n", plot_path))



