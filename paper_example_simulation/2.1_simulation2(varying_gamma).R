# ============================================================
# Simulation Part 2: OOS Validation — Vary Gamma
# SAME DGP as Part 1 (long-dyad-first, cross-dyad dependence)
#
# Goal:
#   Demonstrate the value of temporal dependence (gamma > 0)
#   via out-of-sample prediction on masked dyads in last period.
# ============================================================
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(devtools)
set.seed(42)
save_dir <- "xxx"
load_all(".")

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
n_rep      <- 1000
N          <- 100
M          <- 120
TT         <- 5
K          <- 2
years      <- 2001:(2001 + TT - 1)
row_ids    <- paste0("r", seq_len(N))
col_ids    <- paste0("c", seq_len(M))

# fraction of last-period dyads held out as test set
mask_frac  <- 0.3

# gamma values to evaluate (0 = no temporal smoothing = baseline)
gamma_grid <- c(0, 0.10, 0.30, 0.50, 0.70, 1.00, 1.50, 2)


# ---- DGP parameters (identical to Part 1) ----
sigma_eps      <- 0.35
#sigma_eps      <- 1
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
delta_true     <- c(-0.20, -0.05, 0.00, 0.10, 0.20)
names(delta_true) <- as.character(years)
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

# Build full (N x M) fitted matrix for last period from package output.
# Uses raw (non-demeaned) covariates, consistent with Part 1 prediction code.
# The absolute RMSE level is affected by this, but the comparison across
# gamma values is valid since all gamma share the same prediction formula.
predict_last_period <- function(res_t, coef_df_t, analysis_df, yr) {
  pkg_row_ids <- names(res_t$alpha)   # e.g. "row_r1", "row_r10", ...
  pkg_col_ids <- names(res_t$beta)    # e.g. "col_c1", "col_c10", ...
  N_t <- length(pkg_row_ids)
  M_t <- length(pkg_col_ids)

  # alpha_i + beta_j + U_i V_j'
  fitted_mat <- matrix(res_t$alpha, N_t, M_t, byrow = FALSE) +
    matrix(res_t$beta,  N_t, M_t, byrow = TRUE)  +
    res_t$U %*% t(res_t$V)
  rownames(fitted_mat) <- pkg_row_ids
  colnames(fitted_mat) <- pkg_col_ids

  # strip prefixes to get original node names for covariate lookup
  orig_rows <- sub("^row_", "", pkg_row_ids)
  orig_cols <- sub("^col_", "", pkg_col_ids)

  yr_df <- analysis_df %>% filter(year == yr)

  cr <- setNames(coef_df_t$coefficient[coef_df_t$type == "row"],
                 coef_df_t$covariate[coef_df_t$type == "row"])
  cc <- setNames(coef_df_t$coefficient[coef_df_t$type == "col"],
                 coef_df_t$covariate[coef_df_t$type == "col"])
  cd <- setNames(coef_df_t$coefficient[coef_df_t$type == "dyad"],
                 coef_df_t$covariate[coef_df_t$type == "dyad"])

  # row covariate contribution
  if (length(cr) > 0) {
    row_cov_df  <- yr_df %>% distinct(node_row, x_row_1, x_row_2)
    row_cov_mat <- as.matrix(
      row_cov_df[match(orig_rows, row_cov_df$node_row), names(cr), drop = FALSE]
    )
    fitted_mat <- fitted_mat +
      matrix(as.vector(row_cov_mat %*% cr), N_t, M_t, byrow = FALSE)
  }

  # col covariate contribution
  if (length(cc) > 0) {
    col_cov_df  <- yr_df %>% distinct(node_col, x_col_1, x_col_2)
    col_cov_mat <- as.matrix(
      col_cov_df[match(orig_cols, col_cov_df$node_col), names(cc), drop = FALSE]
    )
    fitted_mat <- fitted_mat +
      matrix(as.vector(col_cov_mat %*% cc), N_t, M_t, byrow = TRUE)
  }

  # dyad covariate contribution
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
oos_rmse_arr <- matrix(NA_real_, n_rep, length(gamma_grid),
                       dimnames = list(NULL, as.character(gamma_grid)))
oos_cor_arr  <- matrix(NA_real_, n_rep, length(gamma_grid),
                       dimnames = list(NULL, as.character(gamma_grid)))
success_arr  <- matrix(FALSE,    n_rep, length(gamma_grid),
                       dimnames = list(NULL, as.character(gamma_grid)))

# ------------------------------------------------------------
# 4. Monte Carlo loop
# ------------------------------------------------------------
for (rr in seq_len(n_rep)) {
  cat(sprintf("\nReplication %d / %d\n", rr, n_rep))

  # ---- 4.1 Generate temporal row/col effects + latent factors ----
  alpha_list <- vector("list", TT)
  beta_list  <- vector("list", TT)
  U_list     <- vector("list", TT)
  V_list     <- vector("list", TT)

  for (tt in seq_len(TT)) {
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
    names(alpha_t) <- row_ids
    names(beta_t)  <- col_ids
    rownames(U_t)  <- row_ids
    rownames(V_t)  <- col_ids
    alpha_list[[tt]] <- alpha_t
    beta_list[[tt]]  <- beta_t
    U_list[[tt]]     <- U_t
    V_list[[tt]]     <- V_t
  }

  # ---- 4.2 Generate long dyadic panel ----
  dyad_panel <- vector("list", TT)

  for (tt in seq_len(TT)) {
    yr        <- years[tt]
    base_grid <- expand.grid(node_row = row_ids, node_col = col_ids,
                             stringsAsFactors = FALSE) %>%
      mutate(year = yr)

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
        alpha_t = alpha_list[[tt]][node_row],
        beta_t  = beta_list[[tt]][node_col],
        u1 = U_list[[tt]][node_row, 1],
        u2 = U_list[[tt]][node_row, 2],
        v1 = V_list[[tt]][node_col, 1],
        v2 = V_list[[tt]][node_col, 2]
      ) %>%
      mutate(
        row_part    = coef_row_true[1] * x_row_1 + coef_row_true[2] * x_row_2,
        col_part    = coef_col_true[1] * x_col_1 + coef_col_true[2] * x_col_2,
        dyad_part   = coef_dyad_true[1] * x_dyad_1,
        year_part   = delta_true[as.character(yr)],
        latent_part = u1 * v1 + u2 * v2,
        base_signal = row_part + col_part + dyad_part + year_part +
          alpha_t + beta_t + latent_part
      )

    # inject cross-dyad dependence
    base_mat  <- long_to_matrix(dyad_df_t, "base_signal", row_ids, col_ids)
    dep_means <- compute_cross_dyad_means(base_mat)
    same_df   <- expand.grid(node_row = row_ids, node_col = col_ids,
                             stringsAsFactors = FALSE) %>%
      mutate(
        same_row_mean = as.vector(dep_means$same_row_mean),
        same_col_mean = as.vector(dep_means$same_col_mean)
      )

    dyad_df_t <- dyad_df_t %>%
      left_join(same_df, by = c("node_row", "node_col")) %>%
      mutate(
        signal = base_signal + rho_row * same_row_mean + rho_col * same_col_mean,
        eps    = rnorm(n(), sd = sigma_eps),
        value  = signal + eps
      )

    dyad_panel[[tt]] <- dyad_df_t
  }

  analysis_df <- bind_rows(dyad_panel)

  # ---- 4.3 Create test mask for last period ----
  last_yr    <- years[TT]
  last_df    <- analysis_df %>% filter(year == last_yr)
  test_idx   <- sample(seq_len(nrow(last_df)), size = floor(mask_frac * nrow(last_df)))
  test_dyads <- last_df[test_idx, ]

  # training = all periods, but last period without test dyads
  train_df <- analysis_df %>%
    anti_join(
      test_dyads %>% dplyr::select(year, node_row, node_col),
      by = c("year", "node_row", "node_col")
    )

  edge_panel_train <- train_df %>% dplyr::select(year, node_row, node_col, value)
  covars           <- make_covariate_dfs_from_long(analysis_df)

  # ---- 4.4 Fit for each gamma ----
  for (gg in seq_along(gamma_grid)) {
    gval <- gamma_grid[gg]

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
        K                     = K,
        lambda                = 0.30,
        gamma                 = gval,
        verbose               = FALSE
      ),
      silent = TRUE
    )

    if (!inherits(pkg_try, "try-error")) {
      last_key  <- as.character(last_yr)
      res_last  <- pkg_try$results[[last_key]]
      coef_last <- pkg_try$coef_df %>% filter(period == last_yr)

      fitted_mat <- predict_last_period(res_last, coef_last, analysis_df, last_yr)

      # look up test positions using internal (prefixed) node names
      test_row_int <- paste0("row_", test_dyads$node_row)
      test_col_int <- paste0("col_", test_dyads$node_col)
      r_idx <- match(test_row_int, rownames(fitted_mat))
      c_idx <- match(test_col_int, colnames(fitted_mat))
      ok    <- !is.na(r_idx) & !is.na(c_idx)

      if (sum(ok) > 10) {
        pred_vals <- fitted_mat[cbind(r_idx[ok], c_idx[ok])]
        true_vals <- test_dyads$value[ok]
        oos_rmse_arr[rr, gg] <- sqrt(mean((pred_vals - true_vals)^2, na.rm = TRUE))
        oos_cor_arr[rr, gg]  <- cor(pred_vals, true_vals, use = "complete.obs")
        success_arr[rr, gg]  <- TRUE
      }
    }
  }
}


# ------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------
n_success_vec <- colSums(success_arr)

rmse_summary <- data.frame(
  gamma     = gamma_grid,
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

# paired差值：每个replication里，每个gamma的RMSE减去gamma=0的RMSE
baseline_rmse <- oos_rmse_arr[, as.character(0)]
rmse_diff_arr <- oos_rmse_arr - baseline_rmse

rmse_summary <- rmse_summary %>%
  mutate(
    rmse_baseline  = mean_rmse[gamma == 0],
    rmse_reduction = rmse_baseline - mean_rmse,
    mean_diff      = colMeans(rmse_diff_arr, na.rm = TRUE),
    sd_diff        = apply(rmse_diff_arr, 2, sd, na.rm = TRUE),
    se_diff        = sd_diff / sqrt(n_success),
    ci_lo_paired   = mean_diff - 1.96 * se_diff,
    ci_hi_paired   = mean_diff + 1.96 * se_diff
  )

cat("\n==============================\n")
cat("OOS Validation Summary: Vary Gamma\n")
cat("==============================\n\n")
print(
  rmse_summary %>% dplyr::select(
    gamma, mean_rmse, sd_rmse, rmse_reduction,
    mean_diff, se_diff, mean_cor, sd_cor, n_success
  ),
  digits = 3
)


# ------------------------------------------------------------
# 6. Visualization
# ------------------------------------------------------------
p_out <- ggplot(rmse_summary, aes(x = gamma, y = mean_rmse)) +
  #geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.20) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  theme_minimal(base_size = 16) +
  labs(
    x     = expression(gamma ~ "(temporal smoothing penalty)"),
    y     = "Out-of-sample RMSE"
    ## title = "OOS Prediction Error vs Temporal Smoothing"
  )

p_out



# ------------------------------------------------------------
# 7. Save results to disk
# ------------------------------------------------------------
save_dir <- "xxx"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 7a. Save all arrays + summary table as .RData
save_path <- file.path(save_dir, "part2a_varyGamma.RData")
save(oos_rmse_arr, oos_cor_arr, success_arr,
     rmse_summary, rmse_diff_arr,
     gamma_grid, n_rep, mask_frac, TT,
     file = save_path)
cat(sprintf("\nResults saved to:\n  %s\n", save_path))


plot_path <- file.path(save_dir, "part2a_varyGamma.png")
ggsave(plot_path, plot = p_out, width = 7, height = 5)
cat(sprintf("Plot saved to:\n  %s\n", plot_path))

