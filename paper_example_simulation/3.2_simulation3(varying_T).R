# ============================================================
## Experiment A — Vary N (N=M ∈ {25,50,100,200,400}, T=5 fixed)
## Experiment B — Vary T ({2,5,10,20}, N=100 fixed)
## Experiment C — Vary K ({1, 2, 3, 4})
# ============================================================


# ============================================================
# Experiment B: Computation Time — Vary T (number of periods)
#
# DGP: identical to Part 1 (long-dyad-first, cross-dyad dependence)
# Fixed: N = 100, K = 2
# Varied: T ∈ {2, 5, 10, 20}
#   ramen  runs on N × (N+20) = 100 × 120  — all T periods jointly
#   AMEN   runs on N × N      = 100 × 100  — T independent ame() calls
#
# Comparison:
#   ramen  = fit_temporal_bipartite_als (T periods) +
#            bootstrap_als_factorize_joint_cov (all T periods)
#   AMEN   = sum of T × ame() calls (one per period, no temporal structure)
#
# NOTE: ramen handles temporal dependence natively via gamma smoothing.
#       AMEN has no temporal structure — running it T times independently
#       is the only way to process a T-period panel.
# ============================================================
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(devtools)
library(amen)

save_dir <- "xxx"
load_all(".")

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
N_fixed <- 100
M_extra <- 20        # ramen: N×(N+20);  AMEN: N×N
M_ramen_fixed <- N_fixed + M_extra   # 120
M_amen_fixed  <- N_fixed             # 100

K          <- 2
lambda_val <- 0.30
gamma_val  <- 0.70

T_grid <- c(2, 5, 10, 20)

nscan_amen <- 10000
burn_amen  <- 1000

B_boot <- 1000

n_time_reps <- 1

# ---- DGP parameters (identical to Part 1) ----
sigma_eps      <- 0.35
phi_alpha      <- 0.75;  phi_beta  <- 0.75
phi_U          <- 0.92;  phi_V     <- 0.92
sd_alpha0      <- 0.60;  sd_beta0  <- 0.60
sd_U0          <- 1.20;  sd_V0     <- 1.20
sd_alpha_innov <- 0.18;  sd_beta_innov  <- 0.18
sd_U_innov     <- 0.30;  sd_V_innov     <- 0.30
rho_row        <- 0.25;  rho_col        <- 0.25
coef_row_true  <- c(0.8, -0.6)
coef_col_true  <- c(0.5,  0.4)
coef_dyad_true <- c(0.7)

# ------------------------------------------------------------
# 2. DGP generator  (TT and M are arguments, not globals)
# ------------------------------------------------------------
generate_panel <- function(N, M, TT, seed = 42) {
  set.seed(seed)
  years   <- 2001:(2000 + TT)
  row_ids <- paste0("r", seq_len(N))
  col_ids <- paste0("c", seq_len(M))

  delta_true <- rep(0, TT);  names(delta_true) <- as.character(years)

  # ---- Latent factors ----
  alpha_list <- beta_list <- U_list <- V_list <- vector("list", TT)
  for (tt in seq_len(TT)) {
    if (tt == 1) {
      alpha_list[[1]] <- rnorm(N, sd = sd_alpha0)
      beta_list[[1]]  <- rnorm(M, sd = sd_beta0)
      U_list[[1]]     <- matrix(rnorm(N * K, sd = sd_U0), N, K)
      V_list[[1]]     <- matrix(rnorm(M * K, sd = sd_V0), M, K)
    } else {
      alpha_list[[tt]] <- phi_alpha * alpha_list[[tt-1]] + rnorm(N, sd = sd_alpha_innov)
      beta_list[[tt]]  <- phi_beta  * beta_list[[tt-1]]  + rnorm(M, sd = sd_beta_innov)
      U_list[[tt]]     <- phi_U * U_list[[tt-1]] + matrix(rnorm(N*K, sd = sd_U_innov), N, K)
      V_list[[tt]]     <- phi_V * V_list[[tt-1]] + matrix(rnorm(M*K, sd = sd_V_innov), M, K)
    }
    names(alpha_list[[tt]]) <- row_ids;  names(beta_list[[tt]])  <- col_ids
    rownames(U_list[[tt]])  <- row_ids;  rownames(V_list[[tt]])  <- col_ids
  }

  # ---- Long dyadic panel ----
  dyad_panel <- vector("list", TT)
  for (tt in seq_len(TT)) {
    yr        <- years[tt]
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
      prev <- dyad_panel[[tt-1]]
      row_cov_t  <- prev %>% distinct(node_row, x_row_1, x_row_2) %>%
        mutate(x_row_1 = 0.50 * x_row_1 + rnorm(n(), sd = 0.8),
               x_row_2 = 0.50 * x_row_2 + rnorm(n(), sd = 0.8))
      col_cov_t  <- prev %>% distinct(node_col, x_col_1, x_col_2) %>%
        mutate(x_col_1 = 0.50 * x_col_1 + rnorm(n(), sd = 0.8),
               x_col_2 = 0.50 * x_col_2 + rnorm(n(), sd = 0.8))
      dyad_cov_t <- prev %>% distinct(node_row, node_col, x_dyad_1) %>%
        mutate(x_dyad_1 = 0.40 * x_dyad_1 + rnorm(n(), sd = 0.9))
    }

    dyad_df_t <- base_grid %>%
      left_join(row_cov_t,  by = "node_row") %>%
      left_join(col_cov_t,  by = "node_col") %>%
      left_join(dyad_cov_t, by = c("node_row", "node_col")) %>%
      mutate(
        alpha_t     = alpha_list[[tt]][node_row],
        beta_t      = beta_list[[tt]][node_col],
        u1 = U_list[[tt]][node_row, 1], u2 = U_list[[tt]][node_row, 2],
        v1 = V_list[[tt]][node_col, 1], v2 = V_list[[tt]][node_col, 2],
        row_part    = coef_row_true[1] * x_row_1  + coef_row_true[2] * x_row_2,
        col_part    = coef_col_true[1] * x_col_1  + coef_col_true[2] * x_col_2,
        dyad_part   = coef_dyad_true[1] * x_dyad_1,
        year_part   = delta_true[as.character(yr)],
        latent_part = u1 * v1 + u2 * v2,
        base_signal = row_part + col_part + dyad_part + year_part +
          alpha_t + beta_t + latent_part
      )

    base_mat <- matrix(NA_real_, N, M, dimnames = list(row_ids, col_ids))
    base_mat[cbind(match(dyad_df_t$node_row, row_ids),
                   match(dyad_df_t$node_col, col_ids))] <- dyad_df_t$base_signal
    row_sum <- rowSums(base_mat);  col_sum <- colSums(base_mat)
    same_row_mean <- (row_sum - base_mat) / (M - 1)
    same_col_mean <- t((col_sum - t(base_mat)) / (N - 1))

    same_df <- expand.grid(node_row = row_ids, node_col = col_ids,
                           stringsAsFactors = FALSE) %>%
      mutate(same_row_mean = as.vector(same_row_mean),
             same_col_mean = as.vector(same_col_mean))

    dyad_df_t <- dyad_df_t %>%
      left_join(same_df, by = c("node_row", "node_col")) %>%
      mutate(signal = base_signal + rho_row * same_row_mean + rho_col * same_col_mean,
             eps    = rnorm(n(), sd = sigma_eps),
             value  = signal + eps)

    dyad_panel[[tt]] <- dyad_df_t
  }

  analysis_df <- bind_rows(dyad_panel)

  edge_panel  <- analysis_df %>% dplyr::select(year, node_row, node_col, value)
  row_cov_df  <- analysis_df %>% distinct(year, node_row, x_row_1, x_row_2)
  col_cov_df  <- analysis_df %>% distinct(year, node_col, x_col_1, x_col_2)
  dyad_cov_df <- analysis_df %>% distinct(year, node_row, node_col, x_dyad_1)

  Y_list <- Xrow_list <- Xcol_list <- vector("list", TT)
  for (tt in seq_len(TT)) {
    df_t  <- dyad_panel[[tt]]
    Y_t   <- matrix(NA_real_, N, M, dimnames = list(row_ids, col_ids))
    Y_t[cbind(match(df_t$node_row, row_ids),
              match(df_t$node_col, col_ids))] <- df_t$value
    rc <- df_t %>% distinct(node_row, x_row_1, x_row_2)
    cc <- df_t %>% distinct(node_col, x_col_1, x_col_2)
    Xr <- as.matrix(rc[match(row_ids, rc$node_row), c("x_row_1","x_row_2")])
    Xc <- as.matrix(cc[match(col_ids, cc$node_col), c("x_col_1","x_col_2")])
    rownames(Xr) <- row_ids;  rownames(Xc) <- col_ids
    Y_list[[tt]] <- Y_t;  Xrow_list[[tt]] <- Xr;  Xcol_list[[tt]] <- Xc
  }

  list(N=N, M=M, TT=TT, years=years,
       edge_panel=edge_panel, row_cov_df=row_cov_df,
       col_cov_df=col_cov_df, dyad_cov_df=dyad_cov_df,
       Y_list=Y_list, Xrow_list=Xrow_list, Xcol_list=Xcol_list)
}

# ------------------------------------------------------------
# 3. Timing helpers
# ------------------------------------------------------------

# ramen: fit all T periods jointly + bootstrap ALL T periods
# Symmetric with AMEN: both produce uncertainty estimates for every period.
# For tt=1: gamma=0, no prev anchors
# For tt>1: gamma=gamma_val, prev period estimates used as anchors
time_ramen_boot <- function(dat, B = B_boot) {
  TT    <- dat$TT
  years <- dat$years
  N_d   <- dat$N;  M_d <- dat$M

  # Step 1: full panel point estimate (all T periods)
  t_fit <- system.time(
    fit_panel <- fit_temporal_bipartite_als(
      edge_panel            = dat$edge_panel,
      row_cov_df            = dat$row_cov_df,
      col_cov_df            = dat$col_cov_df,
      dyad_cov_df           = dat$dyad_cov_df,
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
      lambda                = lambda_val,
      gamma                 = gamma_val,
      verbose               = FALSE
    )
  )["elapsed"]

  # Step 2: bootstrap each period, accumulate time
  t_boot_total <- 0

  for (tt in seq_len(TT)) {
    yr_tt  <- as.character(years[tt])
    res_tt <- fit_panel$results[[yr_tt]]
    r_names <- rownames(res_tt$U)   # "row_r1", ...
    c_names <- names(res_tt$beta)   # "col_c1", ...

    # Prep Y, Xrow, Xcol with prefixed names for period tt
    Y_b <- dat$Y_list[[tt]]
    rownames(Y_b) <- r_names;  colnames(Y_b) <- c_names

    Xr_b <- dat$Xrow_list[[tt]];  rownames(Xr_b) <- r_names
    Xc_b <- dat$Xcol_list[[tt]];  rownames(Xc_b) <- c_names

    X_dyad_b <- array(NA_real_, dim = c(N_d, M_d, 1),
                      dimnames = list(r_names, c_names, "x_dyad_1"))
    dd <- dat$dyad_cov_df %>% filter(year == years[tt])
    ri <- match(paste0("row_", dd$node_row), r_names)
    ci <- match(paste0("col_", dd$node_col), c_names)
    ok <- !is.na(ri) & !is.na(ci)
    X_dyad_b[cbind(ri[ok], ci[ok], 1)] <- dd$x_dyad_1[ok]

    # Previous period anchors
    if (tt > 1) {
      yr_prev  <- as.character(years[tt - 1])
      res_prev <- fit_panel$results[[yr_prev]]
      U_prev_b     <- res_prev$U;      rownames(U_prev_b) <- r_names
      V_prev_b     <- res_prev$V;      rownames(V_prev_b) <- c_names
      alpha_prev_b <- res_prev$alpha;  names(alpha_prev_b) <- r_names
      beta_prev_b  <- res_prev$beta;   names(beta_prev_b)  <- c_names
      gamma_b      <- gamma_val
    } else {
      U_prev_b <- V_prev_b <- alpha_prev_b <- beta_prev_b <- NULL
      gamma_b  <- 0
    }

    t_boot_tt <- system.time(
      bootstrap_als_factorize_joint_cov(
        Y              = Y_b,
        B              = B,
        K              = K,
        lambda         = lambda_val,
        gamma          = gamma_b,
        U_prev         = U_prev_b,
        V_prev         = V_prev_b,
        alpha_prev     = alpha_prev_b,
        beta_prev      = beta_prev_b,
        U_init         = res_tt$U,
        V_init         = res_tt$V,
        X_row          = Xr_b,
        X_col          = Xc_b,
        X_dyad         = X_dyad_b,
        coef_row_init  = res_tt$coef_row,
        coef_col_init  = res_tt$coef_col,
        coef_dyad_init = res_tt$coef_dyad,
        bootstrap_type       = "wild",
        use_original_as_init = TRUE
      )
    )["elapsed"]

    t_boot_total <- t_boot_total + t_boot_tt
  }

  t_fit + t_boot_total
}

# AMEN: T independent ame() calls (one per period), total time summed
# This is the only option — ame() has no temporal structure
time_amen <- function(dat) {
  TT <- dat$TT
  total <- 0
  for (tt in seq_len(TT)) {
    t_t <- system.time(
      suppressMessages(
        amen::ame(
          Y         = dat$Y_list[[tt]],
          Xrow      = dat$Xrow_list[[tt]],
          Xcol      = dat$Xcol_list[[tt]],
          R         = K,
          family = "nrm",
          symmetric = FALSE,
          nscan     = nscan_amen,
          burn      = burn_amen,
          odens     = nscan_amen,
          plot      = FALSE,
          print     = FALSE
        )
      )
    )["elapsed"]
    total <- total + t_t
  }
  total
}

median_elapsed <- function(f, dat) {
  times <- numeric(n_time_reps)
  for (r in seq_len(n_time_reps)) times[r] <- f(dat)
  median(times)
}

# ------------------------------------------------------------
# 4. Experiment B: loop over T
# ------------------------------------------------------------
cat("\n============================================================\n")
cat(sprintf("Experiment B: Vary T  (N=%d, K=%d, B=%d, nscan_amen=%d)\n",
            N_fixed, K, B_boot, nscan_amen))
cat(sprintf("  ramen: %d\u00d7%d  |  AMEN: %d\u00d7%d (per period)\n",
            N_fixed, M_ramen_fixed, N_fixed, M_amen_fixed))
cat("============================================================\n\n")

timing_B <- data.frame(
  TT           = integer(0),
  t_ramen_boot = numeric(0),
  t_amen       = numeric(0),
  amen_ok      = logical(0)
)

for (TT_val in T_grid) {
  cat(sprintf("── T=%d ──\n", TT_val))

  dat_ramen <- generate_panel(N = N_fixed, M = M_ramen_fixed, TT = TT_val, seed = 42)
  dat_amen  <- generate_panel(N = N_fixed, M = M_amen_fixed,  TT = TT_val, seed = 42)

  # ramen + bootstrap (last period)
  t_rb <- tryCatch(
    median_elapsed(time_ramen_boot, dat_ramen),
    error = function(e) { cat("  ramen boot FAILED:", conditionMessage(e), "\n"); NA }
  )
  cat(sprintf("  ramen+boot (T=%d): %.2f sec  (B=%d)\n", TT_val, t_rb, B_boot))

  # AMEN: T × ame()
  t_a <- NA_real_;  ok <- FALSE
  tryCatch({
    t_a <- median_elapsed(time_amen, dat_amen)
    ok  <- TRUE
    cat(sprintf("  AMEN       (T=%d): %.2f sec  (%d × nscan=%d)\n",
                TT_val, t_a, TT_val, nscan_amen))
  }, error = function(e) {
    cat(sprintf("  AMEN: FAILED — %s\n", conditionMessage(e)))
  })

  timing_B <- rbind(timing_B, data.frame(
    TT = TT_val, t_ramen_boot = t_rb, t_amen = t_a, amen_ok = ok
  ))
}

timing_B <- timing_B %>%
  mutate(speedup = ifelse(amen_ok, t_amen / t_ramen_boot, NA_real_))

cat("\n── Summary ──\n")
print(timing_B %>% dplyr::select(TT, t_ramen_boot, t_amen, speedup), digits = 3)

# ------------------------------------------------------------
# 5. Visualisation
# ------------------------------------------------------------
plot_df <- timing_B %>%
  tidyr::pivot_longer(c(t_ramen_boot, t_amen),
                      names_to  = "method",
                      values_to = "elapsed") %>%
  mutate(
    method = dplyr::recode(method,
                           t_ramen_boot = sprintf("ramen (point est. + bootstrap B=%d)", B_boot),
                           t_amen       = sprintf("AMEN (%d \u00d7 nscan=%d)", 1L, nscan_amen)
    ),
    method = factor(method, levels = c(
      sprintf("ramen (point est. + bootstrap B=%d)", B_boot),
      sprintf("AMEN (%d \u00d7 nscan=%d)", 1L, nscan_amen)
    ))
  ) %>%
  filter(!is.na(elapsed))

p_loglog <- ggplot(plot_df, aes(x = TT, y = elapsed,
                                color = method, shape = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_x_log10(breaks = T_grid, labels = paste0("T=", T_grid)) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.01)) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  scale_shape_manual(values = c(16, 17)) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(
    x        = "Number of periods T  (log scale)",
    y        = "Elapsed time in seconds  (log scale)",
    color    = NULL, shape = NULL,
    title    = sprintf("Computation time vs number of periods  (N=%d, K=%d)", N_fixed, K),
    subtitle = sprintf("ramen: %d\u00d7%d (joint, \u00d7T bootstrap), B=%d  |  AMEN: %d\u00d7%d (\u00d7T), nscan=%d",
                       N_fixed, M_ramen_fixed, B_boot,
                       N_fixed, M_amen_fixed, nscan_amen)
  )
print(p_loglog)

p_bar <- ggplot(plot_df, aes(x = factor(TT), y = elapsed, fill = method)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x = "Number of periods T", y = "Elapsed time (seconds)", fill = NULL,
       title = "Computation time — bar chart")
print(p_bar)



# ------------------------------------------------------------
# 6. Save results
# ------------------------------------------------------------
save_dir  <- "xxx"
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

save_path <- file.path(save_dir, sprintf("part3b.RData", timestamp))
save(timing_B, T_grid, N_fixed, K, lambda_val, gamma_val,
     nscan_amen, burn_amen, B_boot, n_time_reps,
     file = save_path)
cat(sprintf("\nResults saved to:\n  %s\n", save_path))

for (pp in list(list(p_loglog, "loglog"), list(p_bar, "bar"))) {
  path_p <- file.path(save_dir,
                      paste0("part3b_", pp[[2]],".png"))
  ggsave(path_p, plot = pp[[1]], width = 8, height = 5)
  cat(sprintf("Plot saved to:\n  %s\n", path_p))
}

