# Demonstration of fit_dynamic_ame() and the two bootstrap designs.
#
# A small panel is simulated from the model itself, so the estimates can be
# compared with values that are actually known. Runs in well under a minute.

library(RAMEN)

set.seed(2024)

# ── Simulate a bipartite panel from the model ──────────────────────────────
N <- 20   # senders
M <- 15   # receivers
Tn <- 6   # periods
K <- 2    # latent dimension

# One covariate on the sender side, with a coefficient that drifts over time.
z <- matrix(rnorm(N * Tn), N, Tn)
beta_true <- numeric(Tn)
beta_true[1] <- 0.8
for (t in 2:Tn) beta_true[t] <- beta_true[t - 1] + rnorm(1, sd = 0.10)

# States follow Gaussian random walks, as the model assumes.
a <- matrix(0, N, Tn); b <- matrix(0, M, Tn)
U <- V <- vector("list", Tn)
a[, 1] <- rnorm(N, sd = 1.0)
b[, 1] <- rnorm(M, sd = 0.8)
U[[1]] <- matrix(rnorm(N * K, sd = 0.9), N, K)
V[[1]] <- matrix(rnorm(M * K, sd = 0.7), M, K)
for (t in 2:Tn) {
  a[, t] <- a[, t - 1] + rnorm(N, sd = 0.30)
  b[, t] <- b[, t - 1] + rnorm(M, sd = 0.25)
  U[[t]] <- U[[t - 1]] + matrix(rnorm(N * K, sd = 0.25), N, K)
  V[[t]] <- V[[t - 1]] + matrix(rnorm(M * K, sd = 0.22), M, K)
}

sigma_eps <- 0.5
rows <- sprintf("%03d", seq_len(N))
cols <- sprintf("%03d", seq_len(M))
years <- 2015:(2014 + Tn)

# Long-format edge panel, with 10% of dyad-periods unobserved.
edge_panel <- do.call(rbind, lapply(seq_len(Tn), function(t) {
  mu <- matrix(a[, t], N, M) +
    matrix(b[, t], N, M, byrow = TRUE) +
    beta_true[t] * matrix(z[, t], N, M) +
    U[[t]] %*% t(V[[t]])
  y <- mu + matrix(rnorm(N * M, sd = sigma_eps), N, M)
  g <- expand.grid(i = seq_len(N), j = seq_len(M))
  df <- data.frame(year = years[t],
                   node_row = rows[g$i], node_col = cols[g$j],
                   value = y[cbind(g$i, g$j)],
                   stringsAsFactors = FALSE)
  df[sample(nrow(df), round(0.90 * nrow(df))), ]
}))

row_cov <- do.call(rbind, lapply(seq_len(Tn), function(t)
  data.frame(year = years[t], node_row = rows, z = z[, t],
             stringsAsFactors = FALSE)))

cat(sprintf("Simulated panel: %d senders x %d receivers, %d periods, %d edges\n\n",
            N, M, Tn, nrow(edge_panel)))


# ── Fit ────────────────────────────────────────────────────────────────────
fit <- fit_dynamic_ame(
  edge_panel = edge_panel,
  row_cov_df = row_cov,
  row_covar_names = "z",
  row_cov_prefix = "row_",
  K = K,
  n_starts = 3,
  panel_id = "demo",
  verbose = TRUE
)

print(fit)
cat("\n")
print(summary(fit))


# ── How well were the known values recovered? ──────────────────────────────
cat("\n\n===== Recovery against the simulated truth =====\n")

cat("\nCoefficient trajectory:\n")
print(round(rbind(true = beta_true,
                  estimated = as.numeric(coef_trajectory(fit))), 3))

uv_hat <- unlist(latent_trajectory(fit))
uv_true <- unlist(lapply(seq_len(Tn), function(t) U[[t]] %*% t(V[[t]])))
cat(sprintf("\ncor(U V' estimated, U V' true) = %.4f\n", cor(uv_hat, uv_true)))
cat(sprintf("cor(a estimated, a true)       = %.4f\n",
            cor(as.vector(fit$a), as.vector(a))))
cat(sprintf("\nsigma_eps^2: true %.4f, estimated %.4f\n",
            sigma_eps^2, fit$sigma_eps2))
cat("  (the observation variance is a plug-in estimate and is biased downward,\n")
cat("   because the fitted parameters absorb part of the noise)\n")


# ── Bootstrap: both designs in one call ────────────────────────────────────
# The conditional design holds the fitted trajectories fixed and regenerates
# only the observation errors, giving standard errors for this network. The
# full-model design redraws every trajectory from the priors, so each replicate
# has a different but known truth, and the estimator is judged on the
# difference between each estimate and the value that generated it.
#
# B is small here so the demo stays quick; use B = 1000 for real work.
cat("\n\n===== Bootstrap: both designs =====\n")
bt <- bootstrap_ame(fit, design = "both", B = 40, B_full = 20, seed = 1)
print(bt)

cat("\n\nCoefficient trajectory with intervals (conditional design):\n")
print(bt$conditional$summaries$beta, row.names = FALSE, digits = 3)

cat("\nSign probability of the interaction matrix, period 3 (first 5 x 5):\n")
print(round(bt$conditional$latent$sign_prob[[3]][1:5, 1:5], 2))
cat("  values near 0 or 1 indicate dyads whose sign is consistent across\n")
cat("  replicates; values near 0.5 indicate the direction is undetermined\n")

cat("\nLatent-structure recovery (full-model design):\n")
print(bt$full$latent_recovery, row.names = FALSE, digits = 3)

cat("\nEither design can be run on its own:\n")
cat("  bootstrap_ame(fit, design = \"conditional\", B = 1000)\n")
cat("  bootstrap_ame(fit, design = \"full\", B_full = 500)\n")
