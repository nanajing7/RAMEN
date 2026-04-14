library(devtools)
load_all()   # 开发阶段用

set.seed(123)

# Dimensions
N <- 20   # rows
M <- 15   # columns
K <- 2    # latent dimension

# Covariates
X_row <- matrix(rnorm(N * 2), N, 2)
X_col <- matrix(rnorm(M * 2), M, 2)
X_dyad <- array(rnorm(N * M * 1), dim = c(N, M, 1))

# True parameters
beta_row <- c(0.5, -0.3)
beta_col <- c(0.4, 0.2)
beta_dyad <- 0.6

alpha <- rnorm(N, sd = 0.5)
beta <- rnorm(M, sd = 0.5)

U <- matrix(rnorm(N * K), N, K)
V <- matrix(rnorm(M * K), M, K)

# Generate mean
mu <-
  X_row %*% beta_row %*% t(rep(1, M)) +
  rep(1, N) %*% t(X_col %*% beta_col) +
  beta_dyad * X_dyad[, , 1] +
  outer(alpha, rep(1, M)) +
  outer(rep(1, N), beta) +
  U %*% t(V)

# Add noise
Y <- mu + matrix(rnorm(N * M, sd = 0.5), N, M)


fit <- als_factorize_joint_cov(
  Y = Y,
  K = 2,
  lambda = 0.1,
  gamma = 1,
  X_row = X_row,
  X_col = X_col,
  X_dyad = X_dyad)



boot_out <- bootstrap_als_factorize_joint_cov(
  Y = Y,
  B = 500,
  K = 2,
  lambda = 0.1,
  gamma = 1,
  X_row = X_row,
  X_col = X_col,
  X_dyad = X_dyad,
  seed = 123
)



print(boot_out)
summary(boot_out)


coef_summary_bootstrap_als(boot_out)
latent_mean_bootstrap_als(boot_out)
latent_sd_bootstrap_als(boot_out)
latent_sign_prob_bootstrap_als(boot_out)
bootstrap_diagnostics_als(boot_out)


# Example: latent uncertainty (standard deviation)
latent_sd <- latent_sd_bootstrap_als(boot_out)

# Higher values indicate greater uncertainty in latent interactions
summary(as.vector(latent_sd))
