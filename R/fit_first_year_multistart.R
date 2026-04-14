
#' Multi-start fit for the first period only
#'
#' Uses one SVD start plus random perturbation starts and keeps the
#' solution with the lowest objective.
#'
#' @param Y Outcome matrix.
#' @param K Latent dimension.
#' @param lambda Ridge penalty.
#' @param X_row,X_col,X_dyad Optional covariates.
#' @param max_iter Maximum iterations.
#' @param tol Relative tolerance.
#' @param n_starts Number of starts.
#' @param random_sd SD of random perturbation.
#' @param verbose Logical.
#'
#' @return Best fit object.
fit_first_period_multistart <- function(Y,
                                        K = 2,
                                        lambda = 0.1,
                                        X_row = NULL,
                                        X_col = NULL,
                                        X_dyad = NULL,
                                        max_iter = 1000,
                                        tol = 1e-5,
                                        n_starts = 5,
                                        random_sd = 0.1,
                                        verbose = FALSE) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  N <- nrow(Y)
  M <- ncol(Y)

  if (n_starts < 1) stop("n_starts must be >= 1.", call. = FALSE)

  best_fit <- NULL
  best_obj <- Inf

  init_svd <- svd_init(Y, K = K)

  fit <- als_factorize_joint_cov(
    Y = Y,
    K = K,
    lambda = lambda,
    gamma = 0,
    U_init = init_svd$U,
    V_init = init_svd$V,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad,
    coef_row_prev = NULL,
    coef_col_prev = NULL,
    coef_dyad_prev = NULL,
    coef_row_init = NULL,
    coef_col_init = NULL,
    coef_dyad_init = NULL,
    max_iter = max_iter,
    tol = tol
  )

  if (verbose) {
    cat(sprintf("  start 1/%d [svd] objective = %.6f\n", n_starts, fit$objective))
  }

  if (is.finite(fit$objective) && fit$objective < best_obj) {
    best_fit <- fit
    best_obj <- fit$objective
  }

  if (n_starts >= 2) {
    for (ss in 2:n_starts) {
      U_init <- init_svd$U + matrix(stats::rnorm(N * K, mean = 0, sd = random_sd), nrow = N, ncol = K)
      V_init <- init_svd$V + matrix(stats::rnorm(M * K, mean = 0, sd = random_sd), nrow = M, ncol = K)

      rownames(U_init) <- rownames(Y)
      rownames(V_init) <- colnames(Y)
      colnames(U_init) <- paste0("k", seq_len(K))
      colnames(V_init) <- paste0("k", seq_len(K))

      fit_try <- try(
        als_factorize_joint_cov(
          Y = Y,
          K = K,
          lambda = lambda,
          gamma = 0,
          U_init = U_init,
          V_init = V_init,
          X_row = X_row,
          X_col = X_col,
          X_dyad = X_dyad,
          coef_row_prev = NULL,
          coef_col_prev = NULL,
          coef_dyad_prev = NULL,
          coef_row_init = NULL,
          coef_col_init = NULL,
          coef_dyad_init = NULL,
          max_iter = max_iter,
          tol = tol
        ),
        silent = TRUE
      )

      if (inherits(fit_try, "try-error")) {
        if (verbose) {
          cat(sprintf("  start %d/%d [random] failed\n", ss, n_starts))
        }
        next
      }

      if (verbose) {
        cat(sprintf("  start %d/%d [random] objective = %.6f\n", ss, n_starts, fit_try$objective))
      }

      if (is.finite(fit_try$objective) && fit_try$objective < best_obj) {
        best_fit <- fit_try
        best_obj <- fit_try$objective
      }
    }
  }

  if (is.null(best_fit)) {
    stop("All first-period multi-start attempts failed.", call. = FALSE)
  }

  best_fit
}
