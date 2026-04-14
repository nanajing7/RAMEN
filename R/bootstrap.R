#' Generate a bootstrap outcome matrix from fitted values and residuals
#'
#' @param fitted_mat Fitted value matrix.
#' @param resid_mat Residual matrix.
#' @param bootstrap_type One of "residual", "wild", or "parametric".
#'
#' @return A bootstrap outcome matrix.
.generate_bootstrap_Y_als <- function(
    fitted_mat,
    resid_mat,
    bootstrap_type = c("residual", "wild", "parametric")
) {
  bootstrap_type <- match.arg(bootstrap_type)

  if (!is.matrix(fitted_mat)) {
    stop("fitted_mat must be a matrix.", call. = FALSE)
  }
  if (!is.matrix(resid_mat)) {
    stop("resid_mat must be a matrix.", call. = FALSE)
  }
  if (!all(dim(fitted_mat) == dim(resid_mat))) {
    stop("fitted_mat and resid_mat must have the same dimensions.", call. = FALSE)
  }

  N <- nrow(fitted_mat)
  M <- ncol(fitted_mat)

  resid_vec <- as.vector(resid_mat)

  # Recenter residuals so bootstrap samples are centered at fitted values
  resid_vec <- resid_vec - mean(resid_vec, na.rm = TRUE)

  if (bootstrap_type == "residual") {
    e_star <- sample(resid_vec, size = length(resid_vec), replace = TRUE)
  } else if (bootstrap_type == "wild") {
    w <- sample(c(-1, 1), size = length(resid_vec), replace = TRUE)
    e_star <- resid_vec * w
  } else if (bootstrap_type == "parametric") {
    sigma_hat <- sqrt(mean(resid_vec^2, na.rm = TRUE))
    e_star <- stats::rnorm(length(resid_vec), mean = 0, sd = sigma_hat)
  }

  Y_star <- fitted_mat + matrix(e_star, nrow = N, ncol = M)

  rownames(Y_star) <- rownames(fitted_mat)
  colnames(Y_star) <- colnames(fitted_mat)

  Y_star
}


#' Summarize bootstrap coefficient draws
#'
#' @param coef_draws A B x P matrix of bootstrap coefficient draws.
#' @param coef_hat Named numeric vector of original estimates.
#' @param conf_level Confidence level.
#'
#' @return A data.frame with bootstrap summaries.
summarize_bootstrap_coefficients_als <- function(
    coef_draws,
    coef_hat,
    conf_level = 0.95
) {
  if (!is.matrix(coef_draws)) {
    stop("coef_draws must be a matrix.", call. = FALSE)
  }

  if (!is.numeric(coef_hat) || is.null(names(coef_hat))) {
    stop("coef_hat must be a named numeric vector.", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  if (ncol(coef_draws) != length(coef_hat)) {
    stop(
      "ncol(coef_draws) must equal length(coef_hat).",
      call. = FALSE
    )
  }

  if (is.null(colnames(coef_draws))) {
    stop("coef_draws must have column names.", call. = FALSE)
  }

  if (!identical(colnames(coef_draws), names(coef_hat))) {
    stop(
      "Column names of coef_draws must exactly match names(coef_hat).",
      call. = FALSE
    )
  }

  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  se_boot <- apply(coef_draws, 2, stats::sd, na.rm = TRUE)

  ci_boot <- t(apply(
    coef_draws,
    2,
    stats::quantile,
    probs = probs,
    na.rm = TRUE
  ))

  data.frame(
    term = names(coef_hat),
    estimate = as.numeric(coef_hat),
    std.error = as.numeric(se_boot),
    conf.low = ci_boot[, 1],
    conf.high = ci_boot[, 2],
    row.names = NULL
  )
}


#' Summarize bootstrap latent interaction draws
#'
#' @param latent_draws A list of bootstrap latent interaction matrices.
#' @param conf_level Confidence level.
#'
#' @return A list of matrices: mean, sd, conf.low, conf.high, sign.prob.
summarize_bootstrap_latent_als <- function(
    latent_draws,
    conf_level = 0.95
) {
  if (!is.list(latent_draws) || length(latent_draws) == 0) {
    stop("latent_draws must be a non-empty list of matrices.", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  first_mat <- latent_draws[[1]]
  if (!is.matrix(first_mat)) {
    stop("Each element of latent_draws must be a matrix.", call. = FALSE)
  }

  N <- nrow(first_mat)
  M <- ncol(first_mat)
  B <- length(latent_draws)

  first_rownames <- rownames(first_mat)
  first_colnames <- colnames(first_mat)

  arr <- array(NA_real_, dim = c(N, M, B))

  for (b in seq_len(B)) {
    Mb <- latent_draws[[b]]

    if (!is.matrix(Mb) || !all(dim(Mb) == c(N, M))) {
      stop("All latent draw matrices must have the same dimensions.", call. = FALSE)
    }

    if (!identical(rownames(Mb), first_rownames) ||
        !identical(colnames(Mb), first_colnames)) {
      stop(
        "All latent draw matrices must have identical row and column names.",
        call. = FALSE
      )
    }

    arr[, , b] <- Mb
  }

  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  mean_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  sd_mat <- apply(arr, c(1, 2), stats::sd, na.rm = TRUE)
  low_mat <- apply(arr, c(1, 2), stats::quantile, probs = probs[1], na.rm = TRUE)
  high_mat <- apply(arr, c(1, 2), stats::quantile, probs = probs[2], na.rm = TRUE)
  sign_prob_mat <- apply(arr, c(1, 2), function(x) mean(x > 0, na.rm = TRUE))

  rownames(mean_mat) <- first_rownames
  colnames(mean_mat) <- first_colnames

  rownames(sd_mat) <- first_rownames
  colnames(sd_mat) <- first_colnames

  rownames(low_mat) <- first_rownames
  colnames(low_mat) <- first_colnames

  rownames(high_mat) <- first_rownames
  colnames(high_mat) <- first_colnames

  rownames(sign_prob_mat) <- first_rownames
  colnames(sign_prob_mat) <- first_colnames

  list(
    mean = mean_mat,
    sd = sd_mat,
    conf.low = low_mat,
    conf.high = high_mat,
    sign.prob = sign_prob_mat
  )
}


#' Bootstrap inference for als_factorize_joint_cov
#'
#' @param Y Outcome matrix.
#' @param B Number of bootstrap replications.
#' @param K Latent dimension.
#' @param lambda Ridge penalty on U and V.
#' @param gamma Temporal smoothing penalty.
#' @param U_prev,V_prev Previous-period latent factors.
#' @param alpha_prev,beta_prev Previous-period additive effects.
#' @param U_init,V_init Initial values for latent factors.
#' @param X_row Optional row covariate matrix.
#' @param X_col Optional column covariate matrix.
#' @param X_dyad Optional dyad covariate array.
#' @param coef_row_prev,coef_col_prev,coef_dyad_prev Previous-period coefficients.
#' @param coef_row_init,coef_col_init,coef_dyad_init Initial coefficient values.
#' @param max_iter Maximum ALS iterations.
#' @param tol Relative tolerance.
#' @param bootstrap_type One of "residual", "wild", or "parametric".
#' @param include_alpha Logical; include row effects in coefficient summaries.
#' @param include_beta Logical; include column effects in coefficient summaries.
#' @param use_original_as_init Logical; if TRUE, bootstrap replications use original fit as initialization.
#' @param conf_level Confidence level for summaries.
#' @param seed Optional random seed.
#'
#' @return A list containing the original fit, bootstrap draws, and summaries.
#' @export
bootstrap_als_factorize_joint_cov <- function(
    Y,
    B = 200,
    K = 2,
    lambda = 0.1,
    gamma = 1,
    U_prev = NULL,
    V_prev = NULL,
    alpha_prev = NULL,
    beta_prev = NULL,
    U_init = NULL,
    V_init = NULL,
    X_row = NULL,
    X_col = NULL,
    X_dyad = NULL,
    coef_row_prev = NULL,
    coef_col_prev = NULL,
    coef_dyad_prev = NULL,
    coef_row_init = NULL,
    coef_col_init = NULL,
    coef_dyad_init = NULL,
    max_iter = 1000,
    tol = 1e-5,
    bootstrap_type = c("residual", "wild", "parametric"),
    include_alpha = TRUE,
    include_beta = TRUE,
    use_original_as_init = TRUE,
    conf_level = 0.95,
    seed = NULL
) {
  bootstrap_type <- match.arg(bootstrap_type)

  if (!is.numeric(B) || length(B) != 1L || is.na(B) ||
      B < 1 || B != as.integer(B)) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  # Original fit
  original_fit <- als_factorize_joint_cov(
    Y = Y,
    K = K,
    lambda = lambda,
    gamma = gamma,
    U_prev = U_prev,
    V_prev = V_prev,
    alpha_prev = alpha_prev,
    beta_prev = beta_prev,
    U_init = U_init,
    V_init = V_init,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad,
    coef_row_prev = coef_row_prev,
    coef_col_prev = coef_col_prev,
    coef_dyad_prev = coef_dyad_prev,
    coef_row_init = coef_row_init,
    coef_col_init = coef_col_init,
    coef_dyad_init = coef_dyad_init,
    max_iter = max_iter,
    tol = tol
  )

  fitted_mat <- fitted_als_factorize_joint_cov(
    object = original_fit,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad
  )

  resid_mat <- residuals_als_factorize_joint_cov(
    object = original_fit,
    Y = Y,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad
  )

  coef_hat <- extract_coefficients_als(
    object = original_fit,
    include_alpha = include_alpha,
    include_beta = include_beta
  )

  coef_draws_list <- vector("list", B)
  latent_draws_list <- vector("list", B)
  success <- logical(B)

  for (b in seq_len(B)) {
    Y_star <- .generate_bootstrap_Y_als(
      fitted_mat = fitted_mat,
      resid_mat = resid_mat,
      bootstrap_type = bootstrap_type
    )

    fit_b <- try(
      als_factorize_joint_cov(
        Y = Y_star,
        K = K,
        lambda = lambda,
        gamma = gamma,
        U_prev = U_prev,
        V_prev = V_prev,
        alpha_prev = alpha_prev,
        beta_prev = beta_prev,
        U_init = if (use_original_as_init) original_fit$U else U_init,
        V_init = if (use_original_as_init) original_fit$V else V_init,
        X_row = X_row,
        X_col = X_col,
        X_dyad = X_dyad,
        coef_row_prev = coef_row_prev,
        coef_col_prev = coef_col_prev,
        coef_dyad_prev = coef_dyad_prev,
        coef_row_init = if (use_original_as_init) original_fit$coef_row else coef_row_init,
        coef_col_init = if (use_original_as_init) original_fit$coef_col else coef_col_init,
        coef_dyad_init = if (use_original_as_init) original_fit$coef_dyad else coef_dyad_init,
        max_iter = max_iter,
        tol = tol
      ),
      silent = TRUE
    )

    if (!inherits(fit_b, "try-error")) {
      coef_b <- extract_coefficients_als(
        object = fit_b,
        include_alpha = include_alpha,
        include_beta = include_beta
      )

      if (length(coef_b) == length(coef_hat) &&
          identical(names(coef_b), names(coef_hat))) {
        coef_draws_list[[b]] <- coef_b
        latent_draws_list[[b]] <- extract_latent_matrix_als(fit_b)
        success[b] <- TRUE
      }
    }
  }

  coef_draws_list <- coef_draws_list[success]
  latent_draws_list <- latent_draws_list[success]

  if (length(coef_draws_list) == 0) {
    stop("All bootstrap replications failed.", call. = FALSE)
  }

  coef_draws <- matrix(
    unlist(coef_draws_list, use.names = FALSE),
    nrow = length(coef_draws_list),
    ncol = length(coef_hat),
    byrow = TRUE
  )
  colnames(coef_draws) <- names(coef_hat)
  rownames(coef_draws) <- paste0("boot_", seq_len(nrow(coef_draws)))

  coef_summary <- summarize_bootstrap_coefficients_als(
    coef_draws = coef_draws,
    coef_hat = coef_hat,
    conf_level = conf_level
  )

  latent_summary <- summarize_bootstrap_latent_als(
    latent_draws = latent_draws_list,
    conf_level = conf_level
  )

  out <- list(
    original_fit = original_fit,
    coef_hat = coef_hat,
    coef_draws = coef_draws,
    latent_draws = latent_draws_list,
    coef_summary = coef_summary,
    latent_summary = latent_summary,
    bootstrap_type = bootstrap_type,
    B_requested = B,
    B_successful = sum(success),
    success = success
  )

  class(out) <- "bootstrap_als_factorize_joint_cov"
  out
}
