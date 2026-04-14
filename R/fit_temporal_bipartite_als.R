#' Fit temporal ALS for one weighted bipartite panel
#'
#' This function is package-ready and does not assume file paths or
#' project-specific object names.
#'
#' @param edge_panel Main edge-panel data frame.
#' @param row_cov_df Optional row covariate data frame.
#' @param col_cov_df Optional column covariate data frame.
#' @param dyad_cov_df Optional dyadic covariate data frame.
#' @param row_covar_names Character vector of row covariate names.
#' @param col_covar_names Character vector of column covariate names.
#' @param dyad_covar_names Character vector of dyadic covariate names.
#' @param time_col Name of the time column in `edge_panel`.
#' @param row_col Name of the row-node column in `edge_panel`.
#' @param col_col Name of the column-node column in `edge_panel`.
#' @param value_col Name of the value column in `edge_panel`.
#' @param row_prefix Prefix applied to row-node IDs.
#' @param col_prefix Prefix applied to column-node IDs.
#' @param row_pad Optional integer width for left-padding row IDs.
#' @param col_pad Optional integer width for left-padding column IDs.
#' @param transform Transformation applied to matrix entries.
#' @param row_cov_time_col Time column in `row_cov_df`.
#' @param row_cov_id_col Row-node ID column in `row_cov_df`.
#' @param row_cov_prefix Prefix applied to `row_cov_id_col`.
#' @param col_cov_time_col Time column in `col_cov_df`.
#' @param col_cov_id_col Column-node ID column in `col_cov_df`.
#' @param col_cov_prefix Prefix applied to `col_cov_id_col`.
#' @param dyad_cov_time_col Time column in `dyad_cov_df`.
#' @param dyad_cov_row_id_col Row-node ID column in `dyad_cov_df`.
#' @param dyad_cov_col_id_col Column-node ID column in `dyad_cov_df`.
#' @param dyad_cov_row_prefix Prefix applied to dyadic row IDs.
#' @param dyad_cov_col_prefix Prefix applied to dyadic column IDs.
#' @param demean_row_covariates Logical. If \code{TRUE}, subtracts each row
#'   node's time-mean from its covariates before fitting (within-node
#'   demeaning). Recommended when \code{row_covar_names} is non-empty and
#'   covariates are time-varying, to avoid collinearity with node fixed
#'   effects (\code{alpha_i}).
#' @param demean_col_covariates Logical. Same as \code{demean_row_covariates}
#'   but applied to column-node covariates.
#' @param K Latent dimension.
#' @param lambda Ridge penalty.
#' @param gamma Temporal smoothing penalty.
#' @param max_iter Maximum iterations.
#' @param tol Relative tolerance.
#' @param first_period_n_starts Number of first-period starts.
#' @param first_period_random_sd SD of first-period random perturbation.
#' @param verbose Logical.
#' @param panel_id Optional identifier stored in output.
#'
#' @return A list of period-specific results plus stacked node and coefficient data.
#' @export
fit_temporal_bipartite_als <- function(edge_panel,
                                       row_cov_df = NULL,
                                       col_cov_df = NULL,
                                       dyad_cov_df = NULL,
                                       row_covar_names = character(0),
                                       col_covar_names = character(0),
                                       dyad_covar_names = character(0),
                                       time_col = "year",
                                       row_col = "node_row",
                                       col_col = "node_col",
                                       value_col = "value",
                                       row_prefix = "row_",
                                       col_prefix = "col_",
                                       row_pad = NULL,
                                       col_pad = NULL,
                                       transform = identity,
                                       row_cov_time_col = time_col,
                                       row_cov_id_col = row_col,
                                       row_cov_prefix = row_prefix,
                                       col_cov_time_col = time_col,
                                       col_cov_id_col = col_col,
                                       col_cov_prefix = col_prefix,
                                       dyad_cov_time_col = time_col,
                                       dyad_cov_row_id_col = row_col,
                                       dyad_cov_col_id_col = col_col,
                                       dyad_cov_row_prefix = row_prefix,
                                       dyad_cov_col_prefix = col_prefix,
                                       K = 2,
                                       lambda = 0.1,
                                       gamma = 1,
                                       max_iter = 1000,
                                       tol = 1e-5,
                                       first_period_n_starts = 5,
                                       first_period_random_sd = 0.1,
                                       demean_row_covariates = FALSE,
                                       demean_col_covariates = FALSE,
                                       verbose = FALSE,
                                       panel_id = NA_character_) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  panel <- build_panel_matrices(
    edge_panel = edge_panel,
    time_col = time_col,
    row_col = row_col,
    col_col = col_col,
    value_col = value_col,
    row_prefix = row_prefix,
    col_prefix = col_prefix,
    row_pad = row_pad,
    col_pad = col_pad,
    transform = transform
  )

  X_row_list <- build_row_covariates(
    row_cov_df = row_cov_df,
    years = panel$years,
    row_ids = panel$row_ids,
    covar_names = row_covar_names,
    time_col = row_cov_time_col,
    row_id_col = row_cov_id_col,
    row_prefix = row_cov_prefix
  )

  X_col_list <- build_col_covariates(
    col_cov_df = col_cov_df,
    years = panel$years,
    col_ids = panel$col_ids,
    covar_names = col_covar_names,
    time_col = col_cov_time_col,
    col_id_col = col_cov_id_col,
    col_prefix = col_cov_prefix
  )

  X_dyad_list <- build_dyad_covariates(
    dyad_cov_df = dyad_cov_df,
    years = panel$years,
    row_ids = panel$row_ids,
    col_ids = panel$col_ids,
    covar_names = dyad_covar_names,
    time_col = dyad_cov_time_col,
    row_id_col = dyad_cov_row_id_col,
    col_id_col = dyad_cov_col_id_col,
    row_prefix = dyad_cov_row_prefix,
    col_prefix = dyad_cov_col_prefix
  )

  # ── Within-node time demeaning ─────────────────────────────────────────────
  # Subtracts each node's time-mean from its covariate values across periods.
  # This removes the component of X_row / X_col that is collinear with the
  # node fixed effects (alpha_i, beta_j), making coef_row and coef_col
  # identified from within-node time variation only.
  # Only applied when there are >= 2 periods and covariates are present.
  TT <- length(panel$years)

  if (demean_row_covariates &&
      length(row_covar_names) > 0 &&
      TT > 1 &&
      !is.null(X_row_list[[1]])) {

    X_row_mean <- Reduce("+", X_row_list) / TT
    X_row_list <- lapply(X_row_list, function(X) {
      out <- X - X_row_mean
      rownames(out) <- rownames(X)
      colnames(out) <- colnames(X)
      out
    })

    if (verbose) cat("  [within-demeaning applied to row covariates]\n")
  }

  if (demean_col_covariates &&
      length(col_covar_names) > 0 &&
      TT > 1 &&
      !is.null(X_col_list[[1]])) {

    X_col_mean <- Reduce("+", X_col_list) / TT
    X_col_list <- lapply(X_col_list, function(X) {
      out <- X - X_col_mean
      rownames(out) <- rownames(X)
      colnames(out) <- colnames(X)
      out
    })

    if (verbose) cat("  [within-demeaning applied to col covariates]\n")
  }

  out <- vector("list", length(panel$years))
  names(out) <- as.character(panel$years)

  prev_U <- NULL
  prev_V <- NULL
  prev_alpha <- NULL
  prev_beta <- NULL

  prev_coef_row <- NULL
  prev_coef_col <- NULL
  prev_coef_dyad <- NULL

  for (tt in seq_along(panel$years)) {
    period_now <- panel$years[tt]
    Y <- panel$Y_list[[tt]]
    X_row_t <- X_row_list[[tt]]
    X_col_t <- X_col_list[[tt]]
    X_dyad_t <- X_dyad_list[[tt]]

    if (tt == 1) {
      if (verbose) {
        cat(sprintf("  period %s: first-period multi-start (%d starts)\n",
                    period_now, first_period_n_starts))
      }

      fit <- fit_first_period_multistart(
        Y = Y,
        K = K,
        lambda = lambda,
        X_row = X_row_t,
        X_col = X_col_t,
        X_dyad = X_dyad_t,
        max_iter = max_iter,
        tol = tol,
        n_starts = first_period_n_starts,
        random_sd = first_period_random_sd,
        verbose = verbose
      )
    } else {
      fit <- als_factorize_joint_cov(
        Y = Y,
        K = K,
        lambda = lambda,
        gamma = gamma,
        U_prev = prev_U,
        V_prev = prev_V,
        alpha_prev = prev_alpha,
        beta_prev = prev_beta,
        U_init = prev_U,
        V_init = prev_V,
        X_row = X_row_t,
        X_col = X_col_t,
        X_dyad = X_dyad_t,
        coef_row_prev = prev_coef_row,
        coef_col_prev = prev_coef_col,
        coef_dyad_prev = prev_coef_dyad,
        coef_row_init = prev_coef_row,
        coef_col_init = prev_coef_col,
        coef_dyad_init = prev_coef_dyad,
        max_iter = max_iter,
        tol = tol
      )
    }

    prev_U <- fit$U
    prev_V <- fit$V
    prev_alpha <- fit$alpha
    prev_beta <- fit$beta
    prev_coef_row <- fit$coef_row
    prev_coef_col <- fit$coef_col
    prev_coef_dyad <- fit$coef_dyad

    row_latent_df <- as.data.frame(fit$U, stringsAsFactors = FALSE)
    names(row_latent_df) <- paste0("X", seq_len(K))

    col_latent_df <- as.data.frame(fit$V, stringsAsFactors = FALSE)
    names(col_latent_df) <- paste0("X", seq_len(K))

    row_df <- data.frame(
      panel_id = panel_id,
      period = period_now,
      node = rownames(fit$U),
      mode = "row",
      alpha = fit$alpha,
      beta = NA_real_,
      objective = fit$objective,
      iterations = fit$iterations,
      stringsAsFactors = FALSE
    ) |>
      dplyr::bind_cols(row_latent_df)

    col_df <- data.frame(
      panel_id = panel_id,
      period = period_now,
      node = rownames(fit$V),
      mode = "col",
      alpha = NA_real_,
      beta = fit$beta,
      objective = fit$objective,
      iterations = fit$iterations,
      stringsAsFactors = FALSE
    ) |>
      dplyr::bind_cols(col_latent_df)

    coef_row_df <- if (length(fit$coef_row) > 0) {
      data.frame(
        panel_id = panel_id,
        period = period_now,
        covariate = row_covar_names,
        coefficient = fit$coef_row,
        type = "row",
        stringsAsFactors = FALSE
      )
    } else NULL

    coef_col_df <- if (length(fit$coef_col) > 0) {
      data.frame(
        panel_id = panel_id,
        period = period_now,
        covariate = col_covar_names,
        coefficient = fit$coef_col,
        type = "col",
        stringsAsFactors = FALSE
      )
    } else NULL

    coef_dyad_df <- if (length(fit$coef_dyad) > 0) {
      data.frame(
        panel_id = panel_id,
        period = period_now,
        covariate = dyad_covar_names,
        coefficient = fit$coef_dyad,
        type = "dyad",
        stringsAsFactors = FALSE
      )
    } else NULL

    out[[tt]] <- list(
      period = period_now,
      U = fit$U,
      V = fit$V,
      alpha = fit$alpha,
      beta = fit$beta,
      coef_row = fit$coef_row,
      coef_col = fit$coef_col,
      coef_dyad = fit$coef_dyad,
      node_df = dplyr::bind_rows(row_df, col_df),
      coef_df = dplyr::bind_rows(coef_row_df, coef_col_df, coef_dyad_df),
      objective = fit$objective,
      iterations = fit$iterations
    )
  }

  node_df <- dplyr::bind_rows(lapply(out, function(x) x$node_df))
  coef_df <- dplyr::bind_rows(lapply(out, function(x) x$coef_df))

  list(
    panel_id = panel_id,
    results = out,
    node_df = node_df,
    coef_df = coef_df,
    params = list(
      K = K,
      lambda = lambda,
      gamma = gamma,
      max_iter = max_iter,
      tol = tol,
      row_covar_names = row_covar_names,
      col_covar_names = col_covar_names,
      dyad_covar_names = dyad_covar_names,
      time_col = time_col,
      row_col = row_col,
      col_col = col_col,
      value_col = value_col,
      row_prefix = row_prefix,
      col_prefix = col_prefix,
      row_pad = row_pad,
      col_pad = col_pad,
      first_period_n_starts = first_period_n_starts,
      first_period_random_sd = first_period_random_sd,
      demean_row_covariates = demean_row_covariates,
      demean_col_covariates = demean_col_covariates
    )
  )
}
