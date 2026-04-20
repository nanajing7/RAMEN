
#' Bootstrap inference for a fitted temporal bipartite ALS panel
#'
#' Runs \code{bootstrap_als_factorize_joint_cov} for every period in the
#' panel, automatically threading the previous-period anchors
#' (\code{U_prev}, \code{V_prev}, \code{alpha_prev}, \code{beta_prev},
#' \code{coef_row_prev}, \code{coef_col_prev}, \code{coef_dyad_prev})
#' across time.  All arguments mirror \code{fit_temporal_bipartite_als}
#' exactly so you can copy-paste your fitting call and just add the
#' bootstrap-specific parameters.
#'
#' @param fit A fitted object returned by \code{fit_temporal_bipartite_als}.
#' @param edge_panel Main edge-panel data frame (same one used to produce \code{fit}).
#' @param row_cov_df Optional row covariate data frame.
#' @param col_cov_df Optional column covariate data frame.
#' @param dyad_cov_df Optional dyadic covariate data frame.
#' @param row_covar_names Character vector of row covariate names.
#' @param col_covar_names Character vector of column covariate names.
#' @param dyad_covar_names Character vector of dyadic covariate names.
#' @param time_col Name of the time column in \code{edge_panel}.
#' @param row_col Name of the row-node column in \code{edge_panel}.
#' @param col_col Name of the column-node column in \code{edge_panel}.
#' @param value_col Name of the value column in \code{edge_panel}.
#' @param row_prefix Prefix applied to row-node IDs.
#' @param col_prefix Prefix applied to column-node IDs.
#' @param row_pad Optional integer width for left-padding row IDs.
#' @param col_pad Optional integer width for left-padding column IDs.
#' @param transform Transformation applied to matrix entries.
#' @param row_cov_time_col Time column in \code{row_cov_df}.
#' @param row_cov_id_col Row-node ID column in \code{row_cov_df}.
#' @param row_cov_prefix Prefix applied to \code{row_cov_id_col}.
#' @param col_cov_time_col Time column in \code{col_cov_df}.
#' @param col_cov_id_col Column-node ID column in \code{col_cov_df}.
#' @param col_cov_prefix Prefix applied to \code{col_cov_id_col}.
#' @param dyad_cov_time_col Time column in \code{dyad_cov_df}.
#' @param dyad_cov_row_id_col Row-node ID column in \code{dyad_cov_df}.
#' @param dyad_cov_col_id_col Column-node ID column in \code{dyad_cov_df}.
#' @param dyad_cov_row_prefix Prefix applied to dyadic row IDs.
#' @param dyad_cov_col_prefix Prefix applied to dyadic column IDs.
#' @param demean_row_covariates Logical. Apply within-node time demeaning to
#'   row covariates before bootstrap (should match the original \code{fit}).
#' @param demean_col_covariates Logical. Same for column covariates.
#' @param K Latent dimension (must match \code{fit}).
#' @param lambda Ridge penalty (must match \code{fit}).
#' @param gamma Temporal smoothing penalty (must match \code{fit}).
#' @param max_iter Maximum ALS iterations per bootstrap replication.
#' @param tol Convergence tolerance.
#' @param B Number of bootstrap replications per period.
#' @param bootstrap_type One of \code{"residual"}, \code{"wild"}, or
#'   \code{"parametric"}.  \code{"wild"} is recommended for heteroskedastic
#'   trade data.
#' @param include_alpha Logical; include row fixed effects in coefficient
#'   summaries.
#' @param include_beta Logical; include column fixed effects in coefficient
#'   summaries.
#' @param use_original_as_init Logical; if \code{TRUE} each bootstrap
#'   replication is initialised at the original period estimate (faster
#'   convergence, recommended).
#' @param conf_level Confidence level for bootstrap intervals.
#' @param seed Optional integer random seed.  Period \emph{tt} uses
#'   \code{seed + tt - 1} so results are reproducible yet vary across periods.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{period_results}}{Named list (one entry per period) of
#'     \code{bootstrap_als_factorize_joint_cov} objects.}
#'   \item{\code{coef_df}}{Long data frame with columns
#'     \code{period}, \code{term}, \code{estimate}, \code{std.error},
#'     \code{conf.low}, \code{conf.high} — stacked across all periods.}
#'   \item{\code{success_df}}{Data frame reporting \code{B_requested} and
#'     \code{B_successful} for each period.}
#'   \item{\code{params}}{List of all hyperparameters used.}
#' }
#' @export
bootstrap_temporal_bipartite_als <- function(
    fit,
    edge_panel,
    row_cov_df          = NULL,
    col_cov_df          = NULL,
    dyad_cov_df         = NULL,
    row_covar_names     = character(0),
    col_covar_names     = character(0),
    dyad_covar_names    = character(0),
    time_col            = "year",
    row_col             = "node_row",
    col_col             = "node_col",
    value_col           = "value",
    row_prefix          = "row_",
    col_prefix          = "col_",
    row_pad             = NULL,
    col_pad             = NULL,
    transform           = identity,
    row_cov_time_col    = time_col,
    row_cov_id_col      = row_col,
    row_cov_prefix      = row_prefix,
    col_cov_time_col    = time_col,
    col_cov_id_col      = col_col,
    col_cov_prefix      = col_prefix,
    dyad_cov_time_col   = time_col,
    dyad_cov_row_id_col = row_col,
    dyad_cov_col_id_col = col_col,
    dyad_cov_row_prefix = row_prefix,
    dyad_cov_col_prefix = col_prefix,
    demean_row_covariates = FALSE,
    demean_col_covariates = FALSE,
    K                   = 2,
    lambda              = 0.1,
    gamma               = 1,
    max_iter            = 1000,
    tol                 = 1e-5,
    B                   = 200,
    bootstrap_type      = c("wild", "residual", "parametric"),
    include_alpha       = TRUE,
    include_beta        = TRUE,
    use_original_as_init = TRUE,
    conf_level          = 0.95,
    seed                = NULL,
    verbose             = TRUE
) {
  bootstrap_type <- match.arg(bootstrap_type)

  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.", call. = FALSE)

  # ── 1. Rebuild panel matrices and covariate lists ────────────────────────────
  panel <- build_panel_matrices(
    edge_panel = edge_panel,
    time_col   = time_col,
    row_col    = row_col,
    col_col    = col_col,
    value_col  = value_col,
    row_prefix = row_prefix,
    col_prefix = col_prefix,
    row_pad    = row_pad,
    col_pad    = col_pad,
    transform  = transform
  )

  X_row_list <- build_row_covariates(
    row_cov_df  = row_cov_df,
    years       = panel$years,
    row_ids     = panel$row_ids,
    covar_names = row_covar_names,
    time_col    = row_cov_time_col,
    row_id_col  = row_cov_id_col,
    row_prefix  = row_cov_prefix
  )

  X_col_list <- build_col_covariates(
    col_cov_df  = col_cov_df,
    years       = panel$years,
    col_ids     = panel$col_ids,
    covar_names = col_covar_names,
    time_col    = col_cov_time_col,
    col_id_col  = col_cov_id_col,
    col_prefix  = col_cov_prefix
  )

  X_dyad_list <- build_dyad_covariates(
    dyad_cov_df = dyad_cov_df,
    years       = panel$years,
    row_ids     = panel$row_ids,
    col_ids     = panel$col_ids,
    covar_names = dyad_covar_names,
    time_col    = dyad_cov_time_col,
    row_id_col  = dyad_cov_row_id_col,
    col_id_col  = dyad_cov_col_id_col,
    row_prefix  = dyad_cov_row_prefix,
    col_prefix  = dyad_cov_col_prefix
  )

  # ── 2. Within-node time demeaning (must mirror the original fit) ─────────────
  TT <- length(panel$years)

  if (demean_row_covariates && length(row_covar_names) > 0 &&
      TT > 1 && !is.null(X_row_list[[1]])) {
    X_row_mean <- Reduce("+", X_row_list) / TT
    X_row_list <- lapply(X_row_list, function(X) {
      out <- X - X_row_mean
      rownames(out) <- rownames(X)
      colnames(out) <- colnames(X)
      out
    })
    if (verbose) cat("  [within-demeaning applied to row covariates]\n")
  }

  if (demean_col_covariates && length(col_covar_names) > 0 &&
      TT > 1 && !is.null(X_col_list[[1]])) {
    X_col_mean <- Reduce("+", X_col_list) / TT
    X_col_list <- lapply(X_col_list, function(X) {
      out <- X - X_col_mean
      rownames(out) <- rownames(X)
      colnames(out) <- colnames(X)
      out
    })
    if (verbose) cat("  [within-demeaning applied to col covariates]\n")
  }

  # ── 3. Period-by-period bootstrap ────────────────────────────────────────────
  period_results <- vector("list", TT)
  names(period_results) <- as.character(panel$years)

  for (tt in seq_len(TT)) {
    period_now <- panel$years[tt]
    if (verbose) cat(sprintf("  Bootstrap period %s (%d / %d) ...\n",
                             period_now, tt, TT))

    res_t    <- fit$results[[tt]]
    res_prev <- if (tt == 1) NULL else fit$results[[tt - 1]]

    period_seed <- if (!is.null(seed)) seed + tt - 1L else NULL

    period_results[[tt]] <- bootstrap_als_factorize_joint_cov(
      Y      = panel$Y_list[[tt]],
      B      = B,
      K      = K,
      lambda = lambda,
      gamma  = gamma,

      # Previous-period anchors (NULL for the first period)
      U_prev        = res_prev$U,
      V_prev        = res_prev$V,
      alpha_prev    = res_prev$alpha,
      beta_prev     = res_prev$beta,
      coef_row_prev = res_prev$coef_row,
      coef_col_prev = res_prev$coef_col,
      coef_dyad_prev = res_prev$coef_dyad,

      # Initialise at the original period estimate
      U_init         = res_t$U,
      V_init         = res_t$V,
      coef_row_init  = res_t$coef_row,
      coef_col_init  = res_t$coef_col,
      coef_dyad_init = res_t$coef_dyad,

      # Covariates for this period
      X_row  = X_row_list[[tt]],
      X_col  = X_col_list[[tt]],
      X_dyad = X_dyad_list[[tt]],

      max_iter             = max_iter,
      tol                  = tol,
      bootstrap_type       = bootstrap_type,
      include_alpha        = include_alpha,
      include_beta         = include_beta,
      use_original_as_init = use_original_as_init,
      conf_level           = conf_level,
      seed                 = period_seed
    )
  }

  # ── 4. Stack coefficient summaries across all periods ───────────────────────
  coef_df <- dplyr::bind_rows(
    lapply(names(period_results), function(yr) {
      s <- period_results[[yr]]$coef_summary
      if (is.null(s) || nrow(s) == 0) return(NULL)
      dplyr::mutate(s, period = yr, .before = 1)
    })
  )

  # ── 5. Success rate table ────────────────────────────────────────────────────
  success_df <- data.frame(
    period      = as.character(panel$years),
    B_requested = sapply(period_results, function(x) x$B_requested),
    B_successful = sapply(period_results, function(x) x$B_successful),
    row.names   = NULL
  )

  list(
    period_results = period_results,
    coef_df        = coef_df,
    success_df     = success_df,
    params = list(
      K                     = K,
      lambda                = lambda,
      gamma                 = gamma,
      B                     = B,
      bootstrap_type        = bootstrap_type,
      conf_level            = conf_level,
      demean_row_covariates = demean_row_covariates,
      demean_col_covariates = demean_col_covariates,
      row_covar_names       = row_covar_names,
      col_covar_names       = col_covar_names,
      dyad_covar_names      = dyad_covar_names
    )
  )
}
