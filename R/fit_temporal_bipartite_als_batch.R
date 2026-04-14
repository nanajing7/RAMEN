
#' Fit temporal ALS across a batch of weighted bipartite panels
#'
#' @param panel_list A named or unnamed list of edge-panel data frames.
#' @param row_cov_list Optional list of row covariate data frames.
#' @param col_cov_list Optional list of column covariate data frames.
#' @param dyad_cov_list Optional list of dyad covariate data frames.
#' @param row_covar_names Character vector of row covariate names.
#' @param col_covar_names Character vector of column covariate names.
#' @param dyad_covar_names Character vector of dyadic covariate names.
#' @param time_col Name of the time column in each edge panel.
#' @param row_col Name of the row-node column in each edge panel.
#' @param col_col Name of the column-node column in each edge panel.
#' @param value_col Name of the value column in each edge panel.
#' @param row_prefix Prefix applied to row-node IDs.
#' @param col_prefix Prefix applied to column-node IDs.
#' @param row_pad Optional integer width for left-padding row IDs.
#' @param col_pad Optional integer width for left-padding column IDs.
#' @param transform Transformation applied to matrix entries.
#' @param row_cov_time_col Time column in row covariates.
#' @param row_cov_id_col Row-node ID column in row covariates.
#' @param row_cov_prefix Prefix applied to row covariate IDs.
#' @param col_cov_time_col Time column in column covariates.
#' @param col_cov_id_col Column-node ID column in column covariates.
#' @param col_cov_prefix Prefix applied to column covariate IDs.
#' @param dyad_cov_time_col Time column in dyad covariates.
#' @param dyad_cov_row_id_col Row-node ID column in dyad covariates.
#' @param dyad_cov_col_id_col Column-node ID column in dyad covariates.
#' @param dyad_cov_row_prefix Prefix applied to dyad row IDs.
#' @param dyad_cov_col_prefix Prefix applied to dyad column IDs.
#' @param K Latent dimension.
#' @param lambda Ridge penalty.
#' @param gamma Temporal smoothing penalty.
#' @param max_iter Maximum iterations.
#' @param tol Relative tolerance.
#' @param first_period_n_starts Number of first-period starts.
#' @param first_period_random_sd SD of first-period random perturbation.
#' @param verbose Logical.
#'
#' @return A list of fitted panel objects.
#' @export
fit_temporal_bipartite_als_batch <- function(panel_list,
                                             row_cov_list = NULL,
                                             col_cov_list = NULL,
                                             dyad_cov_list = NULL,
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
                                             transform = identity(),
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
                                             verbose = FALSE) {
  if (!is.list(panel_list) || length(panel_list) == 0) {
    stop("panel_list must be a non-empty list.", call. = FALSE)
  }

  results <- vector("list", length(panel_list))
  names(results) <- names(panel_list)

  for (i in seq_along(panel_list)) {
    dat <- panel_list[[i]]

    if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0) {
      results[[i]] <- NULL
      next
    }

    row_cov_df_i  <- if (!is.null(row_cov_list))  row_cov_list[[i]]  else NULL
    col_cov_df_i  <- if (!is.null(col_cov_list))  col_cov_list[[i]]  else NULL
    dyad_cov_df_i <- if (!is.null(dyad_cov_list)) dyad_cov_list[[i]] else NULL

    this_id <- if (!is.null(names(panel_list)) && nzchar(names(panel_list)[i])) {
      names(panel_list)[i]
    } else {
      as.character(i)
    }

    if (verbose) {
      cat("Fitting panel:", this_id, "(", i, "/", length(panel_list), ")\n")
    }

    results[[i]] <- fit_temporal_bipartite_als(
      edge_panel = dat,
      row_cov_df = row_cov_df_i,
      col_cov_df = col_cov_df_i,
      dyad_cov_df = dyad_cov_df_i,
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
      transform = transform,
      row_cov_time_col = row_cov_time_col,
      row_cov_id_col = row_cov_id_col,
      row_cov_prefix = row_cov_prefix,
      col_cov_time_col = col_cov_time_col,
      col_cov_id_col = col_cov_id_col,
      col_cov_prefix = col_cov_prefix,
      dyad_cov_time_col = dyad_cov_time_col,
      dyad_cov_row_id_col = dyad_cov_row_id_col,
      dyad_cov_col_id_col = dyad_cov_col_id_col,
      dyad_cov_row_prefix = dyad_cov_row_prefix,
      dyad_cov_col_prefix = dyad_cov_col_prefix,
      K = K,
      lambda = lambda,
      gamma = gamma,
      max_iter = max_iter,
      tol = tol,
      first_period_n_starts = first_period_n_starts,
      first_period_random_sd = first_period_random_sd,
      verbose = verbose,
      panel_id = this_id
    )
  }

  results
}

