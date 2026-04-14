
#' Build time-varying row covariate matrices
#'
#' Required columns in `row_cov_df`:
#'   time_col, row_id_col, <covariates...>
#'
#' @param row_cov_df Data frame of row covariates, or NULL.
#' @param years Vector of years.
#' @param row_ids Full row ID universe.
#' @param covar_names Character vector of row covariate names.
#' @param time_col Name of the time column in `row_cov_df`.
#' @param row_id_col Name of the row-node ID column in `row_cov_df`.
#' @param row_prefix Prefix applied to `row_id_col` values before matching to
#'   `row_ids`. Use `""` if IDs are already aligned.
#'
#' @return A named list of matrices, one per year.
build_row_covariates <- function(row_cov_df,
                                 years,
                                 row_ids,
                                 covar_names,
                                 time_col = "year",
                                 row_id_col = "node_row",
                                 row_prefix = "") {
  if (is.null(row_cov_df) || length(covar_names) == 0) {
    out <- vector("list", length(years))
    names(out) <- as.character(years)
    return(out)
  }

  check_required_columns(
    row_cov_df,
    c(time_col, row_id_col, covar_names),
    "row_cov_df"
  )

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  row_cov_df <- row_cov_df |>
    dplyr::mutate(
      .time = .data[[time_col]],
      node_row = paste0(row_prefix, as.character(.data[[row_id_col]]))
    )

  out <- vector("list", length(years))
  names(out) <- as.character(years)

  for (tt in seq_along(years)) {
    yy <- years[tt]

    tmp <- row_cov_df |>
      dplyr::filter(.data$.time == yy) |>
      dplyr::distinct(.data$node_row, .keep_all = TRUE) |>
      dplyr::select(.data$node_row, dplyr::all_of(covar_names))

    tmp <- data.frame(
      node_row = row_ids,
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(tmp, by = "node_row")

    X_row <- as.matrix(tmp[, covar_names, drop = FALSE])
    storage.mode(X_row) <- "double"
    rownames(X_row) <- tmp$node_row

    out[[tt]] <- X_row
  }

  out
}


#' Build time-varying column covariate matrices
#'
#' Required columns in `col_cov_df`:
#'   time_col, col_id_col, <covariates...>
#'
#' @param col_cov_df Data frame of column covariates, or NULL.
#' @param years Vector of years.
#' @param col_ids Full column ID universe.
#' @param covar_names Character vector of column covariate names.
#' @param time_col Name of the time column in `col_cov_df`.
#' @param col_id_col Name of the column-node ID column in `col_cov_df`.
#' @param col_prefix Prefix applied to `col_id_col` values before matching to
#'   `col_ids`. Use `""` if IDs are already aligned.
#'
#' @return A named list of matrices, one per year.
build_col_covariates <- function(col_cov_df,
                                 years,
                                 col_ids,
                                 covar_names,
                                 time_col = "year",
                                 col_id_col = "node_col",
                                 col_prefix = "") {
  if (is.null(col_cov_df) || length(covar_names) == 0) {
    out <- vector("list", length(years))
    names(out) <- as.character(years)
    return(out)
  }

  check_required_columns(
    col_cov_df,
    c(time_col, col_id_col, covar_names),
    "col_cov_df"
  )

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  col_cov_df <- col_cov_df |>
    dplyr::mutate(
      .time = .data[[time_col]],
      node_col = paste0(col_prefix, as.character(.data[[col_id_col]]))
    )

  out <- vector("list", length(years))
  names(out) <- as.character(years)

  for (tt in seq_along(years)) {
    yy <- years[tt]

    tmp <- col_cov_df |>
      dplyr::filter(.data$.time == yy) |>
      dplyr::distinct(.data$node_col, .keep_all = TRUE) |>
      dplyr::select(.data$node_col, dplyr::all_of(covar_names))

    tmp <- data.frame(
      node_col = col_ids,
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(tmp, by = "node_col")

    X_col <- as.matrix(tmp[, covar_names, drop = FALSE])
    storage.mode(X_col) <- "double"
    rownames(X_col) <- tmp$node_col

    out[[tt]] <- X_col
  }

  out
}


#' Build time-varying dyadic covariate arrays
#'
#' Required columns in `dyad_cov_df`:
#'   time_col, row_id_col, col_id_col, <covariates...>
#'
#' @param dyad_cov_df Data frame of dyadic covariates, or NULL.
#' @param years Vector of years.
#' @param row_ids Full row ID universe.
#' @param col_ids Full column ID universe.
#' @param covar_names Character vector of dyadic covariate names.
#' @param time_col Name of the time column.
#' @param row_id_col Name of the row-node ID column.
#' @param col_id_col Name of the column-node ID column.
#' @param row_prefix Prefix applied to `row_id_col` before matching.
#' @param col_prefix Prefix applied to `col_id_col` before matching.
#'
#' @return A named list of arrays, one per year.
build_dyad_covariates <- function(dyad_cov_df,
                                  years,
                                  row_ids,
                                  col_ids,
                                  covar_names,
                                  time_col = "year",
                                  row_id_col = "node_row",
                                  col_id_col = "node_col",
                                  row_prefix = "",
                                  col_prefix = "") {
  if (is.null(dyad_cov_df) || length(covar_names) == 0) {
    return(stats::setNames(vector("list", length(years)), as.character(years)))
  }

  check_required_columns(
    dyad_cov_df,
    c(time_col, row_id_col, col_id_col, covar_names),
    "dyad_cov_df"
  )

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required.")

  dyad_cov_df <- dyad_cov_df |>
    dplyr::mutate(
      .time = .data[[time_col]],
      node_row = paste0(row_prefix, as.character(.data[[row_id_col]])),
      node_col = paste0(col_prefix, as.character(.data[[col_id_col]]))
    )

  out <- vector("list", length(years))
  names(out) <- as.character(years)

  for (tt in seq_along(years)) {
    yy <- years[tt]
    df_t <- dyad_cov_df |>
      dplyr::filter(.data$.time == yy)

    Pd <- length(covar_names)
    arr <- array(
      0,
      dim = c(length(row_ids), length(col_ids), Pd),
      dimnames = list(row_ids, col_ids, covar_names)
    )

    for (p in seq_along(covar_names)) {
      vv <- covar_names[p]

      tmp <- df_t |>
        dplyr::select(.data$node_row, .data$node_col, value = dplyr::all_of(vv)) |>
        dplyr::mutate(
          node_row = factor(.data$node_row, levels = row_ids),
          node_col = factor(.data$node_col, levels = col_ids),
          value = as.numeric(.data$value)
        ) |>
        tidyr::complete(.data$node_row, .data$node_col, fill = list(value = 0)) |>
        tidyr::pivot_wider(names_from = .data$node_col, values_from = .data$value) |>
        tibble::column_to_rownames("node_row") |>
        as.matrix()

      storage.mode(tmp) <- "double"
      arr[, , p] <- tmp
    }

    out[[tt]] <- arr
  }

  out
}

