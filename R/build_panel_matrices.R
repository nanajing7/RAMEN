
#' Build yearly outcome matrices from a weighted bipartite panel
#'
#' Required columns in `edge_panel` are specified by `time_col`, `row_col`,
#' `col_col`, and `value_col`.
#'
#' @param edge_panel Data frame containing a weighted bipartite panel.
#' @param time_col Name of the time column.
#' @param row_col Name of the row-node column.
#' @param col_col Name of the column-node column.
#' @param value_col Name of the edge-weight/value column.
#' @param row_prefix Prefix added to row node IDs in the output.
#' @param col_prefix Prefix added to column node IDs in the output.
#' @param row_pad Optional integer width for left-padding row node IDs. Use
#'   `NULL` for no padding.
#' @param col_pad Optional integer width for left-padding column node IDs. Use
#'   `NULL` for no padding.
#' @param transform Function applied elementwise to the completed matrix.
#'   Default is `log1p`.
#'
#' @return A list with years, Y_list, row_ids, and col_ids.
#' @export
build_panel_matrices <- function(edge_panel,
                                 time_col = "year",
                                 row_col = "node_row",
                                 col_col = "node_col",
                                 value_col = "value",
                                 row_prefix = "row_",
                                 col_prefix = "col_",
                                 row_pad = NULL,
                                 col_pad = NULL,
                                 transform = log1p) {
  check_required_columns(
    edge_panel,
    c(time_col, row_col, col_col, value_col),
    "edge_panel"
  )

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")

  df0 <- edge_panel |>
    dplyr::mutate(
      .time_raw = .data[[time_col]],
      .row_raw  = as.character(.data[[row_col]]),
      .col_raw  = as.character(.data[[col_col]]),
      .value    = as.numeric(.data[[value_col]])
    )

  if (!is.null(row_pad)) {
    df0 <- df0 |>
      dplyr::mutate(.row_raw = stringr::str_pad(.data$.row_raw, width = row_pad, side = "left", pad = "0"))
  }

  if (!is.null(col_pad)) {
    df0 <- df0 |>
      dplyr::mutate(.col_raw = stringr::str_pad(.data$.col_raw, width = col_pad, side = "left", pad = "0"))
  }

  df0 <- df0 |>
    dplyr::mutate(
      node_row = paste0(row_prefix, .data$.row_raw),
      node_col = paste0(col_prefix, .data$.col_raw)
    ) |>
    dplyr::filter(
      ## is.finite(.data$.value), .data$.value >= 0,
      is.finite(.data$.value),
      !is.na(.data$.time_raw), !is.na(.data$node_row), !is.na(.data$node_col)
    ) |>
    dplyr::distinct(.data$.time_raw, .data$node_row, .data$node_col, .data$.value)

  if (nrow(df0) == 0) {
    stop("No valid observations remain after preprocessing.", call. = FALSE)
  }

  years <- sort(unique(df0$.time_raw))
  row_nodes_all <- sort(unique(df0$node_row))
  col_nodes_all <- sort(unique(df0$node_col))

  Y_list <- lapply(years, function(tt) {
    df_t <- df0 |>
      dplyr::filter(.data$.time_raw == tt) |>
      dplyr::group_by(.data$node_row, .data$node_col) |>
      dplyr::summarise(value = sum(.data$.value), .groups = "drop")

    Ymat <- df_t |>
      dplyr::mutate(
        node_row = factor(.data$node_row, levels = row_nodes_all),
        node_col = factor(.data$node_col, levels = col_nodes_all)
      ) |>
      tidyr::complete(.data$node_row, .data$node_col, fill = list(value = 0)) |>
      tidyr::pivot_wider(names_from = .data$node_col, values_from = .data$value) |>
      tibble::column_to_rownames("node_row") |>
      as.matrix()

    Y <- transform(Ymat)
    storage.mode(Y) <- "double"
    Y
  })

  list(
    years = years,
    Y_list = Y_list,
    row_ids = row_nodes_all,
    col_ids = col_nodes_all
  )
}
