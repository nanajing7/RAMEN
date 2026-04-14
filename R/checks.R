
#' Check required columns in a data frame
#'
#' @param data A data frame.
#' @param required Character vector of required column names.
#' @param data_name Name used in error messages.
#'
#' @return Invisibly returns TRUE if all required columns are present.
check_required_columns <- function(data, required, data_name = deparse(substitute(data))) {
  if (!is.data.frame(data)) {
    stop(sprintf("%s must be a data.frame.", data_name), call. = FALSE)
  }

  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        data_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
