#' Null coalescing operator
#'
#' Returns the left-hand side if not NULL, otherwise the right-hand side.
#'
#' @param a First value
#' @param b Fallback value if a is NULL
#'
#' @return Either a or b
#' @export
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}
