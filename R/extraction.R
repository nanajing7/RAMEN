#' Extract coefficient vector for inference
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @param include_alpha Logical; include row effects.
#' @param include_beta Logical; include column effects.
#'
#' @return A named numeric vector of coefficients.
#' @export
extract_coefficients_als <- function(
    object,
    include_alpha = TRUE,
    include_beta = TRUE
) {

  .validate_als_fit_object(object)

  out <- numeric(0)

  # Row covariates
  if (length(object$coef_row) > 0) {
    vals <- as.numeric(object$coef_row)
    names(vals) <- paste0("row:", names(object$coef_row) %||% seq_along(vals))
    out <- c(out, vals)
  }

  # Column covariates
  if (length(object$coef_col) > 0) {
    vals <- as.numeric(object$coef_col)
    names(vals) <- paste0("col:", names(object$coef_col) %||% seq_along(vals))
    out <- c(out, vals)
  }

  # Dyad covariates
  if (length(object$coef_dyad) > 0) {
    vals <- as.numeric(object$coef_dyad)
    names(vals) <- paste0("dyad:", names(object$coef_dyad) %||% seq_along(vals))
    out <- c(out, vals)
  }

  # Row effects
  if (include_alpha) {
    vals <- as.numeric(object$alpha)
    names(vals) <- paste0("alpha:", names(object$alpha) %||% seq_along(vals))
    out <- c(out, vals)
  }

  # Column effects
  if (include_beta) {
    vals <- as.numeric(object$beta)
    names(vals) <- paste0("beta:", names(object$beta) %||% seq_along(vals))
    out <- c(out, vals)
  }

  out
}



#' Extract latent interaction matrix U V^T
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#'
#' @return An N x M matrix of latent interactions.
#' @export
extract_latent_matrix_als <- function(object) {
  .validate_als_fit_object(object)
  latent_matrix_als_factorize_joint_cov(object)
}


