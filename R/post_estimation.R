#' Validate fitted ALS object
#'
#' Internal helper to check that the fitted object contains the expected fields.
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @return Invisibly returns TRUE if validation passes.
.validate_als_fit_object <- function(object) {
  required_names <- c(
    "U", "V", "alpha", "beta",
    "coef_row", "coef_col", "coef_dyad"
  )

  missing_names <- setdiff(required_names, names(object))
  if (length(missing_names) > 0) {
    stop(
      sprintf(
        "object is missing required component(s): %s",
        paste(missing_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.matrix(object$U) || !is.matrix(object$V)) {
    stop("object$U and object$V must both be matrices.", call. = FALSE)
  }

  if (ncol(object$U) != ncol(object$V)) {
    stop("object$U and object$V must have the same number of columns.", call. = FALSE)
  }

  if (length(object$alpha) != nrow(object$U)) {
    stop("length(object$alpha) must equal nrow(object$U).", call. = FALSE)
  }

  if (length(object$beta) != nrow(object$V)) {
    stop("length(object$beta) must equal nrow(object$V).", call. = FALSE)
  }

  invisible(TRUE)
}


#' Construct the latent interaction matrix U V^T
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @return An N x M matrix of latent interactions.
#' @export
latent_matrix_als_factorize_joint_cov <- function(object) {
  .validate_als_fit_object(object)

  out <- object$U %*% t(object$V)

  rownames(out) <- rownames(object$U)
  colnames(out) <- rownames(object$V)

  out
}


#' Construct the additive component matrix
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @return An N x M matrix containing alpha_i + beta_j.
#' @export
additive_matrix_als_factorize_joint_cov <- function(object) {
  .validate_als_fit_object(object)

  N <- nrow(object$U)
  M <- nrow(object$V)

  out <- matrix(object$alpha, N, M, byrow = FALSE) +
    matrix(object$beta, N, M, byrow = TRUE)

  rownames(out) <- rownames(object$U)
  colnames(out) <- rownames(object$V)

  out
}


#' Construct the covariate component matrix
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @param X_row Optional row covariate matrix.
#' @param X_col Optional column covariate matrix.
#' @param X_dyad Optional dyad covariate array.
#' @return An N x M matrix containing the fitted covariate contribution.
#' @export
covariate_matrix_als_factorize_joint_cov <- function(
    object,
    X_row = NULL,
    X_col = NULL,
    X_dyad = NULL) {

  .validate_als_fit_object(object)

  N <- nrow(object$U)
  M <- nrow(object$V)

  out <- matrix(0, N, M)

  coef_row <- as.numeric(object$coef_row)
  coef_col <- as.numeric(object$coef_col)
  coef_dyad <- as.numeric(object$coef_dyad)

  # Row covariates
  if (!is.null(X_row)) {
    X_row <- as.matrix(X_row)
    if (nrow(X_row) != N) {
      stop("X_row must have nrow(object$U) rows.", call. = FALSE)
    }
    if (ncol(X_row) != length(coef_row)) {
      stop("ncol(X_row) must match length(object$coef_row).", call. = FALSE)
    }

    storage.mode(X_row) <- "double"
    X_row[is.na(X_row)] <- 0

    row_term <- as.vector(X_row %*% coef_row)
    out <- out + matrix(row_term, N, M, byrow = FALSE)
  } else if (length(coef_row) > 0) {
    warning(
      "object contains row coefficients but X_row is NULL; row covariate contribution set to zero.",
      call. = FALSE
    )
  }

  # Column covariates
  if (!is.null(X_col)) {
    X_col <- as.matrix(X_col)
    if (nrow(X_col) != M) {
      stop("X_col must have nrow(object$V) rows.", call. = FALSE)
    }
    if (ncol(X_col) != length(coef_col)) {
      stop("ncol(X_col) must match length(object$coef_col).", call. = FALSE)
    }

    storage.mode(X_col) <- "double"
    X_col[is.na(X_col)] <- 0

    col_term <- as.vector(X_col %*% coef_col)
    out <- out + matrix(col_term, N, M, byrow = TRUE)
  } else if (length(coef_col) > 0) {
    warning(
      "object contains column coefficients but X_col is NULL; column covariate contribution set to zero.",
      call. = FALSE
    )
  }

  # Dyad covariates
  if (!is.null(X_dyad)) {
    if (length(dim(X_dyad)) != 3) {
      stop("X_dyad must be a 3D array.", call. = FALSE)
    }
    if (!all(dim(X_dyad)[1:2] == c(N, M))) {
      stop("The first two dimensions of X_dyad must be N x M.", call. = FALSE)
    }
    if (dim(X_dyad)[3] != length(coef_dyad)) {
      stop("The third dimension of X_dyad must match length(object$coef_dyad).", call. = FALSE)
    }

    storage.mode(X_dyad) <- "double"
    X_dyad[is.na(X_dyad)] <- 0

    for (p in seq_along(coef_dyad)) {
      out <- out + coef_dyad[p] * X_dyad[, , p]
    }
  } else if (length(coef_dyad) > 0) {
    warning(
      "object contains dyad coefficients but X_dyad is NULL; dyad covariate contribution set to zero.",
      call. = FALSE
    )
  }

  rownames(out) <- rownames(object$U)
  colnames(out) <- rownames(object$V)

  out
}


#' Construct the full fitted value matrix
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @param X_row Optional row covariate matrix.
#' @param X_col Optional column covariate matrix.
#' @param X_dyad Optional dyad covariate array.
#' @return An N x M matrix of fitted values.
#' @export
fitted_als_factorize_joint_cov <- function(
    object,
    X_row = NULL,
    X_col = NULL,
    X_dyad = NULL) {

  .validate_als_fit_object(object)

  additive_part <- additive_matrix_als_factorize_joint_cov(object)
  covariate_part <- covariate_matrix_als_factorize_joint_cov(
    object = object,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad)
  latent_part <- latent_matrix_als_factorize_joint_cov(object)

  out <- additive_part + covariate_part + latent_part

  rownames(out) <- rownames(object$U)
  colnames(out) <- rownames(object$V)

  out
}


#' Compute residual matrix from fitted ALS object
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @param Y Outcome matrix used in estimation.
#' @param X_row Optional row covariate matrix.
#' @param X_col Optional column covariate matrix.
#' @param X_dyad Optional dyad covariate array.
#' @return An N x M residual matrix, Y - fitted.
#' @export
residuals_als_factorize_joint_cov <- function(
    object,
    Y,
    X_row = NULL,
    X_col = NULL,
    X_dyad = NULL) {

  .validate_als_fit_object(object)

  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  if (!all(dim(Y) == c(nrow(object$U), nrow(object$V)))) {
    stop("Y must have dimensions nrow(object$U) x nrow(object$V).", call. = FALSE)
  }

  fitted_vals <- fitted_als_factorize_joint_cov(
    object = object,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad
  )

  out <- Y - fitted_vals
  rownames(out) <- rownames(object$U)
  colnames(out) <- rownames(object$V)

  out
}


#' Decompose a fitted ALS object into additive, covariate, latent, and total parts
#'
#' @param object A fitted object returned by als_factorize_joint_cov().
#' @param X_row Optional row covariate matrix.
#' @param X_col Optional column covariate matrix.
#' @param X_dyad Optional dyad covariate array.
#' @return A list with additive, covariate, latent, and fitted matrices.
#' @export
decompose_als_factorize_joint_cov <- function(
    object,
    X_row = NULL,
    X_col = NULL,
    X_dyad = NULL) {

  .validate_als_fit_object(object)

  additive_part <- additive_matrix_als_factorize_joint_cov(object)
  covariate_part <- covariate_matrix_als_factorize_joint_cov(
    object = object,
    X_row = X_row,
    X_col = X_col,
    X_dyad = X_dyad
  )
  latent_part <- latent_matrix_als_factorize_joint_cov(object)

  fitted_part <- additive_part + covariate_part + latent_part

  list(
    additive = additive_part,
    covariate = covariate_part,
    latent = latent_part,
    fitted = fitted_part
  )
}
