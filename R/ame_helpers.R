# Shared internal helpers for the joint-trajectory MAP estimator
# (fit_dynamic_ame). These are used by the objective, the block updates, and
# the post-estimation accessors, so they live in one place.

#' Covariate contribution matrix for one period
#'
#' Builds the `N x M` matrix of `x_ijt' beta_t` contributions from the stacked
#' coefficient vector. The stacking order is row, then column, then dyadic
#' covariates, matching how `P = Pr + Pc + Pd` is assembled throughout.
#'
#' Row covariates vary over `i` only, so they are constant along each row of the
#' output; column covariates vary over `j` only; dyadic covariates vary freely.
#'
#' **Missing covariates propagate as `NA`, by design.** If node `i` has no
#' recorded row covariate in period `t` then `x_ijt' beta_t` is undefined for
#' every `j`, so the whole row of the output is `NA` and every one of those
#' cells is excluded from the likelihood downstream — the residual for a cell
#' with an unknown covariate contribution cannot be formed. A missing dyadic
#' covariate knocks out only its own cell. This is the same treatment `Y`'s own
#' missing cells receive, and it makes no imputation assumption.
#'
#' Note this differs from the legacy `als_factorize_joint_cov()`, which silently
#' replaced missing covariates with zero.
#'
#' @param beta_t Stacked coefficient vector of length `P` (may be length 0).
#' @param X_row `N x Pr` matrix, or `NULL`.
#' @param X_col `M x Pc` matrix, or `NULL`.
#' @param X_dyad `N x M x Pd` array, or `NULL`.
#' @param N,M Output dimensions.
#'
#' @return An `N x M` numeric matrix, `NA` wherever a contributing covariate is
#'   missing. All zeros when `P = 0`.
#' @keywords internal
#' @noRd
.cov_term <- function(beta_t, X_row = NULL, X_col = NULL, X_dyad = NULL, N, M) {
  out <- matrix(0, N, M)
  if (length(beta_t) == 0L) return(out)

  idx <- 0L

  if (!is.null(X_row) && ncol(X_row) > 0L) {
    Pr <- ncol(X_row)
    b <- beta_t[(idx + 1L):(idx + Pr)]
    out <- out + matrix(as.vector(X_row %*% b), N, M, byrow = FALSE)
    idx <- idx + Pr
  }

  if (!is.null(X_col) && ncol(X_col) > 0L) {
    Pc <- ncol(X_col)
    b <- beta_t[(idx + 1L):(idx + Pc)]
    out <- out + matrix(as.vector(X_col %*% b), N, M, byrow = TRUE)
    idx <- idx + Pc
  }

  if (!is.null(X_dyad) && dim(X_dyad)[3] > 0L) {
    Pd <- dim(X_dyad)[3]
    b <- beta_t[(idx + 1L):(idx + Pd)]
    for (p in seq_len(Pd)) {
      out <- out + b[p] * X_dyad[, , p]
    }
  }

  out
}


#' Stacked design matrix of covariates for one period
#'
#' Returns the `(N*M) x P` design matrix whose columns line up with the stacked
#' coefficient vector used by `.cov_term()`, in the same row-major-free ordering
#' as `as.vector()` of an `N x M` matrix (column-major: cell `(i, j)` is at
#' position `i + (j - 1) * N`).
#'
#' Row covariates are recycled down the rows, column covariates across the
#' columns, and dyadic covariates are flattened directly.
#'
#' **Missing covariates are kept as `NA`**, matching `.cov_term()`. The caller
#' must drop any row of the design that contains an `NA` — those cells have an
#' undefined covariate contribution and are excluded from the likelihood. Since
#' the residual for such a cell is itself `NA`, masking on the residual removes
#' exactly the same rows.
#'
#' @param X_row `N x Pr` matrix, or `NULL`.
#' @param X_col `M x Pc` matrix, or `NULL`.
#' @param X_dyad `N x M x Pd` array, or `NULL`.
#' @param N,M Panel dimensions.
#'
#' @return An `(N*M) x P` numeric matrix, `NA` where a covariate is missing;
#'   `(N*M) x 0` when there are no covariates.
#' @keywords internal
#' @noRd
.cov_design <- function(X_row = NULL, X_col = NULL, X_dyad = NULL, N, M) {
  blocks <- list()

  if (!is.null(X_row) && ncol(X_row) > 0L) {
    # cell (i, j) at position i + (j-1)*N  ->  row index i repeats M times
    blocks <- c(blocks, list(X_row[rep(seq_len(N), times = M), , drop = FALSE]))
  }

  if (!is.null(X_col) && ncol(X_col) > 0L) {
    # cell (i, j) -> column index j repeats N times consecutively
    blocks <- c(blocks, list(X_col[rep(seq_len(M), each = N), , drop = FALSE]))
  }

  if (!is.null(X_dyad) && dim(X_dyad)[3] > 0L) {
    Pd <- dim(X_dyad)[3]
    Zd <- matrix(0, N * M, Pd)
    for (p in seq_len(Pd)) Zd[, p] <- as.vector(X_dyad[, , p])
    blocks <- c(blocks, list(Zd))
  }

  if (length(blocks) == 0L) return(matrix(0, N * M, 0))

  Z <- do.call(cbind, blocks)
  storage.mode(Z) <- "double"
  Z
}


#' Number of covariates implied by the design pieces
#'
#' @param X_row,X_col,X_dyad Covariate pieces, any of which may be `NULL`.
#'
#' @return A list with `Pr`, `Pc`, `Pd`, and their sum `P`.
#' @keywords internal
#' @noRd
.cov_dims <- function(X_row = NULL, X_col = NULL, X_dyad = NULL) {
  Pr <- if (!is.null(X_row))  ncol(X_row)      else 0L
  Pc <- if (!is.null(X_col))  ncol(X_col)      else 0L
  Pd <- if (!is.null(X_dyad)) dim(X_dyad)[3]   else 0L
  list(Pr = as.integer(Pr), Pc = as.integer(Pc), Pd = as.integer(Pd),
       P = as.integer(Pr + Pc + Pd))
}


#' Fitted-value matrix for one period
#'
#' Assembles `a_it + b_jt + x_ijt' beta_t + U_it' V_jt` for a single period.
#'
#' @param a_t Length-`N` vector of row effects.
#' @param b_t Length-`M` vector of column effects.
#' @param U_t `N x K` latent factor matrix.
#' @param V_t `M x K` latent factor matrix.
#' @param beta_t Stacked coefficient vector of length `P`.
#' @param X_row,X_col,X_dyad Covariate pieces for this period.
#'
#' @return An `N x M` matrix of fitted values.
#' @keywords internal
#' @noRd
.fitted_period <- function(a_t, b_t, U_t, V_t, beta_t = numeric(0),
                           X_row = NULL, X_col = NULL, X_dyad = NULL) {
  N <- length(a_t)
  M <- length(b_t)

  matrix(a_t, N, M, byrow = FALSE) +
    matrix(b_t, N, M, byrow = TRUE) +
    .cov_term(beta_t, X_row, X_col, X_dyad, N, M) +
    U_t %*% t(V_t)
}


#' Residual matrix for one period (missing cells stay `NA`)
#'
#' @inheritParams .fitted_period
#' @param Y_t `N x M` outcome matrix with `NA` in unobserved cells.
#'
#' @return An `N x M` residual matrix, `NA` wherever `Y_t` is `NA`.
#' @keywords internal
#' @noRd
.resid_period <- function(Y_t, a_t, b_t, U_t, V_t, beta_t = numeric(0),
                          X_row = NULL, X_col = NULL, X_dyad = NULL) {
  Y_t - .fitted_period(a_t, b_t, U_t, V_t, beta_t, X_row, X_col, X_dyad)
}
