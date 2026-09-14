# Steps 3 and 4 of the inner block coordinate descent: whole-trajectory updates
# of the additive row effects `a` and column effects `b`.
#
# Holding every other block fixed, the objective separates across nodes: row i
# contributes only through its own trajectory a_{i,1:T}. Differentiating and
# setting the gradient to zero gives, for each node, a scalar tridiagonal
# system over the T periods:
#
#   t = 1      (n_1 + lambda + gamma) a_1 - gamma a_2               = S_1
#   1 < t < T  (n_t + 2 gamma) a_t - gamma a_{t-1} - gamma a_{t+1}   = S_t
#   t = T      (n_T + gamma) a_T - gamma a_{T-1}                     = S_T
#
# where n_t is the number of observed cells for that node in period t and S_t
# the sum of its partial residuals. With T = 1 the coupling vanishes and the
# update reduces to a ridge-penalised mean.

#' Assemble the diagonal of the additive-block tridiagonal system
#'
#' @param n_obs `n x T` matrix of observed-cell counts per node and period.
#' @param lambda Initial-state penalty (applies at `t = 1` only).
#' @param gamma Random-walk penalty (one unit per temporal neighbour).
#'
#' @return An `n x T` matrix of diagonal entries.
#' @keywords internal
#' @noRd
.additive_diagonal <- function(n_obs, lambda, gamma) {
  Tn <- ncol(n_obs)

  # Number of temporal neighbours: interior periods have two, the endpoints
  # one, and a single-period panel none.
  neighbours <- if (Tn == 1L) 0 else c(1, rep(2, max(Tn - 2L, 0L)), 1)

  D <- n_obs + gamma * matrix(neighbours, nrow(n_obs), Tn, byrow = TRUE)
  D[, 1] <- D[, 1] + lambda
  D
}


#' Update the row-effect trajectory `a` (Step 3)
#'
#' Solves the whole `a` trajectory exactly, given the current values of every
#' other block. Missing cells are excluded from both the counts and the
#' residual sums.
#'
#' @param params Current parameter list (`a`, `b`, `beta`, `U`, `V`).
#' @param Y_list List of `T` outcome matrices.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param lambda_a Initial-state penalty for `a`.
#' @param gamma_a Random-walk penalty for `a`.
#'
#' @return `params` with the `a` component replaced.
#' @keywords internal
#' @noRd
.update_a_block <- function(params,
                            Y_list,
                            X_row_list = NULL,
                            X_col_list = NULL,
                            X_dyad_list = NULL,
                            lambda_a,
                            gamma_a) {
  Tn <- length(Y_list)
  N <- nrow(params$a)
  M <- nrow(params$b)

  n_obs <- matrix(0, N, Tn)
  S <- matrix(0, N, Tn)

  for (t in seq_len(Tn)) {
    # Partial residual: everything except a_it.
    R <- Y_list[[t]] -
      matrix(params$b[, t], N, M, byrow = TRUE) -
      .cov_term(
        beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
        X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
        X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
        X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]],
        N = N, M = M
      ) -
      params$U[[t]] %*% t(params$V[[t]])

    n_obs[, t] <- rowSums(!is.na(R))
    S[, t] <- rowSums(R, na.rm = TRUE)
  }

  D <- .additive_diagonal(n_obs, lambda_a, gamma_a)
  a_new <- .solve_scalar_tridiagonal_vec(D, off = -gamma_a, r = S)

  dimnames(a_new) <- dimnames(params$a)
  params$a <- a_new
  params
}


#' Update the column-effect trajectory `b` (Step 4)
#'
#' The column analogue of `.update_a_block()`: identical structure with the
#' roles of rows and columns exchanged, so counts and residual sums are taken
#' down columns instead of across rows.
#'
#' @inheritParams .update_a_block
#' @param lambda_b Initial-state penalty for `b`.
#' @param gamma_b Random-walk penalty for `b`.
#'
#' @return `params` with the `b` component replaced.
#' @keywords internal
#' @noRd
.update_b_block <- function(params,
                            Y_list,
                            X_row_list = NULL,
                            X_col_list = NULL,
                            X_dyad_list = NULL,
                            lambda_b,
                            gamma_b) {
  Tn <- length(Y_list)
  N <- nrow(params$a)
  M <- nrow(params$b)

  n_obs <- matrix(0, M, Tn)
  S <- matrix(0, M, Tn)

  for (t in seq_len(Tn)) {
    # Partial residual: everything except b_jt.
    R <- Y_list[[t]] -
      matrix(params$a[, t], N, M, byrow = FALSE) -
      .cov_term(
        beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
        X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
        X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
        X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]],
        N = N, M = M
      ) -
      params$U[[t]] %*% t(params$V[[t]])

    n_obs[, t] <- colSums(!is.na(R))
    S[, t] <- colSums(R, na.rm = TRUE)
  }

  D <- .additive_diagonal(n_obs, lambda_b, gamma_b)
  b_new <- .solve_scalar_tridiagonal_vec(D, off = -gamma_b, r = S)

  dimnames(b_new) <- dimnames(params$b)
  params$b <- b_new
  params
}
