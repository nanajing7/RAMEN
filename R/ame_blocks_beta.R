# Step 2 of the inner block coordinate descent: whole-trajectory update of the
# regression coefficients.
#
# Unlike the other four blocks, beta_t is shared by every cell of period t, so
# there is a single system of size T*P rather than one per node. Holding the
# remaining parameters fixed, the first-order conditions are
#
#   t = 1      (S_1 + (lambda_b + gamma_b) I_P) beta_1 - gamma_b beta_2   = r_1
#   1 < t < T  -gamma_b beta_{t-1} + (S_t + 2 gamma_b I_P) beta_t
#                                              - gamma_b beta_{t+1}       = r_t
#   t = T      -gamma_b beta_{T-1} + (S_T + gamma_b I_P) beta_T           = r_T
#
# with S_t = sum_ij x_ijt x_ijt' and r_t = sum_ij x_ijt (Y_ijt - a_it - b_jt -
# u_it' v_jt), both over the cells observed in period t. Note the partial
# residual here subtracts the latent term (it is known while beta is updated),
# unlike the latent blocks where it does not.

#' Update the regression-coefficient trajectory `beta` (Step 2)
#'
#' Solves the whole `beta` trajectory exactly, given the current values of every
#' other block. Returns `params` unchanged when there are no covariates.
#'
#' The per-period information matrix and right-hand side are formed from the
#' stacked design matrix produced by `.cov_design()`, restricted to the cells
#' that are observed *and* have a fully defined covariate vector. A cell whose
#' covariates are partly missing has an undefined contribution and is dropped —
#' the residual for such a cell is `NA` for the same reason, so masking on the
#' residual removes exactly those rows.
#'
#' @param params Current parameter list (`a`, `b`, `beta`, `U`, `V`).
#' @param Y_list List of `T` outcome matrices.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param lambda_beta Initial-state penalty for `beta`.
#' @param gamma_beta Random-walk penalty for `beta`.
#'
#' @return `params` with the `beta` component replaced.
#' @keywords internal
#' @noRd
.update_beta_block <- function(params,
                               Y_list,
                               X_row_list = NULL,
                               X_col_list = NULL,
                               X_dyad_list = NULL,
                               lambda_beta,
                               gamma_beta) {
  P <- nrow(params$beta)
  if (P == 0L) return(params)          # no covariates: nothing to update

  Tn <- length(Y_list)
  N <- nrow(params$a)
  M <- nrow(params$b)
  Ip <- diag(P)

  S <- vector("list", Tn)              # information matrices
  r <- vector("list", Tn)              # right-hand sides

  for (t in seq_len(Tn)) {
    X_row  <- if (is.null(X_row_list))  NULL else X_row_list[[t]]
    X_col  <- if (is.null(X_col_list))  NULL else X_col_list[[t]]
    X_dyad <- if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]]

    # Partial residual: everything except the covariate term.
    R <- Y_list[[t]] -
      matrix(params$a[, t], N, M, byrow = FALSE) -
      matrix(params$b[, t], N, M, byrow = TRUE) -
      params$U[[t]] %*% t(params$V[[t]])

    Z <- .cov_design(X_row, X_col, X_dyad, N, M)
    y <- as.vector(R)

    # Keep only cells with an observed outcome and a complete covariate vector.
    keep <- !is.na(y) & stats::complete.cases(Z)
    Zk <- Z[keep, , drop = FALSE]

    S[[t]] <- crossprod(Zk)
    r[[t]] <- as.vector(crossprod(Zk, y[keep]))
  }

  D <- vector("list", Tn)
  for (t in seq_len(Tn)) {
    nb <- (t > 1L) + (t < Tn)
    D[[t]] <- S[[t]] + (gamma_beta * nb + if (t == 1L) lambda_beta else 0) * Ip
  }

  sol <- .solve_block_tridiagonal(D, off = -gamma_beta, r = r)

  beta_new <- matrix(unlist(sol), P, Tn)
  dimnames(beta_new) <- dimnames(params$beta)
  params$beta <- beta_new
  params
}
