# Steps 5 and 6 of the inner block coordinate descent: whole-trajectory updates
# of the sender latent positions `U` and receiver latent positions `V`.
#
# Holding the opposite factor fixed, the objective is quadratic in the block
# being updated and separates across nodes. For sender i the first-order
# conditions give a block-tridiagonal system in u_i = (u_i1', ..., u_iT')':
#
#   t = 1      (G_i1 + (lambda_U + gamma_U) I_K) u_i1 - gamma_U u_i2       = c_i1
#   1 < t < T  -gamma_U u_i,t-1 + (G_it + 2 gamma_U I_K) u_it
#                                            - gamma_U u_i,t+1             = c_it
#   t = T      -gamma_U u_i,T-1 + (G_iT + gamma_U I_K) u_iT                = c_iT
#
# with G_it = sum_j v_jt v_jt' and c_it = sum_j v_jt (Y_ijt - a_it - b_jt -
# x_ijt' beta_t), both taken over the cells observed for that node and period.
# The spec solves this system once per node, which is what the loops below do.
#
# The two factors are updated alternately because the multiplicative term is
# jointly non-convex in (U, V) but conditionally quadratic in either one.

#' Per-node Gram matrices and right-hand sides for a latent block
#'
#' Assembles, for every node on one side of the panel, the `K x K` Gram matrix
#' of the opposite factor over that node's observed cells and the matching
#' `K`-vector of cross-products with the partial residual.
#'
#' The Gram entries are formed with `K(K+1)/2` matrix products rather than a
#' per-node loop: for the sender side, `obs %*% (V[, k1] * V[, k2])` yields
#' entry `(k1, k2)` for all senders at once. This is the same trick the legacy
#' single-period path uses.
#'
#' @param R `N x M` partial residual (`NA` where unobserved) with the latent
#'   term excluded; for the receiver side pass its transpose.
#' @param W The opposite factor for this period: `V_t` (`M x K`) when updating
#'   `U`, `U_t` (`N x K`) when updating `V`.
#'
#' @return A list with `G` (an `n x K x K` array) and `c` (an `n x K` matrix),
#'   where `n` is `nrow(R)`.
#' @keywords internal
#' @noRd
.latent_gram_rhs <- function(R, W) {
  K <- ncol(W)
  n <- nrow(R)

  obs <- !is.na(R)
  R0 <- R
  R0[is.na(R0)] <- 0

  # c_i = sum_{j observed} W_j R_ij
  rhs <- R0 %*% W

  # G_i[k1, k2] = sum_{j observed} W_j[k1] W_j[k2]
  G <- array(0, dim = c(n, K, K))
  for (k1 in seq_len(K)) {
    for (k2 in seq_len(k1)) {
      vals <- as.vector(obs %*% (W[, k1] * W[, k2]))
      G[, k1, k2] <- vals
      G[, k2, k1] <- vals
    }
  }

  list(G = G, c = rhs)
}


#' Solve one node's latent trajectory
#'
#' Builds the node's block-tridiagonal system from its per-period Gram matrices
#' and right-hand sides, then calls the shared block solver.
#'
#' @param G_i List of length `T` of `K x K` Gram matrices for this node.
#' @param c_i List of length `T` of `K`-vectors.
#' @param lambda Initial-state penalty for this block.
#' @param gamma Random-walk penalty for this block.
#'
#' @return A `T x K` matrix; row `t` is the node's position in period `t`.
#' @keywords internal
#' @noRd
.solve_latent_node <- function(G_i, c_i, lambda, gamma) {
  Tn <- length(G_i)
  K <- nrow(G_i[[1]])
  Ik <- diag(K)

  D <- vector("list", Tn)
  for (t in seq_len(Tn)) {
    # Temporal neighbours contribute one gamma each; t = 1 also carries lambda.
    nb <- (t > 1L) + (t < Tn)
    D[[t]] <- G_i[[t]] + (gamma * nb + if (t == 1L) lambda else 0) * Ik
  }

  x <- .solve_block_tridiagonal(D, off = -gamma, r = c_i)
  do.call(rbind, x)
}


#' Update the sender latent-position trajectory `U` (Step 5)
#'
#' @param params Current parameter list (`a`, `b`, `beta`, `U`, `V`).
#' @param Y_list List of `T` outcome matrices.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param lambda_U Initial-state penalty for `U`.
#' @param gamma_U Random-walk penalty for `U`.
#'
#' @return `params` with the `U` component replaced.
#' @keywords internal
#' @noRd
.update_U_block <- function(params,
                            Y_list,
                            X_row_list = NULL,
                            X_col_list = NULL,
                            X_dyad_list = NULL,
                            lambda_U,
                            gamma_U) {
  Tn <- length(Y_list)
  N <- nrow(params$a)
  M <- nrow(params$b)
  K <- ncol(params$U[[1]])

  G_all <- vector("list", Tn)
  c_all <- vector("list", Tn)

  for (t in seq_len(Tn)) {
    # Partial residual: everything except the latent term.
    R <- Y_list[[t]] -
      matrix(params$a[, t], N, M, byrow = FALSE) -
      matrix(params$b[, t], N, M, byrow = TRUE) -
      .cov_term(
        beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
        X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
        X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
        X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]],
        N = N, M = M
      )

    gr <- .latent_gram_rhs(R, params$V[[t]])
    G_all[[t]] <- gr$G
    c_all[[t]] <- gr$c
  }

  U_new <- vector("list", Tn)
  for (t in seq_len(Tn)) U_new[[t]] <- matrix(0, N, K)

  for (i in seq_len(N)) {
    G_i <- lapply(seq_len(Tn), function(t) matrix(G_all[[t]][i, , ], K, K))
    c_i <- lapply(seq_len(Tn), function(t) c_all[[t]][i, ])
    Ui <- .solve_latent_node(G_i, c_i, lambda_U, gamma_U)
    for (t in seq_len(Tn)) U_new[[t]][i, ] <- Ui[t, ]
  }

  for (t in seq_len(Tn)) dimnames(U_new[[t]]) <- dimnames(params$U[[t]])
  params$U <- U_new
  params
}


#' Update the receiver latent-position trajectory `V` (Step 6)
#'
#' The receiver analogue of `.update_U_block()`: identical structure with the
#' roles of senders and receivers exchanged, so the partial residual is
#' transposed and the Gram matrices are built from `U_t`.
#'
#' @inheritParams .update_U_block
#' @param lambda_V Initial-state penalty for `V`.
#' @param gamma_V Random-walk penalty for `V`.
#'
#' @return `params` with the `V` component replaced.
#' @keywords internal
#' @noRd
.update_V_block <- function(params,
                            Y_list,
                            X_row_list = NULL,
                            X_col_list = NULL,
                            X_dyad_list = NULL,
                            lambda_V,
                            gamma_V) {
  Tn <- length(Y_list)
  N <- nrow(params$a)
  M <- nrow(params$b)
  K <- ncol(params$V[[1]])

  G_all <- vector("list", Tn)
  c_all <- vector("list", Tn)

  for (t in seq_len(Tn)) {
    R <- Y_list[[t]] -
      matrix(params$a[, t], N, M, byrow = FALSE) -
      matrix(params$b[, t], N, M, byrow = TRUE) -
      .cov_term(
        beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
        X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
        X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
        X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]],
        N = N, M = M
      )

    # Transposing puts receivers on the rows, so the sender-side code applies.
    gr <- .latent_gram_rhs(t(R), params$U[[t]])
    G_all[[t]] <- gr$G
    c_all[[t]] <- gr$c
  }

  V_new <- vector("list", Tn)
  for (t in seq_len(Tn)) V_new[[t]] <- matrix(0, M, K)

  for (j in seq_len(M)) {
    G_j <- lapply(seq_len(Tn), function(t) matrix(G_all[[t]][j, , ], K, K))
    c_j <- lapply(seq_len(Tn), function(t) c_all[[t]][j, ])
    Vj <- .solve_latent_node(G_j, c_j, lambda_V, gamma_V)
    for (t in seq_len(Tn)) V_new[[t]][j, ] <- Vj[t, ]
  }

  for (t in seq_len(Tn)) dimnames(V_new[[t]]) <- dimnames(params$V[[t]])
  params$V <- V_new
  params
}
