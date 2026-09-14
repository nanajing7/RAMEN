# Numerical primitives for the joint-trajectory MAP estimator (fit_dynamic_ame).
#
# These are internal helpers. Each block update in the inner block coordinate
# descent solves, for a single parameter trajectory {theta_t}_{t=1..T}, a
# symmetric block-tridiagonal linear system arising from the Gaussian
# random-walk prior (two-sided temporal coupling). One routine covers all five
# blocks (beta, a_i, b_j, U_i, V_j); only the block size q differs
# (q = P for beta, q = 1 for a/b, q = K for U/V).

#' Solve a symmetric block-tridiagonal linear system (block Thomas)
#'
#' Solves the system whose diagonal blocks are `D[[t]]` (each `q x q`), whose
#' sub- and super-diagonal blocks are both `off * I_q`, and whose right-hand
#' sides are `r[[t]]` (each of length `q`). This is exactly the structure of
#' every trajectory block update in the joint MAP objective, where
#' `off = -gamma` is the random-walk coupling.
#'
#' The algorithm is the block generalisation of the Thomas algorithm: forward
#' elimination followed by back substitution, costing `O(T q^3)` without ever
#' forming the full `Tq x Tq` matrix.
#'
#' @param D List of length `T` of `q x q` diagonal blocks. A length-one scalar
#'   or `1 x 1` matrix is accepted when `q = 1` (the `a` and `b` blocks).
#' @param off Scalar multiplying the identity on both off-diagonals (`-gamma`).
#' @param r List of length `T` of right-hand-side vectors, each of length `q`.
#'
#' @return A list of length `T` of solution vectors, each of length `q`.
#' @keywords internal
#' @noRd
.solve_block_tridiagonal <- function(D, off, r) {
  Tn <- length(D)

  if (Tn == 0L) {
    stop("D must have at least one block.", call. = FALSE)
  }
  if (length(r) != Tn) {
    stop("D and r must have the same length.", call. = FALSE)
  }
  if (length(off) != 1L || !is.finite(off)) {
    stop("off must be a single finite scalar.", call. = FALSE)
  }

  Dm <- lapply(D, function(M) if (is.matrix(M)) M else matrix(as.numeric(M), 1, 1))
  rv <- lapply(r, as.numeric)

  q <- nrow(Dm[[1]])
  for (t in seq_len(Tn)) {
    if (nrow(Dm[[t]]) != q || ncol(Dm[[t]]) != q) {
      stop("All blocks in D must be square and of the same size.", call. = FALSE)
    }
    if (length(rv[[t]]) != q) {
      stop("Each element of r must have length equal to the block size.", call. = FALSE)
    }
  }

  if (Tn == 1L) {
    return(list(as.numeric(solve(Dm[[1]], rv[[1]]))))
  }

  off2 <- off * off

  # ── Forward elimination ──────────────────────────────────────────────────
  # With sub-diagonal L_t = off * I and super-diagonal U_t = off * I:
  #   D'_1 = D_1,                r'_1 = r_1
  #   D'_t = D_t - off^2 D'_{t-1}^{-1}
  #   r'_t = r_t - off * D'_{t-1}^{-1} r'_{t-1}
  Dp <- vector("list", Tn)
  rp <- vector("list", Tn)
  Dp[[1]] <- Dm[[1]]
  rp[[1]] <- rv[[1]]

  for (t in 2:Tn) {
    Dinv_prev <- solve(Dp[[t - 1]])
    Dp[[t]] <- Dm[[t]] - off2 * Dinv_prev
    rp[[t]] <- rv[[t]] - off * as.numeric(Dinv_prev %*% rp[[t - 1]])
  }

  # ── Back substitution ────────────────────────────────────────────────────
  #   x_T = D'_T^{-1} r'_T
  #   x_t = D'_t^{-1} (r'_t - off * x_{t+1})
  x <- vector("list", Tn)
  x[[Tn]] <- as.numeric(solve(Dp[[Tn]], rp[[Tn]]))

  for (t in (Tn - 1):1) {
    x[[t]] <- as.numeric(solve(Dp[[t]], rp[[t]] - off * x[[t + 1]]))
  }

  x
}


#' Forward-elimination pivots of a symmetric block-tridiagonal matrix
#'
#' Runs only the forward sweep of the block Thomas algorithm, returning the
#' modified diagonal blocks
#' \deqn{D'_1 = D_1, \qquad D'_t = D_t - \mathrm{off}^2 D_{t-1}'^{-1}.}
#' These are the pivots of the block LDL' factorisation, and both the solve and
#' the selected inversion are built from them.
#'
#' @param D List of `T` diagonal blocks (`q x q`).
#' @param off Scalar on both off-diagonals.
#'
#' @return A list of `T` pivot blocks.
#' @keywords internal
#' @noRd
.block_tridiagonal_pivots <- function(D, off) {
  Tn <- length(D)
  Dm <- lapply(D, function(M) if (is.matrix(M)) M else matrix(as.numeric(M), 1, 1))

  Dp <- vector("list", Tn)
  Dp[[1]] <- Dm[[1]]
  if (Tn > 1L) {
    off2 <- off * off
    for (t in 2:Tn) Dp[[t]] <- Dm[[t]] - off2 * solve(Dp[[t - 1]])
  }
  Dp
}


#' Diagonal and first off-diagonal blocks of the inverse (selected inversion)
#'
#' For a symmetric block-tridiagonal `H` with diagonal blocks `D[[t]]` and both
#' off-diagonals equal to `off * I`, returns the blocks of `H^{-1}` that the
#' empirical-Bayes step needs — the diagonal blocks and the first sub-diagonal —
#' without forming the full inverse. This is the Takahashi/BTA backward
#' recursion, run on the pivots from the forward sweep:
#' \deqn{G_T = D_T'^{-1}, \qquad
#'       G_t = D_t'^{-1}\left(I + \mathrm{off}^2 G_{t+1} D_t'^{-1}\right),}
#' \deqn{(H^{-1})_{t+1,t} = -G_{t+1}\,\mathrm{off}\,D_t'^{-1}.}
#'
#' Cost is `O(T q^3)`, the same order as the solve itself.
#'
#' @param D List of `T` diagonal blocks (`q x q`), or scalars when `q = 1`.
#' @param off Scalar on both off-diagonals (`-gamma`).
#'
#' @return A list with `diag` (list of `T` blocks `(H^{-1})_{tt}`) and `offdiag`
#'   (list of `T - 1` blocks `(H^{-1})_{t+1,t}`; empty when `T = 1`).
#' @keywords internal
#' @noRd
.block_tridiagonal_selected_inverse <- function(D, off) {
  Tn <- length(D)
  Dp <- .block_tridiagonal_pivots(D, off)

  G <- vector("list", Tn)
  G[[Tn]] <- solve(Dp[[Tn]])

  if (Tn == 1L) {
    return(list(diag = G, offdiag = list()))
  }

  off2 <- off * off
  Goff <- vector("list", Tn - 1L)
  Iq <- diag(nrow(G[[Tn]]))

  for (t in (Tn - 1):1) {
    Dinv <- solve(Dp[[t]])
    G[[t]] <- Dinv %*% (Iq + off2 * G[[t + 1]] %*% Dinv)
    # (H^{-1})_{t+1,t}
    Goff[[t]] <- -G[[t + 1]] %*% (off * Dinv)
  }

  list(diag = G, offdiag = Goff)
}


#' Solve many scalar tridiagonal systems at once
#'
#' Specialisation of `.solve_block_tridiagonal()` to block size one, vectorised
#' across independent systems. The additive blocks `a` and `b` give one such
#' system per node, all sharing the same off-diagonal constant, so the scalar
#' Thomas recursion can be run on whole columns at a time: `O(T)` vector
#' operations instead of `n * T` one-by-one matrix solves.
#'
#' Row `i` of the inputs describes system `i`; column `t` is period `t`.
#'
#' @param D `n x T` matrix of diagonal entries.
#' @param off Scalar on both off-diagonals (`-gamma`).
#' @param r `n x T` matrix of right-hand sides.
#'
#' @return An `n x T` matrix of solutions.
#' @keywords internal
#' @noRd
.solve_scalar_tridiagonal_vec <- function(D, off, r) {
  D <- as.matrix(D)
  r <- as.matrix(r)

  if (!all(dim(D) == dim(r))) {
    stop("D and r must have the same dimensions.", call. = FALSE)
  }
  if (length(off) != 1L || !is.finite(off)) {
    stop("off must be a single finite scalar.", call. = FALSE)
  }

  Tn <- ncol(D)

  # A zero pivot means the system is singular: a node with no observations in
  # some period and no penalty tying it to anything.
  if (any(abs(D) < .Machine$double.eps)) {
    stop(
      "Singular additive block: a node has no observations and no penalty ",
      "anchoring it. Ensure the initial or random-walk penalty is positive.",
      call. = FALSE
    )
  }

  if (Tn == 1L) {
    return(r / D)
  }

  off2 <- off * off

  Dp <- matrix(0, nrow(D), Tn)
  rp <- matrix(0, nrow(D), Tn)
  Dp[, 1] <- D[, 1]
  rp[, 1] <- r[, 1]

  for (t in 2:Tn) {
    Dp[, t] <- D[, t] - off2 / Dp[, t - 1]
    rp[, t] <- r[, t] - off * rp[, t - 1] / Dp[, t - 1]
  }

  x <- matrix(0, nrow(D), Tn)
  x[, Tn] <- rp[, Tn] / Dp[, Tn]

  for (t in (Tn - 1):1) {
    x[, t] <- (rp[, t] - off * x[, t + 1]) / Dp[, t]
  }

  x
}


#' Selected inverse of many scalar tridiagonal systems at once
#'
#' Block size one specialisation of
#' `.block_tridiagonal_selected_inverse()`, vectorised across independent
#' systems in the same way as `.solve_scalar_tridiagonal_vec()`. Used for the
#' `a` and `b` blocks, which contribute one system per node.
#'
#' @param D `n x T` matrix of diagonal entries.
#' @param off Scalar on both off-diagonals.
#'
#' @return A list with `diag` (`n x T`, the entries `(H^{-1})_{tt}`) and
#'   `offdiag` (`n x (T-1)`, the entries `(H^{-1})_{t+1,t}`).
#' @keywords internal
#' @noRd
.scalar_tridiagonal_selected_inverse_vec <- function(D, off) {
  D <- as.matrix(D)
  n <- nrow(D)
  Tn <- ncol(D)

  Dp <- matrix(0, n, Tn)
  Dp[, 1] <- D[, 1]
  if (Tn > 1L) {
    off2 <- off * off
    for (t in 2:Tn) Dp[, t] <- D[, t] - off2 / Dp[, t - 1]
  }

  G <- matrix(0, n, Tn)
  G[, Tn] <- 1 / Dp[, Tn]

  if (Tn == 1L) {
    return(list(diag = G, offdiag = matrix(0, n, 0)))
  }

  off2 <- off * off
  Goff <- matrix(0, n, Tn - 1L)
  for (t in (Tn - 1):1) {
    G[, t] <- (1 + off2 * G[, t + 1] / Dp[, t]) / Dp[, t]
    Goff[, t] <- -G[, t + 1] * off / Dp[, t]
  }

  list(diag = G, offdiag = Goff)
}


#' Single common orthogonal transform aligning one latent set to another
#'
#' Finds the orthogonal matrix `R` (`K x K`) minimising
#' `sum_t || A_t R - B_t ||_F^2`. With `M = sum_t A_t' B_t = U S V'`, the
#' orthogonal Procrustes solution is `R = U V'`. No determinant constraint is
#' imposed, so reflections are allowed: `U_t V_t'` is invariant under any common
#' orthogonal transform of `(U_t, V_t)`, and permitting reflections gives the
#' better alignment.
#'
#' Used both for the Step-1 sequential alignment (one pair at a time) and for
#' Step-7 identification against the reference trajectory (all `T` at once).
#'
#' @param Alist List of matrices with `K` columns (the set being transformed).
#' @param Blist List of target matrices of matching shapes.
#'
#' @return A `K x K` orthogonal matrix `R`.
#' @keywords internal
#' @noRd
.orthogonal_procrustes <- function(Alist, Blist) {
  if (length(Alist) != length(Blist)) {
    stop("Alist and Blist must have the same length.", call. = FALSE)
  }
  if (length(Alist) == 0L) {
    stop("Alist must have at least one element.", call. = FALSE)
  }

  K <- ncol(Alist[[1]])
  M <- matrix(0, K, K)

  for (t in seq_along(Alist)) {
    if (ncol(Alist[[t]]) != K || ncol(Blist[[t]]) != K) {
      stop("All matrices must have the same number of columns.", call. = FALSE)
    }
    if (nrow(Alist[[t]]) != nrow(Blist[[t]])) {
      stop("Alist[[t]] and Blist[[t]] must have the same number of rows.", call. = FALSE)
    }
    M <- M + crossprod(Alist[[t]], Blist[[t]])   # A_t' B_t
  }

  s <- svd(M)
  s$u %*% t(s$v)
}


#' Whole-trajectory reciprocal scale normalisation of latent factors
#'
#' Step 8: fixes the multiplicative scale gauge of the latent factors. Since
#' `U_t V_t'` is invariant under `U_t -> c U_t`, `V_t -> V_t / c` for any
#' `c > 0`, the observation model cannot say how the magnitude of the
#' interaction is split between the two node sets — and an unchecked split lets
#' one factor shrink while the other grows across outer iterations, which would
#' corrupt the plug-in variance-component updates.
#'
#' The scale is pinned by the average squared coordinate magnitudes
#' \deqn{S_U = \frac{1}{NKT}\sum_t \|U_t\|_F^2, \qquad
#'       S_V = \frac{1}{MKT}\sum_t \|V_t\|_F^2,}
#' taking `c = (S_V / S_U)^{1/4}`, so that the normalised trajectories satisfy
#' `c^2 S_U = S_V / c^2`. Dividing by `N` and `M` respectively is what makes
#' this well behaved when the two node sets differ in size.
#'
#' The convention fixes only the relative unit of measurement: `sigma_U^2`,
#' `sigma_V^2`, `tau_U^2`, and `tau_V^2` remain distinct and are estimated
#' separately from the normalised trajectories.
#'
#' @section Which component ends up carrying the free coordinate:
#' One factor is applied to every period, so the increments are scaled
#' uniformly — `U_t - U_{t-1}` becomes `c (U_t - U_{t-1})` exactly — and the
#' convention cannot manufacture or erase temporal movement. What has to be
#' chosen is which magnitudes fix `c`, and that choice decides which of the
#' four latent variance components is left holding the undetermined direction.
#'
#' The rescaling carries `(sigma_U^2, tau_U^2, sigma_V^2, tau_V^2)` to
#' `(c^2 sigma_U^2, c^2 tau_U^2, c^-2 sigma_V^2, c^-2 tau_V^2)`, so only three
#' functions of the four are determined by the data. The one that carries the
#' asymmetry between the two node sets, and the only one interpretable without
#' reference to any convention, is the ratio of dimensionless drifts
#'
#'   delta = (tau_U^2 / sigma_U^2) / (tau_V^2 / sigma_V^2),
#'
#' in which `c^2` cancels within each side. It reads as how far the senders
#' move from period to period relative to their own dispersion, against the
#' same quantity for the receivers.
#'
#' \describe{
#'   \item{`"innovation"`}{Equalises the mean squared increments, so that
#'     `tau_U^2 = tau_V^2` after normalisation and `delta` is read off
#'     `sigma_V^2 / sigma_U^2`.}
#'   \item{`"pooled"`}{Equalises the magnitudes averaged over all periods,
#'     leaving `delta` in `tau_U^2 / tau_V^2`. A node set that genuinely moves
#'     more also spreads further, so the pooled magnitude picks the spread up
#'     and the factor absorbs part of the asymmetry as well as holding it.}
#'   \item{`"first"`}{Equalises the magnitudes of the initial states, which
#'     also leaves `delta` in `tau_U^2 / tau_V^2` but without that
#'     contamination.}
#' }
#'
#' The three are one model in different coordinates and `delta` is the same
#' number under all of them. They are not equally well behaved. The
#' empirical-Bayes loop derives each smoothing penalty from the innovation
#' variance it has just measured, `gamma_U = sigma_eps^2 / tau_U^2`, so a
#' smaller `tau_U^2` produces a heavier penalty on `U`'s increments, a smoother
#' `U`, and a smaller `tau_U^2` again. Any convention that leaves the free
#' coordinate in `tau_U^2 / tau_V^2` exposes it to that loop: on panels
#' generated with a true ratio of four the fitted ratio reached 68 and was
#' still rising when the loop stopped — it stopped because the fitted values
#' had settled, the drift running along a direction the convergence criterion
#' cannot see — and on others it fell below 1e-3, while the product
#' `tau_U^2 tau_V^2` was recovered to within twenty percent throughout.
#' `"innovation"` puts the free coordinate in `sigma_V^2 / sigma_U^2` instead,
#' which no penalty on the increments feeds back into.
#'
#' Two conventions cannot be imposed at once: there is one factor and one
#' degree of freedom. Equalising the increments and the magnitudes together
#' would drive `delta` to one and destroy the quantity being reported.
#'
#' The scalar `c` covers only the scale part of the gauge. `U_t V_t'` is also
#' invariant under `U_t A_t`, `V_t A_t^{-T}` for any invertible `A_t`, and
#' nothing here constrains the rest.
#'
#' If either magnitude is zero the latent block has collapsed; the scale is
#' left unchanged (`c = 1`) and a warning is issued rather than producing
#' `NaN`. With `T = 1` there are no increments, and `"innovation"` falls back
#' to `"first"` — on a single period the two are the same thing.
#'
#' @param Ulist List of length `T` of `N x K` latent-factor matrices.
#' @param Vlist List of length `T` of `M x K` latent-factor matrices.
#' @param anchor Which magnitudes fix the factor: the mean squared increments
#'   (`"innovation"`), the average over all periods (`"pooled"`), or the
#'   initial states alone (`"first"`).
#'
#' @return A list with the rescaled `U` and `V` lists, the scalar `c` applied,
#'   and the pre-normalisation magnitudes `S_U` and `S_V` it came from.
#' @keywords internal
#' @noRd
.scale_normalize_UV <- function(Ulist, Vlist,
                                anchor = c("pooled", "first", "innovation")) {
  anchor <- match.arg(anchor)
  if (length(Ulist) != length(Vlist)) {
    stop("Ulist and Vlist must have the same length.", call. = FALSE)
  }

  Tn <- length(Ulist)
  N <- nrow(Ulist[[1]])
  M <- nrow(Vlist[[1]])
  K <- ncol(Ulist[[1]])

  if (anchor == "innovation" && Tn == 1L) anchor <- "first"

  # Per-coordinate quantities, not raw totals: the sender and receiver node
  # sets differ in size, and it is the per-coordinate magnitudes the variance
  # components are defined on.
  if (anchor == "innovation") {
    S_U <- sum(vapply(2:Tn, function(t) sum((Ulist[[t]] - Ulist[[t - 1L]])^2),
                      numeric(1))) / (N * K * (Tn - 1L))
    S_V <- sum(vapply(2:Tn, function(t) sum((Vlist[[t]] - Vlist[[t - 1L]])^2),
                      numeric(1))) / (M * K * (Tn - 1L))
  } else {
    sq_U <- vapply(Ulist, function(U) sum(U^2), numeric(1)) / (N * K)
    sq_V <- vapply(Vlist, function(V) sum(V^2), numeric(1)) / (M * K)
    S_U <- if (anchor == "first") sq_U[1] else sum(sq_U) / Tn
    S_V <- if (anchor == "first") sq_V[1] else sum(sq_V) / Tn
  }

  if (S_U > 0 && S_V > 0) {
    # c^2 S_U = S_V / c^2 makes the two sides agree on whichever quantity the
    # anchor selected.
    cc <- (S_V / S_U)^(1 / 4)
  } else {
    warning(
      "Latent factors are degenerate (",
      if (anchor == "innovation") "no temporal movement in U or V" else
        "all zero in U or V",
      "); scale normalization skipped (c = 1). Consider a smaller K.",
      call. = FALSE
    )
    cc <- 1
  }

  list(
    U = lapply(Ulist, function(U) cc * U),
    V = lapply(Vlist, function(V) V / cc),
    c = cc,
    S_U = S_U,
    S_V = S_V
  )
}
