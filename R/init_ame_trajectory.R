# Step 1 of the joint-trajectory MAP estimator (fit_dynamic_ame): build the
# initial parameter trajectory, fix the Step-7 alignment reference, and generate
# the multistart candidates.

#' Initialise a parameter trajectory for the joint MAP estimator
#'
#' Implements Step 1 of the estimation procedure:
#' \enumerate{
#'   \item Per-period starting values `a_t^0 = rowMeans(Y_t)`,
#'     `b_t^0 = colMeans(Y_t) - mean(Y_t)` (all over observed cells) and
#'     `beta^0 = 0`.
#'   \item A balanced SVD of the residualised matrix
#'     `R_t = Y_t - a_t^0 - b_t^0`, giving `U_t^0 = P D^{1/2}`,
#'     `V_t^0 = Q D^{1/2}` so that `U_t^0 V_t^{0'}` is the rank-`K` truncation
#'     of `R_t`. Missing cells are set to zero in `R_t` before the SVD; they
#'     are starting values only and every later block update masks them out.
#'   \item A sequential orthogonal alignment for `t = 2, ..., T`: a single
#'     common transform is applied to `(U_t, V_t)` to bring them as close as
#'     possible to `(U_{t-1}, V_{t-1})`. This leaves each `U_t V_t'` unchanged
#'     while removing the arbitrary per-period rotation left by the SVD.
#'   \item Retention of the aligned trajectory as the reference
#'     `{U^ref, V^ref}` used later by the Step-7 identification.
#'   \item Generation of `n_starts` candidates: the unperturbed trajectory plus
#'     `n_starts - 1` copies with `N(0, perturb_sd^2)` noise added to `U` and
#'     `V` (the additive effects and regression coefficients are not perturbed;
#'     conditional on the latent factors they enter a convex problem).
#' }
#'
#' @param Y_list List of length `T` of `N x M` outcome matrices. Missing cells
#'   must be `NA`.
#' @param K Latent dimension.
#' @param P Length of the stacked regression coefficient vector. `beta^0` is a
#'   `P x T` matrix of zeros; use `0` when there are no covariates.
#' @param n_starts Number of multistart candidates (`>= 1`). The first is always
#'   the unperturbed trajectory.
#' @param perturb_sd Standard deviation of the perturbation applied to `U` and
#'   `V` in candidates `2, ..., n_starts`.
#' @param seed Random seed used for the perturbations. Defaults to a fixed value
#'   so that results replicate; pass another value to vary the starts. The
#'   caller's random-number state is restored on exit.
#'
#' @return A list with:
#'   \describe{
#'     \item{reference}{List with `U` and `V`, the aligned unperturbed latent
#'       trajectories, used as the Step-7 alignment target.}
#'     \item{candidates}{List of length `n_starts`. Each element has `a`
#'       (`N x T`), `b` (`M x T`), `beta` (`P x T`), and `U`, `V` (lists of
#'       length `T`).}
#'     \item{dims}{List with `N`, `M`, `K`, `P`, `T`.}
#'   }
#' @keywords internal
#' @noRd
.init_ame_trajectory <- function(Y_list,
                                 K = 2,
                                 P = 0,
                                 n_starts = 5,
                                 perturb_sd = 0.1,
                                 seed = 1) {
  Tn <- length(Y_list)
  if (Tn < 1L) {
    stop("Y_list must contain at least one period.", call. = FALSE)
  }
  if (length(K) != 1 || K < 1 || K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (length(n_starts) != 1 || n_starts < 1 || n_starts != as.integer(n_starts)) {
    stop("n_starts must be a positive integer.", call. = FALSE)
  }

  Y_list <- lapply(Y_list, function(Y) {
    Y <- as.matrix(Y)
    storage.mode(Y) <- "double"
    Y
  })

  N <- nrow(Y_list[[1]])
  M <- ncol(Y_list[[1]])

  for (t in seq_len(Tn)) {
    if (nrow(Y_list[[t]]) != N || ncol(Y_list[[t]]) != M) {
      stop("All matrices in Y_list must have the same dimensions.", call. = FALSE)
    }
  }
  if (K > min(N, M)) {
    stop(
      sprintf("K = %d exceeds min(N, M) = %d.", K, min(N, M)),
      call. = FALSE
    )
  }

  row_ids <- rownames(Y_list[[1]])
  col_ids <- colnames(Y_list[[1]])
  latent_ids <- paste0("k", seq_len(K))

  # Any statistic taken over an entirely missing row or column is undefined;
  # zero is the neutral starting value.
  zap <- function(x) {
    x[!is.finite(x)] <- 0
    x
  }

  a0 <- matrix(0, N, Tn)
  b0 <- matrix(0, M, Tn)
  U0 <- vector("list", Tn)
  V0 <- vector("list", Tn)

  # ── Per-period starting values and balanced SVD ──────────────────────────
  for (t in seq_len(Tn)) {
    Y <- Y_list[[t]]

    grand <- mean(Y, na.rm = TRUE)
    if (!is.finite(grand)) grand <- 0

    a_t <- zap(rowMeans(Y, na.rm = TRUE))
    b_t <- zap(colMeans(Y, na.rm = TRUE) - grand)

    a0[, t] <- a_t
    b0[, t] <- b_t

    # Residualise, then treat unobserved cells as zero for the SVD only.
    Rt <- Y - matrix(a_t, N, M, byrow = FALSE) - matrix(b_t, N, M, byrow = TRUE)
    Rt[is.na(Rt)] <- 0

    s <- svd(Rt, nu = K, nv = K)
    dK <- sqrt(s$d[seq_len(K)])
    Dh <- diag(dK, K, K)

    U0[[t]] <- s$u[, seq_len(K), drop = FALSE] %*% Dh
    V0[[t]] <- s$v[, seq_len(K), drop = FALSE] %*% Dh
  }

  # ── Sequential orthogonal alignment (t = 2, ..., T) ──────────────────────
  # One common transform per period, applied to U_t and V_t together, so that
  # U_t V_t' is unchanged. Chaining onto the already-aligned t-1 propagates a
  # single consistent frame along the whole trajectory.
  if (Tn > 1L) {
    for (t in 2:Tn) {
      Rt <- .orthogonal_procrustes(
        Alist = list(U0[[t]],     V0[[t]]),
        Blist = list(U0[[t - 1]], V0[[t - 1]])
      )
      U0[[t]] <- U0[[t]] %*% Rt
      V0[[t]] <- V0[[t]] %*% Rt
    }
  }

  for (t in seq_len(Tn)) {
    dimnames(U0[[t]]) <- list(row_ids, latent_ids)
    dimnames(V0[[t]]) <- list(col_ids, latent_ids)
  }

  rownames(a0) <- row_ids
  rownames(b0) <- col_ids
  colnames(a0) <- names(Y_list)
  colnames(b0) <- names(Y_list)

  beta0 <- matrix(0, P, Tn)

  # ── Multistart candidates ────────────────────────────────────────────────
  # Restore the caller's RNG stream so seeding here has no side effects.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(
        suppressWarnings(rm(".Random.seed", envir = globalenv())),
        add = TRUE
      )
    }
    set.seed(seed)
  }

  base_candidate <- list(a = a0, b = b0, beta = beta0, U = U0, V = V0)
  candidates <- vector("list", n_starts)
  candidates[[1]] <- base_candidate

  if (n_starts >= 2L) {
    for (s_i in 2:n_starts) {
      Us <- lapply(U0, function(U) {
        out <- U + matrix(stats::rnorm(N * K, sd = perturb_sd), N, K)
        dimnames(out) <- list(row_ids, latent_ids)
        out
      })
      Vs <- lapply(V0, function(V) {
        out <- V + matrix(stats::rnorm(M * K, sd = perturb_sd), M, K)
        dimnames(out) <- list(col_ids, latent_ids)
        out
      })
      candidates[[s_i]] <- list(a = a0, b = b0, beta = beta0, U = Us, V = Vs)
    }
  }

  list(
    reference = list(U = U0, V = V0),
    candidates = candidates,
    dims = list(
      N = as.integer(N), M = as.integer(M), K = as.integer(K),
      P = as.integer(P), T = as.integer(Tn)
    )
  )
}
