# Gibbs sampler for the dynamic bipartite AME model.
#
# The same model fit_dynamic_ame() maximises, sampled instead of optimised. It
# exists so the cost of the two inference paths can be compared on equal terms:
# same likelihood, same priors, same data, same machine, only the algorithm
# differing. A sampler for a *different* model would not answer that question —
# amen needs a square sociomatrix and has no temporal structure, and a
# general-purpose engine cannot exploit the conjugacy and mixes badly on the
# rotational ridge of u'v.
#
# Every conditional is closed form, so this is a plain Gibbs sampler with no
# Metropolis step and nothing to tune:
#
#   beta, a, b, U, V | rest   Gaussian, with the SAME block-tridiagonal system
#                             the block updates build. The optimiser takes that
#                             system's mean; the sampler draws from it. One
#                             sweep therefore costs about one inner iteration
#                             of the block coordinate descent.
#   the eleven variances      inverse-gamma, from the sums of squares of the
#                             *sampled* trajectories.
#
# Note the last point. The empirical-Bayes update in .eb_update_variances()
# adds a posterior-variance trace to the squared increments because it plugs in
# a point estimate that is itself smoothed. A Gibbs step conditions on a drawn
# trajectory rather than a smoothed one, so no trace correction belongs here;
# adding it would count the same dispersion twice.
#
# ── The data are taken already prepared ────────────────────────────────────
#
# The entry point takes `Y_list` and the covariate lists rather than an edge
# panel. Building those is shared work that belongs to neither inference path,
# and a timing comparison that charged it to one of them would be measuring the
# wrong thing. `gibbs_from_fit()` pulls them out of a fitted object.
#
# ── Identification ─────────────────────────────────────────────────────────
#
# `U_t V_t'` is unchanged by `(U_t A, V_t A^-T)` for constant invertible `A`,
# so the likelihood cannot see that direction. Two halves of it are handled
# differently:
#
#   the reciprocal scale   fixed after every sweep, by the same normalisation
#                          .ame_identify() applies to the estimator. Without
#                          it one factor grows while the other shrinks, and
#                          sigma_U^2 and tau_U^2 come out describing where the
#                          chain drifted rather than the data.
#   the rotation           left alone. The variance components are sums of
#                          squares over all coordinates and are already
#                          invariant to it, so fixing it would buy nothing
#                          that is reported.
#
# The columns of U and V therefore still wander, and their R-hat never
# approaches one. That remains a property of the parameter and not of the
# sampler, and diagnosing on them would report a failure that is not there.
# Monitor `beta`, the variance components, or the products `U_t V_t'` — all
# identified. The draws kept by default are `beta` and the variances.


#' Draw from a Gaussian with block-tridiagonal precision
#'
#' Returns a draw from `N(H^-1 r, scale * H^-1)`, where `H` is the symmetric
#' block-tridiagonal matrix with diagonal blocks `D` and both off-diagonals
#' `off * I` — the system the block updates solve.
#'
#' The mean comes from the existing solver. The dispersion comes from the LDL'
#' factorisation `H = L D' L'` whose pivots the existing forward sweep already
#' returns. With `L[t, t-1] = off * D'[t-1]^-1`, drawing `z ~ N(0, I)`,
#' whitening by `chol(D'[t])` and back-substituting through `L'` gives
#' covariance `L^-T D'^-1 L^-1 = H^-1`, at `O(T q^3)` — the same order as the
#' solve, which is why a sweep costs about what an inner iteration costs.
#'
#' @param D List of `T` diagonal blocks (`q x q`).
#' @param off Scalar on both off-diagonals.
#' @param r List of `T` right-hand sides of length `q`.
#' @param scale Multiplies the covariance. `sigma_eps^2` here: `H` is built in
#'   the units the objective uses, where the precision is `H / sigma_eps^2`.
#' @return A `q x T` matrix.
#' @keywords internal
#' @noRd
.rtridiag_block <- function(D, off, r, scale) {
  Tn <- length(D)
  q <- nrow(as.matrix(D[[1]]))

  # the solver returns one vector per period; the draw is assembled as q x T
  mu <- matrix(unlist(.solve_block_tridiagonal(D, off, r)), q, Tn)
  Dp <- .block_tridiagonal_pivots(D, off)

  w <- matrix(0, q, Tn)
  for (t in seq_len(Tn)) {
    Ct <- chol(Dp[[t]])                       # upper triangular, Ct'Ct = D'[t]
    w[, t] <- backsolve(Ct, stats::rnorm(q))  # covariance D'[t]^-1
  }

  th <- matrix(0, q, Tn)
  th[, Tn] <- w[, Tn]
  if (Tn > 1L)
    for (t in (Tn - 1L):1L)
      th[, t] <- w[, t] - off * solve(Dp[[t]], th[, t + 1L])

  mu + sqrt(scale) * th
}


#' Draw many scalar tridiagonal systems at once
#'
#' Block size one, vectorised across independent systems, matching
#' `.solve_scalar_tridiagonal_vec()`. The additive blocks give one system per
#' node, all sharing the same off-diagonal, so the recursion runs on columns.
#'
#' @param D `n x T` matrix of diagonal entries.
#' @param off Scalar on both off-diagonals.
#' @param r `n x T` matrix of right-hand sides.
#' @param scale Multiplies the covariance.
#' @return An `n x T` matrix.
#' @keywords internal
#' @noRd
.rtridiag_scalar_vec <- function(D, off, r, scale) {
  D <- as.matrix(D); r <- as.matrix(r)
  n <- nrow(D); Tn <- ncol(D)

  mu <- .solve_scalar_tridiagonal_vec(D, off, r)

  Dp <- matrix(0, n, Tn)
  Dp[, 1] <- D[, 1]
  if (Tn > 1L)
    for (t in 2:Tn) Dp[, t] <- D[, t] - off * off / Dp[, t - 1L]

  w <- matrix(stats::rnorm(n * Tn), n, Tn) / sqrt(Dp)

  th <- matrix(0, n, Tn)
  th[, Tn] <- w[, Tn]
  if (Tn > 1L)
    for (t in (Tn - 1L):1L)
      th[, t] <- w[, t] - off * th[, t + 1L] / Dp[, t]

  mu + sqrt(scale) * th
}


#' Partial residual with one block removed
#'
#' @param params Current parameters.
#' @param Y_list Outcome matrices.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param t Period.
#' @param drop Block to leave out: `"beta"`, `"a"`, `"b"`, `"latent"`, or
#'   `"none"` for the full residual.
#' @return An `N x M` matrix, `NA` where the outcome is unobserved.
#' @keywords internal
#' @noRd
.partial_residual <- function(params, Y_list, X_row_list, X_col_list,
                              X_dyad_list, t, drop) {
  N <- nrow(params$a); M <- nrow(params$b)
  R <- Y_list[[t]]
  if (drop != "beta" && nrow(params$beta) > 0L)
    R <- R - .cov_term(params$beta[, t],
                       if (is.null(X_row_list))  NULL else X_row_list[[t]],
                       if (is.null(X_col_list))  NULL else X_col_list[[t]],
                       if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]],
                       N, M)
  if (drop != "a")      R <- R - matrix(params$a[, t], N, M, byrow = FALSE)
  if (drop != "b")      R <- R - matrix(params$b[, t], N, M, byrow = TRUE)
  if (drop != "latent") R <- R - params$U[[t]] %*% t(params$V[[t]])
  R
}


#' One Gibbs sweep over the five trajectory blocks
#'
#' Each block is drawn from its full conditional given the current value of
#' every other one, in the order the optimiser visits them.
#'
#' @param params Current parameters.
#' @param Y_list,X_row_list,X_col_list,X_dyad_list Data.
#' @param lambda,gamma Named penalties over the five blocks.
#' @param sigma_eps2 Current observation variance.
#' @return `params` with every block replaced by a draw.
#' @keywords internal
#' @noRd
.gibbs_trajectories <- function(params, Y_list, X_row_list, X_col_list,
                                X_dyad_list, lambda, gamma, sigma_eps2) {
  Tn <- length(Y_list)
  N <- nrow(params$a); M <- nrow(params$b)
  P <- nrow(params$beta); K <- ncol(params$U[[1]])
  nb <- (seq_len(Tn) > 1L) + (seq_len(Tn) < Tn)

  # ── beta: one P x P system shared by every cell ─────────────────────────
  if (P > 0L) {
    Ip <- diag(P)
    D <- r <- vector("list", Tn)
    for (t in seq_len(Tn)) {
      R <- .partial_residual(params, Y_list, X_row_list, X_col_list,
                             X_dyad_list, t, "beta")
      Z <- .cov_design(if (is.null(X_row_list))  NULL else X_row_list[[t]],
                       if (is.null(X_col_list))  NULL else X_col_list[[t]],
                       if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]],
                       N, M)
      y <- as.vector(R)
      keep <- !is.na(y) & stats::complete.cases(Z)
      Zk <- Z[keep, , drop = FALSE]
      D[[t]] <- crossprod(Zk) +
        (gamma[["beta"]] * nb[t] + if (t == 1L) lambda[["beta"]] else 0) * Ip
      r[[t]] <- as.vector(crossprod(Zk, y[keep]))
    }
    bd <- dimnames(params$beta)
    params$beta <- .rtridiag_block(D, -gamma[["beta"]], r, sigma_eps2)
    dimnames(params$beta) <- bd
  }

  # ── a and b: one scalar system per node ─────────────────────────────────
  for (blk in c("a", "b")) {
    n_row <- if (blk == "a") N else M
    n_obs <- matrix(0, n_row, Tn)
    rhs   <- matrix(0, n_row, Tn)
    for (t in seq_len(Tn)) {
      R <- .partial_residual(params, Y_list, X_row_list, X_col_list,
                             X_dyad_list, t, blk)
      if (blk == "b") R <- t(R)
      n_obs[, t] <- rowSums(!is.na(R))
      rhs[, t]   <- rowSums(R, na.rm = TRUE)
    }
    Dm <- .additive_diagonal(n_obs, lambda[[blk]], gamma[[blk]])
    dn <- dimnames(params[[blk]])
    params[[blk]] <- .rtridiag_scalar_vec(Dm, -gamma[[blk]], rhs, sigma_eps2)
    dimnames(params[[blk]]) <- dn
  }

  # ── U and V: one K x K system per node ──────────────────────────────────
  # The Gram matrices come from the same builder the block update and the
  # variance step use, so there is one source of truth for them.
  Ik <- diag(K)
  for (blk in c("U", "V")) {
    n_sys <- if (blk == "U") N else M
    Gt <- ct <- vector("list", Tn)
    for (t in seq_len(Tn)) {
      R <- .partial_residual(params, Y_list, X_row_list, X_col_list,
                             X_dyad_list, t, "latent")
      gr <- if (blk == "U") .latent_gram_rhs(R, params$V[[t]]) else
        .latent_gram_rhs(t(R), params$U[[t]])
      Gt[[t]] <- gr$G; ct[[t]] <- gr$c
    }
    new <- lapply(seq_len(Tn), function(t) matrix(0, n_sys, K))
    for (i in seq_len(n_sys)) {
      D <- r <- vector("list", Tn)
      for (t in seq_len(Tn)) {
        D[[t]] <- matrix(Gt[[t]][i, , ], K, K) +
          (gamma[[blk]] * nb[t] + if (t == 1L) lambda[[blk]] else 0) * Ik
        r[[t]] <- ct[[t]][i, ]
      }
      th <- .rtridiag_block(D, -gamma[[blk]], r, sigma_eps2)
      for (t in seq_len(Tn)) new[[t]][i, ] <- th[, t]
    }
    for (t in seq_len(Tn)) {
      dimnames(new[[t]]) <- dimnames(params[[blk]][[t]])
      params[[blk]][[t]] <- new[[t]]
    }
  }

  # ── fix the reciprocal scale gauge ──────────────────────────────────────
  #
  # U_t V_t' is unchanged by (c U_t, V_t / c), so the likelihood cannot see
  # that direction at all and the prior barely constrains it. Left alone the
  # pair random-walks along it: one factor grows while the other shrinks, the
  # product staying exactly put. Measured on a 40 x 50 x 8 panel, four chains
  # ended at sigma_U^2 of 2.83, 1.76, 0.37 and 1.23 while their products
  # agreed to three digits.
  #
  # Nothing identified is harmed by that -- beta, the fitted values and
  # sigma_U^2 * sigma_V^2 were the same in all four -- but sigma_U^2,
  # sigma_V^2, tau_U^2 and tau_V^2 are reported quantities, and a number that
  # differs eightfold between runs of the same data describes where the
  # sampler drifted rather than what the data say. It is also not comparable
  # with the estimator's, which are reported in the normalised gauge.
  #
  # So the same convention is applied here, and for the same reason
  # .ame_identify() applies it there. Only the scale needs fixing: the
  # variance components are sums of squares over all coordinates and are
  # already invariant to the rotation, which is the other half of the gauge
  # freedom and is left alone.
  # Pinned explicitly rather than left to the default, which the estimator has
  # since moved to the initial states. Whether the sampler should follow is a
  # question about the model the two paths are meant to share, and it is not
  # settled by the estimator's reasons alone: the estimator amplifies the split
  # because it derives its penalties from the split it has just measured,
  # whereas a Gibbs step draws the variances from a full conditional and has no
  # such loop. Until the sampler is shown to drift, it keeps the convention its
  # calibration and timings were run under.
  sn <- .scale_normalize_UV(params$U, params$V, anchor = "pooled")
  for (t in seq_along(params$U)) {
    dimnames(sn$U[[t]]) <- dimnames(params$U[[t]])
    dimnames(sn$V[[t]]) <- dimnames(params$V[[t]])
  }
  params$U <- sn$U
  params$V <- sn$V

  params
}


#' Draw the eleven variance components
#'
#' Inverse-gamma throughout, from the sums of squares of the *sampled*
#' trajectories. Conditioning on a draw rather than on a smoothed point
#' estimate is what removes the posterior-variance trace the empirical-Bayes
#' update needs.
#'
#' @param params Current (just drawn) parameters.
#' @param Y_list,X_row_list,X_col_list,X_dyad_list Data.
#' @param prior List with `a0` and `b0`.
#' @return A named vector of eleven variances, in `Omega`'s order.
#' @keywords internal
#' @noRd
.gibbs_omega <- function(params, Y_list, X_row_list, X_col_list, X_dyad_list,
                         prior, tie_latent = TRUE) {
  Tn <- length(Y_list)
  rig <- function(shape, rate) 1 / stats::rgamma(1L, shape = shape, rate = rate)
  # The latent innovation sums, kept aside so they can be pooled into one draw.
  lat_n <- c(U = 0, V = 0)
  lat_d <- c(U = 0, V = 0)

  ssr <- 0; n_obs <- 0L
  for (t in seq_len(Tn)) {
    e <- .partial_residual(params, Y_list, X_row_list, X_col_list,
                           X_dyad_list, t, "none")
    ssr   <- ssr + sum(e^2, na.rm = TRUE)
    n_obs <- n_obs + sum(!is.na(e))
  }

  out <- c(sigma_eps2 = rig(prior$a0 + n_obs / 2, prior$b0 + ssr / 2))

  traj <- list(beta = params$beta, a = params$a, b = params$b,
               U = params$U, V = params$V)
  for (nm in c("beta", "a", "b", "U", "V")) {
    th <- traj[[nm]]
    flat <- if (is.list(th)) th else
      lapply(seq_len(Tn), function(t) th[, t, drop = TRUE])
    n_scalar <- length(flat[[1]])

    if (n_scalar == 0L) {
      out[paste0("sigma_", nm, "2")] <- NA_real_
      out[paste0("tau_", nm, "2")]   <- NA_real_
      next
    }

    out[paste0("sigma_", nm, "2")] <-
      rig(prior$a0 + n_scalar / 2, prior$b0 + sum(flat[[1]]^2) / 2)

    if (Tn > 1L) {
      d <- 0
      for (t in 2:Tn) d <- d + sum((flat[[t]] - flat[[t - 1L]])^2)
      if (nm %in% c("U", "V")) {
        lat_n[nm] <- n_scalar * (Tn - 1L)
        lat_d[nm] <- d
      }
      out[paste0("tau_", nm, "2")] <-
        rig(prior$a0 + n_scalar * (Tn - 1L) / 2, prior$b0 + d / 2)
    } else {
      out[paste0("tau_", nm, "2")] <- NA_real_
    }
  }

  # ── One innovation variance for the latent block ─────────────────────────
  #
  # The estimator constrains these to be equal, only their product being
  # determined by the data, and a timing comparison between the two inference
  # paths is a comparison of nothing unless both are fitting the same model.
  #
  # The sampler does not average two draws: under a shared kappa the full
  # conditional is a single inverse-gamma on the pooled sufficient statistics,
  # which is what the model implies and what keeps the chain a correct sampler
  # rather than an approximation of one.
  #
  # The estimator's reason for the constraint does not apply here — there is no
  # plug-in feedback in a Gibbs step, each variance being drawn from its full
  # conditional rather than derived from the point it just produced — so this is
  # imposed for comparability, not because the sampler was shown to drift.
  if (tie_latent && Tn > 1L && all(lat_n > 0)) {
    kappa <- rig(prior$a0 + sum(lat_n) / 2, prior$b0 + sum(lat_d) / 2)
    out["tau_U2"] <- kappa
    out["tau_V2"] <- kappa
  }

  out[c("sigma_eps2", "sigma_beta2", "sigma_a2", "sigma_b2", "sigma_U2",
        "sigma_V2", "tau_beta2", "tau_a2", "tau_b2", "tau_U2", "tau_V2")]
}


#' Gibbs sampler for a dynamic bipartite AME model
#'
#' Samples the posterior of the model [fit_dynamic_ame()] maximises. Provided
#' for comparing the cost of the two inference paths on equal terms, and for
#' users who want full posterior inference; the package's primary estimator
#' remains [fit_dynamic_ame()].
#'
#' Every conditional is closed form, so there is nothing to tune and no
#' acceptance rate to watch. A sweep solves the same block-tridiagonal systems
#' the optimiser solves, drawing from each rather than taking its mean, and so
#' costs about one inner iteration of the block coordinate descent.
#'
#' The data arrive already prepared. Building them from an edge panel is shared
#' work belonging to neither inference path, and charging it to one of them
#' would measure the wrong thing; [gibbs_from_fit()] takes them from a fitted
#' object.
#'
#' `U_t V_t'` is unchanged by `(U_t A, V_t A^-T)` for constant invertible `A`.
#' The reciprocal-scale half of that freedom is fixed after every sweep, by the
#' same normalisation the estimator applies: left alone one factor grows while
#' the other shrinks with the product exactly fixed, and `sigma_U^2` and
#' `tau_U^2` end up describing where the chain drifted rather than the data.
#' The rotation is left alone, the variance components being invariant to it.
#' The columns of `U` and `V` therefore still wander and their convergence
#' diagnostics never settle, which is a property of the parameter rather than
#' of the sampler; monitor `beta`, the variance components, or the products
#' `U_t V_t'`, all of which are identified.
#'
#' @param Y_list List of `T` outcome matrices, `NA` where unobserved.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists, or `NULL`.
#' @param K Latent dimension.
#' @param n_iter Sweeps to keep after `burn`.
#' @param burn Sweeps to discard.
#' @param thin Keep every `thin`-th sweep.
#' @param store_latent Also keep `U_t V_t'` for every kept sweep. Memory grows
#'   as `n_iter/thin * T * N * M`; off by default.
#' @param prior Inverse-gamma shape `a0` and rate `b0` for every variance.
#' @param latent_innovation Whether the latent block has one innovation
#'   variance or two, matching the argument of the same name in
#'   [fit_dynamic_ame()]. Under `"shared"` a single `kappa` is drawn from the
#'   pooled sufficient statistics of both sides — the full conditional the model
#'   implies, not an average of two draws. The default matches the estimator's,
#'   because a comparison between the two inference paths is a comparison of
#'   nothing unless both are fitting the same model. `"separate"` draws them
#'   independently, as earlier versions did.
#' @param init Optional starting parameters, as returned in `$params`. Chains
#'   of a diagnostic run must start far apart or their agreement carries no
#'   information; when `NULL` the start is drawn from a diffuse Gaussian.
#' @param covariate_names Row names for `beta`; taken from `X_dyad_list` when
#'   available.
#' @param years Period labels used to name the stored columns.
#' @param seed Random seed.
#' @param verbose Report progress every 500 sweeps.
#'
#' @return An object of class `gibbs_dynamic_ame`: `beta` (kept sweeps by
#'   `P * T`, columns named `covariate:period`), `omega` (kept sweeps by 11),
#'   `latent` when requested, `params` (last state, for restarting), `settings`,
#'   and `seconds`.
#'
#' @examples
#' \dontrun{
#' fit <- fit_dynamic_ame(edge_panel = my_panel, K = 2)
#' g <- gibbs_from_fit(fit, n_iter = 5000, burn = 1000)
#' summary(g)
#' }
#' @export
gibbs_dynamic_ame <- function(Y_list,
                              X_row_list = NULL,
                              X_col_list = NULL,
                              X_dyad_list = NULL,
                              K = 2,
                              n_iter = 5000,
                              burn = 1000,
                              thin = 1,
                              store_latent = FALSE,
                              prior = list(a0 = 1e-3, b0 = 1e-3),
                              latent_innovation = c("shared", "separate"),
                              init = NULL,
                              covariate_names = NULL,
                              years = NULL,
                              seed = 1,
                              verbose = FALSE) {

  latent_innovation <- match.arg(latent_innovation)
  t_start <- Sys.time()

  if (!is.list(Y_list) || !length(Y_list))
    stop("Y_list must be a non-empty list of outcome matrices.", call. = FALSE)
  if (n_iter < 1L || burn < 0L || thin < 1L)
    stop("n_iter must be positive, burn non-negative, thin at least one.",
         call. = FALSE)

  Tn <- length(Y_list)
  N <- nrow(Y_list[[1]]); M <- ncol(Y_list[[1]])
  if (is.null(years)) years <- seq_len(Tn)

  pdims <- .cov_dims(
    X_row  = if (is.null(X_row_list))  NULL else X_row_list[[1]],
    X_col  = if (is.null(X_col_list))  NULL else X_col_list[[1]],
    X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[1]])
  P <- pdims$P
  if (is.null(covariate_names) && P > 0L)
    covariate_names <- paste0("x", seq_len(P))

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)

  rn <- rownames(Y_list[[1]]); cn <- colnames(Y_list[[1]])
  params <- if (!is.null(init)) init else list(
    beta = matrix(stats::rnorm(P * Tn), P, Tn,
                  dimnames = list(covariate_names, NULL)),
    a = matrix(stats::rnorm(N * Tn), N, Tn, dimnames = list(rn, NULL)),
    b = matrix(stats::rnorm(M * Tn), M, Tn, dimnames = list(cn, NULL)),
    U = lapply(seq_len(Tn), function(t)
      matrix(stats::rnorm(N * K), N, K, dimnames = list(rn, NULL))),
    V = lapply(seq_len(Tn), function(t)
      matrix(stats::rnorm(M * K), M, K, dimnames = list(cn, NULL))))

  blocks <- c("beta", "a", "b", "U", "V")
  om_names <- c("sigma_eps2", "sigma_beta2", "sigma_a2", "sigma_b2",
                "sigma_U2", "sigma_V2", "tau_beta2", "tau_a2", "tau_b2",
                "tau_U2", "tau_V2")
  Om <- stats::setNames(rep(1, 11), om_names)

  n_keep <- floor(n_iter / thin)
  beta_draws <- if (P > 0L) matrix(NA_real_, n_keep, P * Tn) else NULL
  if (P > 0L)
    colnames(beta_draws) <- paste(rep(covariate_names, times = Tn),
                                  rep(years, each = P), sep = ":")
  omega_draws <- matrix(NA_real_, n_keep, 11L,
                        dimnames = list(NULL, om_names))
  latent_draws <- if (store_latent) vector("list", n_keep) else NULL

  kept <- 0L
  for (it in seq_len(burn + n_iter)) {

    # penalties from the current variance draw, derived exactly as the
    # optimiser derives them: lambda = sigma_eps^2 / sigma^2,
    # gamma = sigma_eps^2 / tau^2
    lam <- stats::setNames(vapply(blocks, function(nm) {
      s2 <- Om[[paste0("sigma_", nm, "2")]]
      if (is.na(s2)) 1 else Om[["sigma_eps2"]] / s2
    }, numeric(1)), blocks)
    gam <- stats::setNames(vapply(blocks, function(nm) {
      t2 <- Om[[paste0("tau_", nm, "2")]]
      if (is.na(t2)) 1 else Om[["sigma_eps2"]] / t2
    }, numeric(1)), blocks)

    params <- .gibbs_trajectories(params, Y_list, X_row_list, X_col_list,
                                  X_dyad_list, lam, gam, Om[["sigma_eps2"]])
    Om <- .gibbs_omega(params, Y_list, X_row_list, X_col_list, X_dyad_list,
                       prior, tie_latent = (latent_innovation == "shared"))

    if (it > burn && ((it - burn) %% thin == 0L)) {
      kept <- kept + 1L
      if (P > 0L) beta_draws[kept, ] <- as.vector(params$beta)
      omega_draws[kept, ] <- Om
      if (store_latent)
        latent_draws[[kept]] <- lapply(seq_len(Tn), function(t)
          params$U[[t]] %*% t(params$V[[t]]))
    }

    if (verbose && it %% 500L == 0L) {
      cat(sprintf("  sweep %d / %d\n", it, burn + n_iter))
      utils::flush.console()
    }
  }

  structure(
    list(beta = beta_draws, omega = omega_draws, latent = latent_draws,
         params = params, years = years, row_ids = rn, col_ids = cn,
         covariate_names = covariate_names,
         settings = list(K = K, n_iter = n_iter, burn = burn, thin = thin,
                         prior = prior, seed = seed, n_kept = kept),
         seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))),
    class = "gibbs_dynamic_ame")
}


#' Run the Gibbs sampler on the data behind a fitted model
#'
#' Convenience wrapper: takes the prepared outcome and covariate arrays out of
#' a `dynamic_ame` object so the sampler sees exactly the panel the estimator
#' saw. Nothing about the fit itself is used, so the two remain independent
#' inference paths over the same data.
#'
#' @param fit A `dynamic_ame` object.
#' @param ... Passed to [gibbs_dynamic_ame()].
#' @return An object of class `gibbs_dynamic_ame`.
#' @export
gibbs_from_fit <- function(fit, ...) {
  .validate_dynamic_ame(fit)
  cv <- fit$covariates
  gibbs_dynamic_ame(
    Y_list = fit$Y_list,
    X_row_list = cv$X_row_list, X_col_list = cv$X_col_list,
    X_dyad_list = cv$X_dyad_list,
    K = ncol(fit$U[[1]]),
    covariate_names = rownames(fit$beta),
    years = fit$years, ...)
}


#' Posterior summaries and within-chain diagnostics
#'
#' Quantiles for every monitored parameter, with the effective sample size the
#' chain reached. Both are reported: `ess_bulk` governs the precision of a
#' posterior mean, `ess_tail` that of the interval endpoints, and since what is
#' reported here is a 95% interval, `ess_tail` is the one that matters.
#'
#' R-hat needs several chains started far apart and is therefore not computed
#' here; see [gibbs_rhat()].
#'
#' @param object A `gibbs_dynamic_ame` object.
#' @param pars Which draws to summarise: `"beta"`, `"omega"`, or `"both"`.
#' @param conf_level Interval level.
#' @param ... Ignored.
#' @return A data frame, one row per parameter.
#' @export
summary.gibbs_dynamic_ame <- function(object, pars = c("beta", "omega", "both"),
                                      conf_level = 0.95, ...) {
  pars <- match.arg(pars)
  mats <- list()
  if (pars %in% c("beta", "both") && !is.null(object$beta))
    mats$beta <- object$beta
  if (pars %in% c("omega", "both")) mats$omega <- object$omega
  if (!length(mats)) stop("nothing to summarise.", call. = FALSE)

  probs <- c((1 - conf_level) / 2, 0.5, 1 - (1 - conf_level) / 2)
  have_posterior <- requireNamespace("posterior", quietly = TRUE)

  out <- lapply(names(mats), function(blk) {
    Mx <- mats[[blk]]
    keep <- vapply(seq_len(ncol(Mx)),
                   function(j) any(!is.na(Mx[, j])), logical(1))
    Mx <- Mx[, keep, drop = FALSE]
    q <- t(apply(Mx, 2, stats::quantile, probs = probs, na.rm = TRUE))
    data.frame(
      block = blk, parameter = colnames(Mx),
      mean = colMeans(Mx, na.rm = TRUE), sd = apply(Mx, 2, stats::sd, na.rm = TRUE),
      conf.low = q[, 1], median = q[, 2], conf.high = q[, 3],
      ess_bulk = if (have_posterior)
        apply(Mx, 2, posterior::ess_bulk) else NA_real_,
      ess_tail = if (have_posterior)
        apply(Mx, 2, posterior::ess_tail) else NA_real_,
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


#' Rank-normalised split R-hat across chains
#'
#' Compares the variance within chains with the variance between them, which is
#' informative only when the chains started far apart: chains launched from the
#' same point agree immediately whether or not they have found the target, and
#' R-hat would report convergence that has not happened.
#'
#' Computed on `beta` by default. `U` and `V` are identified only up to a
#' rotation, so their chains never agree and their R-hat never settles — a
#' property of the parameter rather than a failure to converge, and the reason
#' they are not monitored.
#'
#' The threshold in current use is 1.01 (Vehtari et al. 2021), not the older
#' 1.1. It is also the value consistent with an effective sample size in the
#' hundreds: at R-hat 1.01 the residual between-chain discrepancy is about
#' 0.5% of a posterior standard deviation, an order of magnitude below the
#' Monte Carlo error already being accepted at ESS 500, whereas at 1.1 the two
#' would be the same size and the diagnostic would carry no information.
#'
#' @param chains List of `gibbs_dynamic_ame` objects from dispersed starts.
#' @param pars `"beta"` or `"omega"`.
#' @return A data frame with `parameter`, `rhat`, `ess_bulk`, `ess_tail`.
#' @export
gibbs_rhat <- function(chains, pars = c("beta", "omega")) {
  pars <- match.arg(pars)
  if (!is.list(chains) || length(chains) < 2L)
    stop("at least two chains are needed; one chain cannot be compared with ",
         "itself.", call. = FALSE)
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("package 'posterior' is required for R-hat.", call. = FALSE)

  mats <- lapply(chains, function(ch) ch[[pars]])
  if (any(vapply(mats, is.null, logical(1))))
    stop("some chains have no '", pars, "' draws.", call. = FALSE)

  n_par <- ncol(mats[[1]])
  if (!all(vapply(mats, function(Mx) ncol(Mx) == n_par, logical(1))))
    stop("chains do not share the same parameters.", call. = FALSE)
  nm <- colnames(mats[[1]])
  if (is.null(nm)) nm <- paste0("par", seq_len(n_par))
  if (!all(vapply(mats, function(Mx)
    is.null(colnames(Mx)) || identical(colnames(Mx), nm), logical(1))))
    stop("chains do not share the same parameters.", call. = FALSE)

  n_it <- min(vapply(mats, nrow, integer(1)))
  keep <- which(vapply(seq_len(n_par),
                       function(j) all(vapply(mats, function(Mx)
                         any(!is.na(Mx[, j])), logical(1))), logical(1)))

  # An empty result would be read by any caller taking max() over it as a
  # perfect R-hat rather than as no R-hat at all, and a convergence check
  # that passes because it computed nothing is worse than one that fails.
  if (!length(keep))
    stop("no parameter has usable draws in every chain: R-hat is undefined ",
         "here, and returning nothing would look like convergence.",
         call. = FALSE)

  res <- lapply(keep, function(j) {
    arr <- vapply(mats, function(Mx) Mx[seq_len(n_it), j], numeric(n_it))
    data.frame(parameter = nm[j],
               rhat = posterior::rhat(arr),
               ess_bulk = posterior::ess_bulk(arr),
               ess_tail = posterior::ess_tail(arr),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}


#' Print a Gibbs sampler result
#'
#' @param x A `gibbs_dynamic_ame` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.gibbs_dynamic_ame <- function(x, ...) {
  cat("Gibbs sampler for a dynamic bipartite AME model\n")
  cat(sprintf("Sweeps: %d kept after %d discarded, thin = %d  (%.1f s)\n",
              x$settings$n_kept, x$settings$burn, x$settings$thin, x$seconds))
  cat(sprintf("Panel: %d senders, %d receivers, %d periods, K = %d\n",
              length(x$row_ids), length(x$col_ids), length(x$years),
              x$settings$K))

  if (!is.null(x$beta)) {
    s <- summary(x, pars = "beta")
    cat("\nCoefficient posterior (first six):\n")
    print(utils::head(s[, c("parameter", "mean", "conf.low", "conf.high",
                            "ess_tail")], 6L), row.names = FALSE, digits = 4)
    if (nrow(s) > 6L) cat("  ...\n")
  }

  cat("\nU and V are identified only up to a rotation and a reciprocal\n")
  cat("rescaling; monitor beta, the variances, or U_t V_t' -- not U or V.\n")
  invisible(x)
}
