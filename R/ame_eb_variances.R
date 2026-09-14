# Step 9: empirical Bayes update of the variance components, and the penalties
# derived from them.
#
# This is the M-step of an EM algorithm, not a plug-in. For a variance component
# governing a block theta, the M-step maximises the expected complete-data
# log-likelihood, which requires the conditional expectation of a sum of squares
# rather than the sum of squares of the point estimates:
#
#   E[ ||theta_1||^2 | Y ]            = ||theta_1_hat||^2 + tr(Sigma_11)
#   E[ ||theta_t - theta_{t-1}||^2 ]  = ||d theta_hat||^2
#                                       + tr(Sigma_tt) + tr(Sigma_{t-1,t-1})
#                                       - 2 tr(Sigma_{t,t-1})
#
# where Sigma = Var(theta | Y). The objective Q is the penalised least squares
# criterion, so the negative log posterior is Q / (2 sigma_eps^2) and its Hessian
# is H / sigma_eps^2, giving Sigma = sigma_eps^2 H^{-1} with H the same
# block-tridiagonal matrix the block updates solve. Only the diagonal and first
# sub-diagonal blocks of H^{-1} are needed, and those come from the selected
# inversion at O(T q^3) — the same order as the solve.
#
# Dropping the trace terms is what makes a plug-in update collapse: the point
# estimates are already smoothed by the penalty, so their successive differences
# understate the true innovation, which raises the penalty, which smooths them
# further. The trace terms grow exactly as the data stop pinning the trajectory
# down, and hold the estimate up.
#
# Note that sigma_eps^2 itself remains a plug-in estimate: it has a known
# downward bias because the fitted parameters absorb some of the noise.

#' Number of scalar parameters in one period of each block
#'
#' @param params Parameter list.
#' @return A named integer vector over `beta`, `a`, `b`, `U`, `V`.
#' @keywords internal
#' @noRd
.block_scalar_counts <- function(params) {
  K <- ncol(params$U[[1]])
  c(
    beta = nrow(params$beta),
    a    = nrow(params$a),
    b    = nrow(params$b),
    U    = nrow(params$U[[1]]) * K,
    V    = nrow(params$V[[1]]) * K
  )
}


#' Trace corrections for one block, summed over its independent systems
#'
#' Given the selected inverse of each system's block-tridiagonal matrix, returns
#' the two totals the M-step needs: the trace attached to the initial state and
#' the trace attached to the innovations, both summed over systems (and over
#' periods, for the innovations).
#'
#' @param sel_list List of selected-inverse results, one per system, each with
#'   `diag` and `offdiag`.
#'
#' @return A list with `init` and `innov`.
#' @keywords internal
#' @noRd
.trace_totals <- function(sel_list) {
  init <- 0
  innov <- 0
  for (s in sel_list) {
    init <- init + sum(diag(as.matrix(s$diag[[1]])))
    if (length(s$offdiag) > 0L) {
      for (t in seq_along(s$offdiag)) {
        innov <- innov +
          sum(diag(as.matrix(s$diag[[t + 1]]))) +
          sum(diag(as.matrix(s$diag[[t]]))) -
          2 * sum(diag(as.matrix(s$offdiag[[t]])))
      }
    }
  }
  list(init = init, innov = innov)
}


#' Diagonal blocks of the system solved by each additive-block node
#'
#' @keywords internal
#' @noRd
.additive_systems <- function(n_obs, lambda, gamma) {
  .additive_diagonal(n_obs, lambda, gamma)
}


#' Per-node observation counts and latent Grams at the current parameters
#'
#' Rebuilds the information each block update used, evaluated at the converged
#' parameters, so the posterior covariances can be formed.
#'
#' @keywords internal
#' @noRd
.block_information <- function(params, Y_list, X_row_list, X_col_list,
                               X_dyad_list) {
  Tn <- length(Y_list)
  N <- nrow(params$a)
  M <- nrow(params$b)
  K <- ncol(params$U[[1]])
  P <- nrow(params$beta)

  n_a <- matrix(0, N, Tn)
  n_b <- matrix(0, M, Tn)
  G_U <- vector("list", Tn)
  G_V <- vector("list", Tn)
  S_beta <- vector("list", Tn)

  for (t in seq_len(Tn)) {
    X_row  <- if (is.null(X_row_list))  NULL else X_row_list[[t]]
    X_col  <- if (is.null(X_col_list))  NULL else X_col_list[[t]]
    X_dyad <- if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]]
    beta_t <- if (P > 0L) params$beta[, t] else numeric(0)

    cov_t <- .cov_term(beta_t, X_row, X_col, X_dyad, N, M)

    # Observation mask: a cell counts only if the outcome and every covariate
    # entering it are present.
    obs <- !is.na(Y_list[[t]] - cov_t)

    n_a[, t] <- rowSums(obs)
    n_b[, t] <- colSums(obs)

    R_mask <- ifelse(obs, 0, NA_real_)
    G_U[[t]] <- .latent_gram_rhs(R_mask, params$V[[t]])$G
    G_V[[t]] <- .latent_gram_rhs(t(R_mask), params$U[[t]])$G

    if (P > 0L) {
      Z <- .cov_design(X_row, X_col, X_dyad, N, M)
      keep <- as.vector(obs) & stats::complete.cases(Z)
      S_beta[[t]] <- crossprod(Z[keep, , drop = FALSE])
    }
  }

  list(n_a = n_a, n_b = n_b, G_U = G_U, G_V = G_V, S_beta = S_beta,
       N = N, M = M, K = K, P = P, Tn = Tn)
}


#' Empirical Bayes update of the variance components (Step 9)
#'
#' Runs the EM M-step for the eleven variance components and converts the ten
#' non-observation ones into penalties via `lambda = sigma_eps^2 / sigma^2` and
#' `gamma = sigma_eps^2 / tau^2`.
#'
#' `sigma_eps^2` is the mean squared residual over observed cells. The other ten
#' add the posterior-variance trace terms described in the file header; without
#' them the procedure drives the innovation variances to zero.
#'
#' Undefined components are reported as `NA` in `Omega`: `sigma_beta^2` and
#' `tau_beta^2` when there are no covariates, and every `tau^2` when `T = 1`.
#' Their penalties are set to `1`, the spec's neutral value, and provably cannot
#' influence the fit.
#'
#' Each variance is floored at `eps_var`, and the floored value is what `Omega`
#' reports. Penalties are additionally capped at `max_penalty` for conditioning;
#' the cap never changes the reported variances.
#'
#' @param params Parameter list, already identified and scale-normalised.
#' @param Y_list List of `T` outcome matrices.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param lambda,gamma The penalties in force during the inner loop that
#'   produced `params`; needed to rebuild the posterior covariances.
#' @param eps_var Lower bound applied to every variance component.
#' @param max_penalty Upper bound applied to each derived penalty.
#'
#' @return A list with `Omega`, `lambda`, `gamma`, `sigma_eps2`, `n_obs`,
#'   `capped`, and `trace_share` — the fraction of each variance estimate
#'   contributed by the posterior-variance correction, a useful diagnostic.
#' @keywords internal
#' @noRd
.eb_update_variances <- function(params,
                                 Y_list,
                                 X_row_list = NULL,
                                 X_col_list = NULL,
                                 X_dyad_list = NULL,
                                 lambda,
                                 gamma,
                                 eps_var = 1e-8,
                                 max_penalty = 1e6,
                                 tie_latent = TRUE) {
  Tn <- length(Y_list)
  blocks <- c("beta", "a", "b", "U", "V")
  counts <- .block_scalar_counts(params)
  lambda <- .check_penalty(lambda, "lambda", blocks)
  gamma  <- .check_penalty(gamma,  "gamma",  blocks)

  # ── Observation variance (plug-in) ───────────────────────────────────────
  ssr <- 0
  n_obs <- 0L
  for (t in seq_len(Tn)) {
    r <- .resid_period(
      Y_t    = Y_list[[t]],
      a_t    = params$a[, t], b_t = params$b[, t],
      U_t    = params$U[[t]], V_t = params$V[[t]],
      beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
      X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
      X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
      X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]]
    )
    ssr <- ssr + sum(r^2, na.rm = TRUE)
    n_obs <- n_obs + sum(!is.na(r))
  }
  if (n_obs == 0L) {
    stop("No observed cells: cannot estimate the observation variance.",
         call. = FALSE)
  }
  sigma_eps2 <- max(ssr / n_obs, eps_var)

  # ── Sums of squares of the point estimates ───────────────────────────────
  traj <- list(
    beta = .cols_to_list(params$beta),
    a    = .cols_to_list(params$a),
    b    = .cols_to_list(params$b),
    U    = params$U,
    V    = params$V
  )
  ss_init  <- stats::setNames(numeric(length(blocks)), blocks)
  ss_innov <- stats::setNames(numeric(length(blocks)), blocks)
  for (nm in blocks) {
    th <- traj[[nm]]
    ss_init[nm] <- sum(th[[1]]^2)
    if (Tn > 1L) {
      d <- 0
      for (t in 2:Tn) d <- d + sum((th[[t]] - th[[t - 1]])^2)
      ss_innov[nm] <- d
    }
  }

  # ── Posterior-variance trace corrections ─────────────────────────────────
  info <- .block_information(params, Y_list, X_row_list, X_col_list, X_dyad_list)
  tr_init  <- stats::setNames(numeric(length(blocks)), blocks)
  tr_innov <- stats::setNames(numeric(length(blocks)), blocks)

  # a and b: one scalar system per node, done in one vectorised pass each.
  for (nm in c("a", "b")) {
    n_obs_mat <- if (nm == "a") info$n_a else info$n_b
    D <- .additive_diagonal(n_obs_mat, lambda[[nm]], gamma[[nm]])
    sel <- .scalar_tridiagonal_selected_inverse_vec(D, off = -gamma[[nm]])
    tr_init[nm] <- sigma_eps2 * sum(sel$diag[, 1])
    if (Tn > 1L) {
      tot <- 0
      for (t in 2:Tn) {
        tot <- tot + sum(sel$diag[, t]) + sum(sel$diag[, t - 1]) -
          2 * sum(sel$offdiag[, t - 1])
      }
      tr_innov[nm] <- sigma_eps2 * tot
    }
  }

  # U and V: one K x K system per node.
  for (nm in c("U", "V")) {
    G_all <- if (nm == "U") info$G_U else info$G_V
    n_sys <- if (nm == "U") info$N else info$M
    K <- info$K
    Ik <- diag(K)
    sel_list <- vector("list", n_sys)
    for (i in seq_len(n_sys)) {
      D <- vector("list", Tn)
      for (t in seq_len(Tn)) {
        nb <- (t > 1L) + (t < Tn)
        D[[t]] <- matrix(G_all[[t]][i, , ], K, K) +
          (gamma[[nm]] * nb + if (t == 1L) lambda[[nm]] else 0) * Ik
      }
      sel_list[[i]] <- .block_tridiagonal_selected_inverse(D, off = -gamma[[nm]])
    }
    tt <- .trace_totals(sel_list)
    tr_init[nm]  <- sigma_eps2 * tt$init
    tr_innov[nm] <- sigma_eps2 * tt$innov
  }

  # beta: a single P x P system.
  if (info$P > 0L) {
    Ip <- diag(info$P)
    D <- vector("list", Tn)
    for (t in seq_len(Tn)) {
      nb <- (t > 1L) + (t < Tn)
      D[[t]] <- info$S_beta[[t]] +
        (gamma[["beta"]] * nb + if (t == 1L) lambda[["beta"]] else 0) * Ip
    }
    tt <- .trace_totals(list(
      .block_tridiagonal_selected_inverse(D, off = -gamma[["beta"]])
    ))
    tr_init[["beta"]]  <- sigma_eps2 * tt$init
    tr_innov[["beta"]] <- sigma_eps2 * tt$innov
  }

  # ── M-step ───────────────────────────────────────────────────────────────
  sigma2 <- stats::setNames(rep(NA_real_, length(blocks)), blocks)
  tau2   <- stats::setNames(rep(NA_real_, length(blocks)), blocks)
  share_init  <- stats::setNames(rep(NA_real_, length(blocks)), blocks)
  share_innov <- stats::setNames(rep(NA_real_, length(blocks)), blocks)

  for (nm in blocks) {
    n_scalar <- counts[[nm]]
    if (n_scalar == 0L) next

    num_init <- ss_init[[nm]] + tr_init[[nm]]
    sigma2[nm] <- max(num_init / n_scalar, eps_var)
    share_init[nm] <- if (num_init > 0) tr_init[[nm]] / num_init else NA_real_

    if (Tn > 1L) {
      num_innov <- ss_innov[[nm]] + tr_innov[[nm]]
      tau2[nm] <- max(num_innov / (n_scalar * (Tn - 1L)), eps_var)
      share_innov[nm] <- if (num_innov > 0) tr_innov[[nm]] / num_innov else NA_real_
    }
  }

  # ── One innovation variance for the latent block ─────────────────────────
  #
  # `U` and `V` can trade temporal movement between them while leaving every
  # `U_t V_t'` untouched, so the likelihood says nothing about how it is split;
  # only the product is determined. Estimated separately, the two are not
  # measured by this update but driven apart by it: each penalty is derived
  # from the variance just measured, `gamma_U = sigma_eps^2 / tau_U^2`, so a
  # smaller `tau_U^2` smooths `U` further and lowers `tau_U^2` again. On panels
  # generated with a true ratio of four the fitted ratio passed 68 and was
  # still rising when the loop stopped — it stopped because the FITTED VALUES
  # had settled, the drift running along a direction the convergence criterion
  # cannot see — and on others it fell below 1e-3. The product was recovered to
  # within twenty percent throughout.
  #
  # The geometric mean keeps what is determined and fixes the undetermined
  # ratio at one. Note what this is not: it is not a rescaling of the
  # trajectories. Imposing the same equality by rescaling `U` and `V` until
  # their increments match was tried and fails — the factor does not return to
  # one between passes, so the correction compounds, and after forty passes the
  # trajectories carry a factor of 1e7 and the solve goes singular. Leaving the
  # trajectories alone and constraining only the reported variance has no such
  # mechanism.
  #
  # `tie_latent = FALSE` restores the separate estimates. It is not a supported
  # way to fit; it exists so the drift can be measured.
  latent_ratio <- NA_real_
  if (Tn > 1L && !is.na(tau2[["U"]]) && !is.na(tau2[["V"]])) {
    latent_ratio <- tau2[["U"]] / tau2[["V"]]
    if (tie_latent) {
      kappa <- sqrt(tau2[["U"]] * tau2[["V"]])
      tau2["U"] <- kappa
      tau2["V"] <- kappa
    }
  }

  # ── Penalties ────────────────────────────────────────────────────────────
  lam_new <- stats::setNames(numeric(length(blocks)), blocks)
  gam_new <- stats::setNames(numeric(length(blocks)), blocks)
  for (nm in blocks) {
    lam_new[nm] <- if (is.na(sigma2[nm])) 1 else sigma_eps2 / sigma2[nm]
    gam_new[nm] <- if (is.na(tau2[nm]))   1 else sigma_eps2 / tau2[nm]
  }

  capped <- c(lambda = lam_new > max_penalty, gamma = gam_new > max_penalty)
  lam_new <- pmin(lam_new, max_penalty)
  gam_new <- pmin(gam_new, max_penalty)

  Omega <- c(
    sigma_eps2 = sigma_eps2,
    sigma_beta2 = sigma2[["beta"]], sigma_a2 = sigma2[["a"]],
    sigma_b2 = sigma2[["b"]], sigma_U2 = sigma2[["U"]], sigma_V2 = sigma2[["V"]],
    tau_beta2 = tau2[["beta"]], tau_a2 = tau2[["a"]],
    tau_b2 = tau2[["b"]], tau_U2 = tau2[["U"]], tau_V2 = tau2[["V"]]
  )

  list(
    Omega = Omega,
    lambda = lam_new,
    gamma = gam_new,
    sigma_eps2 = sigma_eps2,
    n_obs = n_obs,
    capped = capped,
    # What the separate estimates would have said, recorded whether or not they
    # were used. Kept because the direction the data push it is worth seeing
    # even when it is not acted on, and because a run of these across passes is
    # how the instability was found in the first place.
    latent_ratio = latent_ratio,
    trace_share = list(initial = share_init, innovation = share_innov)
  )
}
