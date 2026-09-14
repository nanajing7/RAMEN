# Conditional parametric bootstrap for fit_dynamic_ame().
#
# The fitted systematic component is held fixed and only the observation errors
# are regenerated, so every replicate carries exactly the dependence structure
# the model estimated — sender and receiver heterogeneity through the additive
# and latent effects, temporal persistence through the random-walk priors. The
# dependence is preserved by construction, not by the resampling scheme, which
# is why the errors can be drawn independently across cells.
#
# Nodes and periods are conditioned on rather than resampled: the inferential
# target is the precision of the estimated trajectories for the observed
# network, not extrapolation to other node sets or time windows.

#' Running accumulator for a stack of matrices
#'
#' Bootstrap draws of `U_t V_t'` are one `N x M` matrix per period per replicate
#' and cannot be held in memory for large panels. This accumulates the mean, the
#' sum of squared deviations (Welford), and the count of positive draws in a
#' single pass, so summaries are available without storing anything.
#'
#' @param dims Integer vector `c(N, M)`.
#' @param Tn Number of periods.
#' @return A list of accumulator state.
#' @keywords internal
#' @noRd
.new_matrix_accumulator <- function(dims, Tn) {
  blank <- function() lapply(seq_len(Tn), function(t) matrix(0, dims[1], dims[2]))
  list(n = 0L, mean = blank(), m2 = blank(), n_pos = blank(), Tn = Tn)
}

#' Add one draw to a matrix accumulator
#'
#' @param acc Accumulator from `.new_matrix_accumulator()`.
#' @param mats List of `T` matrices for this replicate.
#' @return The updated accumulator.
#' @keywords internal
#' @noRd
.accumulate_matrices <- function(acc, mats) {
  acc$n <- acc$n + 1L
  for (t in seq_len(acc$Tn)) {
    delta <- mats[[t]] - acc$mean[[t]]
    acc$mean[[t]] <- acc$mean[[t]] + delta / acc$n
    acc$m2[[t]] <- acc$m2[[t]] + delta * (mats[[t]] - acc$mean[[t]])
    acc$n_pos[[t]] <- acc$n_pos[[t]] + (mats[[t]] > 0)
  }
  acc
}

#' Finalise a matrix accumulator into mean, sd, and sign probability
#'
#' @param acc Accumulator.
#' @param dn Dimnames to attach.
#' @param periods Names for the period list.
#' @return A list with `mean`, `sd`, and `sign_prob`, each a list of `T` matrices.
#' @keywords internal
#' @noRd
.finalize_matrix_accumulator <- function(acc, dn, periods) {
  nm <- function(x) { dimnames(x) <- dn; x }
  out <- list(
    mean = lapply(acc$mean, nm),
    sd = lapply(acc$m2, function(M) nm(sqrt(M / max(acc$n - 1L, 1L)))),
    sign_prob = lapply(acc$n_pos, function(M) nm(M / acc$n))
  )
  for (k in names(out)) names(out[[k]]) <- periods
  out
}


#' Simulate one conditional bootstrap panel
#'
#' Adds fresh Gaussian observation errors to the fitted values, leaving cells
#' that were unobserved in the data unobserved: each replicate therefore carries
#' exactly the same information as the original panel.
#'
#' @param fitted_list List of `T` fitted-value matrices.
#' @param Y_list Original outcome matrices, used only for their missing pattern.
#' @param sigma_eps2 Observation variance.
#' @return A list of `T` simulated outcome matrices.
#' @keywords internal
#' @noRd
.simulate_conditional_panel <- function(fitted_list, Y_list, sigma_eps2) {
  sd_eps <- sqrt(sigma_eps2)
  lapply(seq_along(fitted_list), function(t) {
    Y <- fitted_list[[t]] +
      matrix(stats::rnorm(length(fitted_list[[t]]), sd = sd_eps),
             nrow(fitted_list[[t]]), ncol(fitted_list[[t]]))
    Y[is.na(Y_list[[t]])] <- NA_real_
    dimnames(Y) <- dimnames(Y_list[[t]])
    Y
  })
}


#' Dispersion the fitted trajectories are missing
#'
#' A fitted trajectory is a posterior mean, and a posterior mean moves less than
#' the parameter it estimates. The empirical-Bayes update recovers `tau^2` only
#' by adding a posterior-variance trace to the squared increments of the fitted
#' path: whatever that trace supplied is dispersion the fitted path does not
#' have. Generating bootstrap panels from the fitted path alone therefore poses
#' an easier problem than the data did, because the estimator is asked to track
#' a trajectory that is already as smooth as it wants trajectories to be. The
#' replicates then measure only how the noise propagates, not the estimator's
#' failure to follow a truth that moves, and the standard errors come out too
#' small.
#'
#' This returns, per scalar, the initial-state and innovation variances that
#' restore the dispersion the model itself attributes to the truth. Both follow
#' from `Omega` and the fitted trajectories, so nothing has to be recomputed:
#' `tau^2` is what the model says the truth moves by, and the squared increments
#' of the fitted path are what it actually moves by.
#'
#' Only `beta`, `a` and `b` are covered. `U` and `V` are left out deliberately:
#' `U_t V_t'` is unchanged by `(U_t A, V_t A^-T)` for any constant invertible
#' `A`, and the identification step fixes only the orthogonal and scalar parts of
#' that freedom. `tau_U^2` and `tau_V^2` can therefore be traded between the two
#' factors without changing the fit, so they describe the chosen representative
#' rather than the data, and injecting them as innovation noise would create
#' movement the data does not support. The intervals are correspondingly still
#' conservative in the part of the estimation error that comes from smoothing the
#' latent factors.
#'
#' @param fit A `dynamic_ame` object.
#' @return A named list over `beta`, `a`, `b`, each with `init` and `innov`
#'   variances, floored at zero.
#' @keywords internal
#' @noRd
.dispersion_scales <- function(fit) {
  Om <- fit$variance_components$Omega
  Tn <- length(fit$years)

  traj <- list(beta = fit$beta, a = fit$a, b = fit$b)
  keys <- list(beta = c("sigma_beta2", "tau_beta2"),
               a    = c("sigma_a2",    "tau_a2"),
               b    = c("sigma_b2",    "tau_b2"))

  lapply(stats::setNames(names(traj), names(traj)), function(nm) {
    th <- traj[[nm]]
    n_scalar <- nrow(th)
    if (is.null(n_scalar) || n_scalar == 0L) return(list(init = 0, innov = 0))

    sigma2 <- Om[[keys[[nm]][1]]]
    ss_init <- sum(th[, 1]^2) / n_scalar
    v_init <- if (is.na(sigma2)) 0 else max(sigma2 - ss_init, 0)

    v_innov <- 0
    if (Tn > 1L) {
      tau2 <- Om[[keys[[nm]][2]]]
      d <- th[, -1, drop = FALSE] - th[, -Tn, drop = FALSE]
      ss_innov <- sum(d^2) / (n_scalar * (Tn - 1L))
      v_innov <- if (is.na(tau2)) 0 else max(tau2 - ss_innov, 0)
    }

    list(init = v_init, innov = v_innov)
  })
}


#' Draw one set of dispersed trajectories
#'
#' Keeps the increments of the fitted path and adds the missing dispersion on
#' top, so the drawn trajectory has initial-state variance `sigma^2` and
#' innovation variance `tau^2` — the model's own account of how the truth
#' behaves — while still passing through the neighbourhood of the estimate.
#'
#' @param fit A `dynamic_ame` object.
#' @param scales Output of `.dispersion_scales()`.
#' @return A list with `beta`, `a` and `b`, each the same shape as in `fit`.
#' @keywords internal
#' @noRd
.draw_dispersed_trajectories <- function(fit, scales) {
  Tn <- length(fit$years)

  lapply(stats::setNames(c("beta", "a", "b"), c("beta", "a", "b")), function(nm) {
    th <- fit[[nm]]
    n_scalar <- nrow(th)
    if (is.null(n_scalar) || n_scalar == 0L) return(th)

    sc <- scales[[nm]]
    new <- th
    new[, 1] <- th[, 1] + stats::rnorm(n_scalar, sd = sqrt(sc$init))
    if (Tn > 1L) {
      for (t in 2:Tn) {
        new[, t] <- new[, t - 1L] + (th[, t] - th[, t - 1L]) +
          stats::rnorm(n_scalar, sd = sqrt(sc$innov))
      }
    }
    new
  })
}


#' Rebuild the fitted values from perturbed additive and coefficient paths
#'
#' `U` and `V` are taken from the fit unchanged; see `.dispersion_scales()`.
#'
#' @param fit A `dynamic_ame` object.
#' @param beta,a,b Trajectories to use in place of the fitted ones.
#' @return A list of `T` fitted-value matrices.
#' @keywords internal
#' @noRd
.fitted_from_params <- function(fit, beta, a, b) {
  Tn <- length(fit$years)
  cv <- fit$covariates
  P <- nrow(beta)

  out <- lapply(seq_len(Tn), function(t) {
    f <- .fitted_period(
      a_t = a[, t], b_t = b[, t],
      U_t = fit$U[[t]], V_t = fit$V[[t]],
      beta_t = if (!is.null(P) && P > 0L) beta[, t] else numeric(0),
      X_row  = if (is.null(cv$X_row_list))  NULL else cv$X_row_list[[t]],
      X_col  = if (is.null(cv$X_col_list))  NULL else cv$X_col_list[[t]],
      X_dyad = if (is.null(cv$X_dyad_list)) NULL else cv$X_dyad_list[[t]]
    )
    dimnames(f) <- list(fit$row_ids, fit$col_ids)
    f
  })
  # named by period, exactly as fitted() returns them, so the two are
  # interchangeable wherever a panel of fitted values is expected
  names(out) <- as.character(fit$years)
  out
}


#' Re-estimate the model on a simulated panel
#'
#' Runs the same estimation procedure as the original fit, at the same settings.
#' With `warm_start = TRUE` the inner loop begins from the original estimate
#' rather than from a fresh multistart: the simulated panel is a small
#' perturbation of the observed one, so this keeps replicates in the same basin
#' of attraction and makes the bootstrap measure sampling variability rather
#' than variability of the optimiser.
#'
#' @param Y_list Simulated outcome matrices.
#' @param fit The original `dynamic_ame` object.
#' @param warm_start Start from the original estimate instead of a multistart.
#' @return A list with `params`, `Omega`, `objective`, and `converged`, or `NULL`
#'   if estimation failed.
#' @keywords internal
#' @noRd
.refit_bootstrap <- function(Y_list, fit, warm_start = TRUE) {
  st <- fit$settings
  cv <- fit$covariates
  one <- stats::setNames(rep(1, 5), c("beta", "a", "b", "U", "V"))

  res <- try(suppressWarnings({
    if (warm_start) {
      start <- list(a = fit$a, b = fit$b, beta = fit$beta, U = fit$U, V = fit$V)
      lam <- fit$variance_components$lambda
      gam <- fit$variance_components$gamma
    } else {
      ini <- .init_ame_trajectory(
        Y_list, K = st$K, P = cv$dims$P,
        n_starts = st$n_starts, perturb_sd = st$perturb_sd, seed = st$seed
      )
      cand <- lapply(ini$candidates, function(cd)
        .ame_inner_bcd(cd, Y_list, cv$X_row_list, cv$X_col_list, cv$X_dyad_list,
                       one, one, max_iter = st$inner_max_iter,
                       eps_Q = st$eps_Q, delta = st$delta))
      start <- cand[[which.min(vapply(cand, function(x) x$objective$Q,
                                      numeric(1)))]]$params
      lam <- one
      gam <- one
    }

    .ame_outer_eb(
      params = start, reference = fit$reference, Y_list = Y_list,
      X_row_list = cv$X_row_list, X_col_list = cv$X_col_list,
      X_dyad_list = cv$X_dyad_list, lambda = lam, gamma = gam,
      outer_max_iter = st$outer_max_iter, inner_max_iter = st$inner_max_iter,
      eps_Q = st$eps_Q, eps_Omega = st$eps_Omega, eps_fit = st$eps_fit,
      eps_var = st$eps_var, max_penalty = st$max_penalty, delta = st$delta,
      # Every replicate is estimated the way the original fit was, acceleration
      # included: a bootstrap run at settings the point estimate was not run at
      # measures the wrong thing. Fits made before this argument existed carry
      # no such setting, and get the plain iteration they were made with.
      accelerate = if (is.null(st$accelerate)) "none" else st$accelerate,
      # Fits made before this argument existed were estimated under the pooled
      # convention, and a replicate re-estimated under a different one would
      # not be a replicate of that fit.
      gauge_anchor = if (is.null(st$gauge)) "pooled" else st$gauge,
      tie_latent = if (is.null(st$latent_innovation)) FALSE else
        st$latent_innovation == "shared"
    )
  }), silent = TRUE)

  if (inherits(res, "try-error")) return(NULL)

  list(params = res$params, Omega = res$Omega,
       objective = utils::tail(res$objective_trace, 1),
       converged = res$converged)
}


#' Conditional parametric bootstrap (design a)
#'
#' Internal engine behind [bootstrap_ame()] with `design = "conditional"`.
#' Quantifies the precision of the estimated parameters for the observed
#' network: the fitted systematic component is held fixed and only the
#' observation errors are regenerated, and each replicate is re-estimated with
#' the full procedure, including the empirical-Bayes variance-component update.
#'
#' @param fit A `dynamic_ame` object.
#' @param B Number of bootstrap replications.
#' @param conf_level Confidence level for percentile intervals.
#' @param warm_start Initialise each replicate at the original estimate.
#' @param dispersion `"posterior"` restores the dispersion the fitted
#'   trajectories are missing before generating each panel; `"none"` generates
#'   from the fitted values themselves. See `.dispersion_scales()`.
#' @param store_latent_draws Keep every draw of `U_t V_t'`.
#' @param n_cores Number of parallel workers.
#' @param seed Base seed.
#' @param verbose Report progress.
#'
#' @return An object of class `bootstrap_ame_conditional`.
#' @keywords internal
#' @noRd
.bootstrap_conditional <- function(fit,
                                   B = 1000,
                                   conf_level = 0.95,
                                   warm_start = TRUE,
                                   dispersion = c("posterior", "none"),
                                   store_latent_draws = FALSE,
                                   n_cores = 1,
                                   seed = 1,
                                   verbose = FALSE) {
  .validate_dynamic_ame(fit)
  dispersion <- match.arg(dispersion)

  if (length(B) != 1L || B < 2 || B != as.integer(B)) {
    stop("B must be an integer of at least 2.", call. = FALSE)
  }
  if (length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be strictly between 0 and 1.", call. = FALSE)
  }
  B <- as.integer(B)

  Tn <- length(fit$years)
  N <- nrow(fit$a)
  M <- nrow(fit$b)
  P <- nrow(fit$beta)
  fitted_list <- fitted(fit)
  disp <- if (dispersion == "posterior") .dispersion_scales(fit) else NULL

  # One seed per replicate, drawn up front, so the result is identical whether
  # the replicates run sequentially or in parallel and in whatever order.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  one_rep <- function(b) {
    set.seed(rep_seeds[b])
    # With dispersion off the RNG goes straight to the errors, so the stream --
    # and the whole result -- is what it was before this option existed.
    mu <- if (dispersion == "posterior") {
      tr <- .draw_dispersed_trajectories(fit, disp)
      .fitted_from_params(fit, tr$beta, tr$a, tr$b)
    } else {
      fitted_list
    }
    Y_star <- .simulate_conditional_panel(mu, fit$Y_list, fit$sigma_eps2)
    out <- .refit_bootstrap(Y_star, fit, warm_start = warm_start)
    if (is.null(out)) return(NULL)
    list(
      beta = out$params$beta,
      a = out$params$a,
      b = out$params$b,
      M = lapply(seq_len(Tn), function(t)
        out$params$U[[t]] %*% t(out$params$V[[t]])),
      Omega = out$Omega,
      objective = out$objective,
      converged = out$converged
    )
  }

  if (verbose) {
    cat(sprintf("Conditional parametric bootstrap: B = %d, %d core(s)\n",
                B, n_cores))
  }

  reps <- if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE) &&
              .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(B), one_rep, mc.cores = n_cores)
  } else {
    if (n_cores > 1L) {
      warning("Parallel execution unavailable here; running sequentially.",
              call. = FALSE)
    }
    lapply(seq_len(B), function(b) {
      if (verbose && b %% 50L == 0L) cat(sprintf("  replicate %d / %d\n", b, B))
      one_rep(b)
    })
  }

  ok <- !vapply(reps, is.null, logical(1))
  if (!any(ok)) stop("All bootstrap replications failed.", call. = FALSE)
  reps <- reps[ok]
  B_ok <- length(reps)

  # ── Draws for the small blocks; accumulator for the interaction matrices ──
  draws <- list(
    beta = if (P > 0L) vapply(reps, function(r) as.vector(r$beta),
                              numeric(P * Tn)) else NULL,
    a = vapply(reps, function(r) as.vector(r$a), numeric(N * Tn)),
    b = vapply(reps, function(r) as.vector(r$b), numeric(M * Tn)),
    Omega = vapply(reps, function(r) r$Omega, numeric(11))
  )

  acc <- .new_matrix_accumulator(c(N, M), Tn)
  for (r in reps) acc <- .accumulate_matrices(acc, r$M)
  latent <- .finalize_matrix_accumulator(
    acc, list(fit$row_ids, fit$col_ids), as.character(fit$years))

  if (store_latent_draws) {
    latent$draws <- lapply(reps, function(r) r$M)
  } else {
    # percentile intervals need order statistics, which the accumulator cannot
    # supply; fall back to a normal approximation and say so
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    latent$conf.low <- lapply(seq_len(Tn), function(t)
      latent$mean[[t]] - z * latent$sd[[t]])
    latent$conf.high <- lapply(seq_len(Tn), function(t)
      latent$mean[[t]] + z * latent$sd[[t]])
    names(latent$conf.low) <- names(latent$conf.high) <- as.character(fit$years)
    latent$interval_type <- "normal approximation"
  }

  if (store_latent_draws) {
    probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
    qmat <- function(t, p) {
      arr <- vapply(latent$draws, function(d) as.vector(d[[t]]), numeric(N * M))
      m <- matrix(apply(arr, 1, stats::quantile, probs = p), N, M)
      dimnames(m) <- list(fit$row_ids, fit$col_ids)
      m
    }
    latent$conf.low <- lapply(seq_len(Tn), qmat, p = probs[1])
    latent$conf.high <- lapply(seq_len(Tn), qmat, p = probs[2])
    names(latent$conf.low) <- names(latent$conf.high) <- as.character(fit$years)
    latent$interval_type <- "percentile"
  }

  structure(
    list(
      design = "conditional",
      fit = fit,
      B_requested = B,
      B_successful = B_ok,
      conf_level = conf_level,
      warm_start = warm_start,
      dispersion = dispersion,
      dispersion_scales = disp,
      summaries = list(
        beta = .summarize_draws(draws$beta, fit$beta, conf_level,
                                rownames(fit$beta), fit$years, "covariate"),
        a = .summarize_draws(draws$a, fit$a, conf_level,
                             fit$row_ids, fit$years, "sender"),
        b = .summarize_draws(draws$b, fit$b, conf_level,
                             fit$col_ids, fit$years, "receiver"),
        Omega = .summarize_omega(draws$Omega,
                                 fit$variance_components$Omega, conf_level)
      ),
      latent = latent,
      diagnostics = list(
        objective = vapply(reps, function(r) r$objective, numeric(1)),
        converged = vapply(reps, function(r) isTRUE(r$converged), logical(1))
      ),
      settings = list(B = B, conf_level = conf_level, seed = seed,
                      warm_start = warm_start, dispersion = dispersion,
                      n_cores = n_cores,
                      store_latent_draws = store_latent_draws)
    ),
    class = "bootstrap_ame_conditional"
  )
}


#' Summarise bootstrap draws of a `q x T` parameter block
#'
#' @param draws A `(q*T) x B` matrix of vectorised draws, or `NULL`.
#' @param est The original `q x T` estimate.
#' @param conf_level Confidence level.
#' @param unit_names Names for the rows of the block.
#' @param periods Period labels.
#' @param unit_label Column name for the unit identifier.
#'
#' @return A long data frame, or `NULL` when the block is empty.
#' @keywords internal
#' @noRd
.summarize_draws <- function(draws, est, conf_level, unit_names, periods,
                             unit_label) {
  if (is.null(draws) || nrow(est) == 0L) return(NULL)

  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  q <- nrow(est)
  Tn <- ncol(est)

  se <- apply(draws, 1, stats::sd)
  lo <- apply(draws, 1, stats::quantile, probs = probs[1])
  hi <- apply(draws, 1, stats::quantile, probs = probs[2])

  out <- data.frame(
    unit = rep(unit_names %||% seq_len(q), times = Tn),
    period = rep(periods, each = q),
    estimate = as.vector(est),
    std.error = as.numeric(se),
    conf.low = as.numeric(lo),
    conf.high = as.numeric(hi),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- unit_label
  out
}


#' Summarise bootstrap draws of the variance components
#'
#' @param draws An `11 x B` matrix of draws.
#' @param est The original eleven-entry `Omega`.
#' @param conf_level Confidence level.
#' @return A data frame with one row per component.
#' @keywords internal
#' @noRd
.summarize_omega <- function(draws, est, conf_level) {
  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  keep <- !is.na(est)

  data.frame(
    component = names(est),
    estimate = as.numeric(est),
    std.error = ifelse(keep, apply(draws, 1, stats::sd, na.rm = TRUE), NA_real_),
    conf.low = ifelse(keep, apply(draws, 1, stats::quantile,
                                  probs = probs[1], na.rm = TRUE), NA_real_),
    conf.high = ifelse(keep, apply(draws, 1, stats::quantile,
                                   probs = probs[2], na.rm = TRUE), NA_real_),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


#' Print a conditional bootstrap result
#'
#' @param x A `bootstrap_ame_conditional` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.bootstrap_ame_conditional <- function(x, ...) {
  cat("Conditional parametric bootstrap for a dynamic bipartite AME model\n")
  cat(sprintf("Replications: %d of %d succeeded (%.1f%%)\n",
              x$B_successful, x$B_requested,
              100 * x$B_successful / x$B_requested))
  cat(sprintf("Confidence level: %.0f%%  |  initialisation: %s\n",
              100 * x$conf_level,
              if (isTRUE(x$warm_start)) "warm start at the original estimate"
              else "fresh multistart"))
  cat(sprintf("Generating trajectories: %s\n",
              if (identical(x$dispersion, "none"))
                "fitted values as they stand"
              else paste("fitted values with the missing posterior dispersion",
                         "restored for beta, a and b")))

  conv <- x$diagnostics$converged
  cat(sprintf("Replicates reaching outer convergence: %d of %d\n",
              sum(conv), length(conv)))

  obj <- x$diagnostics$objective
  cat(sprintf("Converged objective across replicates: min %.4f, median %.4f, max %.4f\n",
              min(obj), stats::median(obj), max(obj)))
  cat("  (a single tight cluster indicates the replicates stayed in one basin)\n")

  if (!is.null(x$summaries$beta)) {
    cat("\nCoefficient trajectories:\n")
    print(utils::head(x$summaries$beta, 10), row.names = FALSE)
    if (nrow(x$summaries$beta) > 10) cat("  ...\n")
  }

  cat("\nVariance components:\n")
  print(x$summaries$Omega, row.names = FALSE)

  cat(sprintf("\nInteraction matrices summarised over %d periods (%s intervals)\n",
              length(x$latent$mean), x$latent$interval_type))

  invisible(x)
}
