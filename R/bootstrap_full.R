# Full-model parametric bootstrap for fit_dynamic_ame().
#
# The estimated model is treated as a data-generating process in full: every
# state trajectory is redrawn from the priors that the fitted variance
# components parameterise, and a complete panel is generated from those states.
# Each replicate is therefore a fresh realization of the model rather than a
# perturbation of the observed network.
#
# Because the generating parameters differ from replicate to replicate, the
# dispersion of the estimates across replicates does not measure estimator
# performance: it is governed by how widely the generating values themselves
# vary, and an estimator that recovered every replicate exactly would show the
# same dispersion. What this design provides instead is knowledge of the values
# behind each replicate, so each estimate can be compared with the value it is
# trying to recover. Everything reported here is such a comparison.

#' Draw one trajectory from the model's prior
#'
#' The initial state is drawn from `N(0, sigma2)` — the model's prior is centred
#' at zero — and propagated forward by a random walk with innovation variance
#' `tau2`. With a single period the random walk plays no part.
#'
#' @param n_row,n_col Trajectory dimensions (`n_col` is the number of periods).
#' @param sigma2 Initial-state variance.
#' @param tau2 Innovation variance; ignored when `n_col == 1`.
#' @return An `n_row x n_col` matrix.
#' @keywords internal
#' @noRd
.draw_rw_trajectory <- function(n_row, n_col, sigma2, tau2) {
  out <- matrix(0, n_row, n_col)
  if (n_row == 0L) return(out)

  out[, 1] <- stats::rnorm(n_row, sd = sqrt(sigma2))
  if (n_col > 1L) {
    for (t in 2:n_col) {
      out[, t] <- out[, t - 1] + stats::rnorm(n_row, sd = sqrt(tau2))
    }
  }
  out
}


#' Simulate one complete panel from the fitted model
#'
#' Redraws every state trajectory from the priors implied by the fitted variance
#' components, then generates the outcome. Covariates are held at their observed
#' values — the model specifies the conditional distribution of the outcome given
#' the covariates and does not model the covariates themselves. Dyad-periods
#' that are unobserved in the data remain unobserved.
#'
#' @param fit A `dynamic_ame` object.
#' @return A list with the simulated `Y_list` and the generating `params`.
#' @keywords internal
#' @noRd
.simulate_full_model_panel <- function(fit) {
  Om <- fit$variance_components$Omega
  Tn <- length(fit$years)
  N <- nrow(fit$a)
  M <- nrow(fit$b)
  P <- nrow(fit$beta)
  K <- ncol(fit$U[[1]])
  cv <- fit$covariates

  # A NA variance means the component is undefined (no covariates, or T = 1);
  # the corresponding trajectory is degenerate and contributes nothing.
  sv <- function(x) if (is.na(Om[[x]])) 0 else Om[[x]]

  gen <- list(
    beta = .draw_rw_trajectory(P, Tn, sv("sigma_beta2"), sv("tau_beta2")),
    a = .draw_rw_trajectory(N, Tn, sv("sigma_a2"), sv("tau_a2")),
    b = .draw_rw_trajectory(M, Tn, sv("sigma_b2"), sv("tau_b2"))
  )

  # U and V are drawn one latent dimension at a time, each dimension an
  # independent random walk with the same variances.
  draw_latent <- function(n, sigma2, tau2) {
    mats <- lapply(seq_len(Tn), function(t) matrix(0, n, K))
    traj <- .draw_rw_trajectory(n * K, Tn, sigma2, tau2)
    for (t in seq_len(Tn)) mats[[t]] <- matrix(traj[, t], n, K)
    mats
  }
  gen$U <- draw_latent(N, sv("sigma_U2"), sv("tau_U2"))
  gen$V <- draw_latent(M, sv("sigma_V2"), sv("tau_V2"))

  sd_eps <- sqrt(Om[["sigma_eps2"]])
  Y_list <- lapply(seq_len(Tn), function(t) {
    Y <- .fitted_period(
      a_t = gen$a[, t], b_t = gen$b[, t],
      U_t = gen$U[[t]], V_t = gen$V[[t]],
      beta_t = if (P > 0L) gen$beta[, t] else numeric(0),
      X_row  = if (is.null(cv$X_row_list))  NULL else cv$X_row_list[[t]],
      X_col  = if (is.null(cv$X_col_list))  NULL else cv$X_col_list[[t]],
      X_dyad = if (is.null(cv$X_dyad_list)) NULL else cv$X_dyad_list[[t]]
    ) + matrix(stats::rnorm(N * M, sd = sd_eps), N, M)
    Y[is.na(fit$Y_list[[t]])] <- NA_real_
    dimnames(Y) <- dimnames(fit$Y_list[[t]])
    Y
  })

  list(Y_list = Y_list, params = gen)
}


#' Gauge-invariant summary of the latent variance components
#'
#' `sigma_U2`, `sigma_V2`, `tau_U2`, and `tau_V2` are not comparable across a
#' rescaling `U -> cU`, `V -> V/c`: the first two scale by `c^2` and the last two
#' by `c^-2`. The generated trajectories are drawn from the priors and do not
#' satisfy the estimator's scale convention, whereas the estimates do, so the
#' individual components would differ by an arbitrary gauge factor. Their
#' products are invariant and are reported instead.
#'
#' @param Omega An eleven-entry variance-component vector.
#' @return A named vector of the seven gauge-free components and the two
#'   invariant products.
#' @keywords internal
#' @noRd
.gauge_invariant_omega <- function(Omega) {
  free <- c("sigma_eps2", "sigma_beta2", "sigma_a2", "sigma_b2",
            "tau_beta2", "tau_a2", "tau_b2")
  c(Omega[free],
    sigma_UV2_product = unname(Omega[["sigma_U2"]] * Omega[["sigma_V2"]]),
    tau_UV2_product = unname(Omega[["tau_U2"]] * Omega[["tau_V2"]]))
}


#' Re-estimate a simulated panel from scratch
#'
#' Unlike the conditional design, no warm start is possible: the generating
#' parameters differ from the original estimate, so initialising there would
#' start the optimiser in the wrong place. Each replicate therefore runs the
#' complete procedure, multistart included.
#'
#' @param Y_list Simulated outcome matrices.
#' @param fit The original `dynamic_ame` object, used for settings and covariates.
#' @return A list with `params`, `Omega`, `objective`, `converged`, or `NULL`.
#' @keywords internal
#' @noRd
.refit_full_model <- function(Y_list, fit) {
  st <- fit$settings
  cv <- fit$covariates
  one <- stats::setNames(rep(1, 5), c("beta", "a", "b", "U", "V"))

  res <- try(suppressWarnings({
    ini <- .init_ame_trajectory(
      Y_list, K = st$K, P = cv$dims$P,
      n_starts = st$n_starts, perturb_sd = st$perturb_sd, seed = st$seed
    )
    cand <- lapply(ini$candidates, function(cd)
      .ame_inner_bcd(cd, Y_list, cv$X_row_list, cv$X_col_list, cv$X_dyad_list,
                     one, one, max_iter = st$inner_max_iter,
                     eps_Q = st$eps_Q, delta = st$delta))
    best <- which.min(vapply(cand, function(x) x$objective$Q, numeric(1)))

    .ame_outer_eb(
      params = cand[[best]]$params, reference = ini$reference,
      Y_list = Y_list, X_row_list = cv$X_row_list,
      X_col_list = cv$X_col_list, X_dyad_list = cv$X_dyad_list,
      lambda = one, gamma = one,
      outer_max_iter = st$outer_max_iter, inner_max_iter = st$inner_max_iter,
      eps_Q = st$eps_Q, eps_Omega = st$eps_Omega, eps_fit = st$eps_fit,
      eps_var = st$eps_var, max_penalty = st$max_penalty, delta = st$delta
    )
  }), silent = TRUE)

  if (inherits(res, "try-error")) return(NULL)

  list(params = res$params, Omega = res$Omega,
       objective = utils::tail(res$objective_trace, 1),
       converged = res$converged)
}


#' Full-model parametric bootstrap (design b)
#'
#' Internal engine behind [bootstrap_ame()] with `design = "full"`. Evaluates
#' how well the estimation procedure recovers parameters when applied to fresh
#' realizations of the fitted model: every state trajectory is redrawn from the
#' priors the fitted variance components parameterise, a complete panel is
#' generated, and the model is re-estimated from scratch.
#'
#' @param fit A `dynamic_ame` object.
#' @param B Number of replications.
#' @param n_cores Number of parallel workers.
#' @param seed Base seed.
#' @param verbose Report progress.
#'
#' @return An object of class `bootstrap_ame_full`.
#' @keywords internal
#' @noRd
.bootstrap_full <- function(fit,
                            B = 1000,
                            n_cores = 1,
                            seed = 1,
                            verbose = FALSE) {
  .validate_dynamic_ame(fit)
  if (length(B) != 1L || B < 2 || B != as.integer(B)) {
    stop("B must be an integer of at least 2.", call. = FALSE)
  }
  B <- as.integer(B)

  Tn <- length(fit$years)
  N <- nrow(fit$a)
  M <- nrow(fit$b)
  P <- nrow(fit$beta)

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  one_rep <- function(b) {
    set.seed(rep_seeds[b])
    sim <- .simulate_full_model_panel(fit)
    est <- .refit_full_model(sim$Y_list, fit)
    if (is.null(est)) return(NULL)

    M_true <- lapply(seq_len(Tn), function(t)
      sim$params$U[[t]] %*% t(sim$params$V[[t]]))
    M_hat <- lapply(seq_len(Tn), function(t)
      est$params$U[[t]] %*% t(est$params$V[[t]]))

    list(
      err_beta = if (P > 0L) est$params$beta - sim$params$beta else NULL,
      err_a = est$params$a - sim$params$a,
      err_b = est$params$b - sim$params$b,
      omega_hat = est$Omega,
      omega_true = fit$variance_components$Omega,
      D = vapply(seq_len(Tn), function(t)
        norm(M_hat[[t]] - M_true[[t]], "F") / norm(M_true[[t]], "F"),
        numeric(1)),
      r = vapply(seq_len(Tn), function(t)
        stats::cor(as.vector(M_hat[[t]]), as.vector(M_true[[t]])),
        numeric(1)),
      objective = est$objective,
      converged = est$converged
    )
  }

  if (verbose) {
    cat(sprintf("Full-model parametric bootstrap: B = %d, %d core(s)\n",
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
      if (verbose && b %% 25L == 0L) cat(sprintf("  replicate %d / %d\n", b, B))
      one_rep(b)
    })
  }

  ok <- !vapply(reps, is.null, logical(1))
  if (!any(ok)) stop("All bootstrap replications failed.", call. = FALSE)
  reps <- reps[ok]
  B_ok <- length(reps)

  # ── Estimation errors for the state blocks ───────────────────────────────
  err <- list(
    beta = if (P > 0L)
      vapply(reps, function(r) as.vector(r$err_beta), numeric(P * Tn)) else NULL,
    a = vapply(reps, function(r) as.vector(r$err_a), numeric(N * Tn)),
    b = vapply(reps, function(r) as.vector(r$err_b), numeric(M * Tn))
  )

  # ── Variance components, on the gauge-invariant scale ────────────────────
  om_true <- .gauge_invariant_omega(fit$variance_components$Omega)
  om_hat <- vapply(reps, function(r) .gauge_invariant_omega(r$omega_hat),
                   numeric(length(om_true)))
  om_err <- om_hat - om_true

  structure(
    list(
      design = "full-model",
      fit = fit,
      B_requested = B,
      B_successful = B_ok,
      accuracy = list(
        beta = .summarize_errors(err$beta, fit$beta, rownames(fit$beta),
                                 fit$years, "covariate"),
        a = .summarize_errors(err$a, fit$a, fit$row_ids, fit$years, "sender"),
        b = .summarize_errors(err$b, fit$b, fit$col_ids, fit$years, "receiver"),
        Omega = data.frame(
          component = names(om_true),
          generating = as.numeric(om_true),
          mean_error = ifelse(is.na(om_true), NA_real_, rowMeans(om_err)),
          rmse = ifelse(is.na(om_true), NA_real_,
                        sqrt(rowMeans(om_err^2))),
          relative_rmse = ifelse(is.na(om_true), NA_real_,
                                 sqrt(rowMeans(om_err^2)) / abs(om_true)),
          row.names = NULL, stringsAsFactors = FALSE
        )
      ),
      latent_recovery = .summarize_latent_recovery(reps, fit$years),
      draws = list(
        omega_full = vapply(reps, function(r) r$omega_hat, numeric(11)),
        omega_error = om_err
      ),
      diagnostics = list(
        objective = vapply(reps, function(r) r$objective, numeric(1)),
        converged = vapply(reps, function(r) isTRUE(r$converged), logical(1))
      ),
      settings = list(B = B, seed = seed, n_cores = n_cores)
    ),
    class = "bootstrap_ame_full"
  )
}


#' Summarise estimation errors for a parameter block
#'
#' @param err A `(q*T) x B` matrix of errors, or `NULL`.
#' @param est The original estimate, used only for its shape and labels.
#' @param unit_names Row labels.
#' @param periods Period labels.
#' @param unit_label Column name for the unit identifier.
#' @return A long data frame with the mean error and RMSE per unit-period.
#' @keywords internal
#' @noRd
.summarize_errors <- function(err, est, unit_names, periods, unit_label) {
  if (is.null(err) || nrow(est) == 0L) return(NULL)

  q <- nrow(est)
  Tn <- ncol(est)
  out <- data.frame(
    unit = rep(unit_names %||% seq_len(q), times = Tn),
    period = rep(periods, each = q),
    mean_error = rowMeans(err),
    rmse = sqrt(rowMeans(err^2)),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- unit_label
  out
}


#' Summarise latent-structure recovery across replicates
#'
#' @param reps List of per-replicate results.
#' @param periods Period labels.
#' @return A data frame with one row per period.
#' @keywords internal
#' @noRd
.summarize_latent_recovery <- function(reps, periods) {
  D <- vapply(reps, function(r) r$D, numeric(length(periods)))
  r <- vapply(reps, function(x) x$r, numeric(length(periods)))
  if (is.null(dim(D))) { D <- matrix(D, nrow = 1); r <- matrix(r, nrow = 1) }

  data.frame(
    period = periods,
    frobenius_mean = rowMeans(D),
    frobenius_sd = apply(D, 1, stats::sd),
    correlation_mean = rowMeans(r),
    correlation_sd = apply(r, 1, stats::sd),
    row.names = NULL, stringsAsFactors = FALSE
  )
}


#' Print a full-model bootstrap result
#'
#' @param x A `bootstrap_ame_full` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.bootstrap_ame_full <- function(x, ...) {
  cat("Full-model parametric bootstrap for a dynamic bipartite AME model\n")
  cat(sprintf("Replications: %d of %d succeeded (%.1f%%)\n",
              x$B_successful, x$B_requested,
              100 * x$B_successful / x$B_requested))
  cat(sprintf("Replicates reaching outer convergence: %d of %d\n",
              sum(x$diagnostics$converged), length(x$diagnostics$converged)))

  cat("\nEstimation accuracy against the generating values.\n")
  cat("(The dispersion of the estimates is not reported: it reflects the\n")
  cat(" variation of the generating values, not estimator performance.)\n")

  for (nm in c("beta", "a", "b")) {
    s <- x$accuracy[[nm]]
    if (is.null(s)) next
    cat(sprintf("\n  %-6s  mean error %+.5f   RMSE %.5f\n",
                nm, mean(s$mean_error), sqrt(mean(s$rmse^2))))
  }

  cat("\nVariance components (gauge-invariant scale):\n")
  print(x$accuracy$Omega, row.names = FALSE, digits = 4)

  cat("\nLatent-structure recovery by period:\n")
  print(x$latent_recovery, row.names = FALSE, digits = 4)
  cat("  relative Frobenius distance near zero and correlation near one\n")
  cat("  indicate accurate recovery of the interaction matrices\n")

  invisible(x)
}
