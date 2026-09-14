# Public entry point for the joint-trajectory MAP + empirical Bayes estimator.

#' Fit a dynamic bipartite AME model by joint-trajectory MAP with empirical Bayes
#'
#' Estimates the dynamic additive and multiplicative effects model
#' \deqn{Y_{ijt} = a_{it} + b_{jt} + \mathbf{x}_{ijt}'\boldsymbol{\beta}_t
#'                 + \mathbf{u}_{it}'\mathbf{v}_{jt} + \epsilon_{ijt}}
#' where every parameter follows a Gaussian random walk over time. All five
#' parameter trajectories are estimated jointly as complete trajectories, rather
#' than period by period, and the eleven variance components are estimated by
#' plug-in empirical Bayes.
#'
#' @section Algorithm:
#' A nested loop. The inner block coordinate descent updates each block —
#' `beta`, each sender effect `a_i`, each receiver effect `b_j`, each sender
#' latent trajectory `u_i`, each receiver latent trajectory `v_j` — as a whole
#' trajectory, by solving the block-tridiagonal system its first-order
#' conditions produce. Each such update is the exact minimiser of the penalised
#' objective over that block, so the objective is non-increasing across the
#' sweep. After the inner loop converges, the latent factors are identified
#' (scale gauge, then rotation to the initialisation reference) and the variance
#' components are refreshed, giving new penalties for the next pass. The outer
#' loop stops when both the variance vector and the fitted values stabilise.
#'
#' The multiplicative term makes the problem non-convex in `(U, V)` jointly, so
#' `n_starts` perturbed initialisations are run through the inner loop at neutral
#' penalties and the one reaching the lowest objective is carried forward.
#'
#' @section Missing data:
#' Cells that are unobserved, or whose covariates are incomplete, are excluded
#' from every sum. With a complete panel the computation reduces exactly to the
#' textbook formulas.
#'
#' @section What is and is not identified:
#' `U_t V_t'` is identified; the individual factors are not, beyond the gauge
#' this function fixes. Interpret rotation-invariant quantities — the fitted
#' multiplicative component, fitted values, the coefficient trajectories.
#' `sigma_U2` and `sigma_V2` read period 1 only and are shrunk by their own
#' penalty, so they are the least precise entries of `variance_components`; the
#' innovation variances are far better determined.
#'
#' @param edge_panel Data frame holding the weighted bipartite panel.
#' @param row_cov_df,col_cov_df,dyad_cov_df Optional covariate data frames.
#' @param row_covar_names,col_covar_names,dyad_covar_names Covariate names.
#' @param time_col,row_col,col_col,value_col Column names in `edge_panel`.
#' @param row_prefix,col_prefix Prefixes applied to node IDs.
#' @param row_pad,col_pad Optional left-padding widths for node IDs.
#' @param transform Function applied elementwise to the outcome matrices.
#' @param row_cov_time_col,row_cov_id_col,row_cov_prefix Mapping for `row_cov_df`.
#' @param col_cov_time_col,col_cov_id_col,col_cov_prefix Mapping for `col_cov_df`.
#' @param dyad_cov_time_col,dyad_cov_row_id_col,dyad_cov_col_id_col,dyad_cov_row_prefix,dyad_cov_col_prefix
#'   Mapping for `dyad_cov_df`.
#' @param demean_row_covariates,demean_col_covariates Subtract each node's time
#'   mean from its covariates before fitting. Recommended for time-varying node
#'   covariates, which are otherwise collinear with the node effects.
#' @param K Latent dimension.
#' @param n_starts Number of multistart candidates; the first is unperturbed.
#' @param perturb_sd Standard deviation of the perturbation applied to the
#'   latent factors in candidates `2, ..., n_starts`.
#' @param inner_max_iter,outer_max_iter Iteration caps.
#' @param eps_Q Relative tolerance on the objective, inner loop. Slack here is
#'   largely absorbed by the outer loop, which re-derives the penalties and
#'   runs the inner loop again; what fixes the answer is `eps_Omega` and
#'   `eps_fit`. Tightening it to `1e-6` on a 50 x 70 x 5 panel moved every
#'   coefficient by under `6e-4` and cost 66% more time, so the default trades
#'   a difference well below Monte Carlo error for a substantially cheaper fit.
#' @param eps_Omega,eps_fit Relative tolerances on the variance vector and the
#'   fitted values, outer loop. Both must be met.
#' @param eps_var Floor applied to every variance component before it is turned
#'   into a penalty, preventing division by a near-zero variance.
#' @param max_penalty Cap on each derived penalty. Needed for conditioning: a
#'   penalty near `1e8` pushes the block-tridiagonal systems past the limits of
#'   double precision and the solve fails. Reaching the cap means that variance
#'   component collapsed toward zero, and a warning says so.
#' @param delta Small constant guarding the denominators of the three
#'   convergence ratios.
#' @param accelerate How to drive the outer empirical-Bayes loop. `"none"` is
#'   the plain fixed-point iteration: it is what defines the estimate, and on a
#'   moderate panel it takes on the order of sixty passes, most of the cost of
#'   a fit. `"squarem"` extrapolates over that sequence (Varadhan and Roland,
#'   2008) to reach the same fixed point in fewer passes, then hands back to the
#'   plain iteration, which is what decides that the answer has been reached.
#'   The two agree to within the convergence tolerance by construction; if they
#'   ever disagree it is the acceleration that is wrong.
#' @param gauge Where to anchor the latent factors' scale convention. `U_t V_t'`
#'   is unchanged by `U_t -> c U_t`, `V_t -> V_t / c`, so only three functions
#'   of the four latent variance components are determined by the data and one
#'   convention has to be supplied. Which one is chosen decides which component
#'   is left holding the undetermined direction.
#'
#'   The asymmetry between the two node sets is carried by the ratio of
#'   dimensionless drifts `(tau_U^2 / sigma_U^2) / (tau_V^2 / sigma_V^2)`,
#'   which the rescaling leaves alone and which reads as how far the senders
#'   move per period relative to their own dispersion, against the same for the
#'   receivers. It is the same number under every setting below; they differ
#'   only in which component reports it.
#'
#'   `"innovation"` equalises the mean squared increments, so `tau_U^2 =
#'   tau_V^2` and the drift ratio is read off `sigma_V^2 / sigma_U^2`.
#'
#'   `"pooled"` and `"first"` equalise magnitudes instead — over all periods,
#'   or at `t = 1` — and leave the drift ratio in `tau_U^2 / tau_V^2`. That is
#'   where the empirical-Bayes loop is unstable: it sets `gamma_U =
#'   sigma_eps^2 / tau_U^2` from the innovation variance it has just measured,
#'   so a smaller `tau_U^2` smooths `U` further and lowers `tau_U^2` again. On
#'   panels generated with a true ratio of four the fitted ratio reached 68 and
#'   was still rising; on others it fell below `1e-3`. Both are kept for
#'   comparison, and `"pooled"` is what fits made before this argument existed
#'   used. See `.scale_normalize_UV()`.
#' @param latent_innovation Whether the latent block gets one innovation
#'   variance or two. `U` and `V` can trade temporal movement between them
#'   while leaving every `U_t V_t'` untouched, so the data determine only the
#'   product `tau_U^2 tau_V^2`; estimating the two separately leaves the
#'   undetermined split exposed to the empirical-Bayes feedback, which drove it
#'   past 68 on panels generated at four and below `1e-3` on others while the
#'   product came back to within twenty percent throughout. `"shared"` reports
#'   the geometric mean for both, which keeps what is determined and fixes the
#'   rest by convention. `"separate"` restores the older behaviour and exists
#'   so the drift can be measured; it is not a recommended way to fit.
#' @param seed Seed for the multistart perturbations; the caller's random-number
#'   state is restored on exit.
#' @param verbose Print progress.
#' @param panel_id Optional identifier stored in the result.
#'
#' @return An object of class `dynamic_ame`:
#'   \describe{
#'     \item{a, b, beta}{Trajectory matrices, `N x T`, `M x T`, `P x T`.}
#'     \item{U, V}{Lists of `T` latent factor matrices.}
#'     \item{results}{Per-period slices, each with `U`, `V`, `alpha` (the row
#'       effects `a_it`), `beta` (the column effects `b_jt`), `coef_row`,
#'       `coef_col`, `coef_dyad`. The field names follow the older per-period
#'       functions so the helpers in `post_estimation.R` apply directly — note
#'       that `beta` there means the column effects, not the regression
#'       coefficients.}
#'     \item{node_df, coef_df}{Tidy data frames, one row per node-period and per
#'       covariate-period.}
#'     \item{variance_components}{The eleven-entry `Omega`, plus the derived
#'       `lambda` and `gamma`.}
#'     \item{convergence}{Outer/inner iteration counts, flags, and traces.}
#'   }
#' @export
fit_dynamic_ame <- function(edge_panel,
                            row_cov_df = NULL,
                            col_cov_df = NULL,
                            dyad_cov_df = NULL,
                            row_covar_names = character(0),
                            col_covar_names = character(0),
                            dyad_covar_names = character(0),
                            time_col = "year",
                            row_col = "node_row",
                            col_col = "node_col",
                            value_col = "value",
                            row_prefix = "row_",
                            col_prefix = "col_",
                            row_pad = NULL,
                            col_pad = NULL,
                            transform = identity,
                            row_cov_time_col = time_col,
                            row_cov_id_col = row_col,
                            row_cov_prefix = row_prefix,
                            col_cov_time_col = time_col,
                            col_cov_id_col = col_col,
                            col_cov_prefix = col_prefix,
                            dyad_cov_time_col = time_col,
                            dyad_cov_row_id_col = row_col,
                            dyad_cov_col_id_col = col_col,
                            dyad_cov_row_prefix = row_prefix,
                            dyad_cov_col_prefix = col_prefix,
                            demean_row_covariates = FALSE,
                            demean_col_covariates = FALSE,
                            K = 2,
                            n_starts = 5,
                            perturb_sd = 0.1,
                            inner_max_iter = 500,
                            outer_max_iter = 100,
                            eps_Q = 1e-5,
                            eps_Omega = 1e-4,
                            eps_fit = 1e-4,
                            eps_var = 1e-8,
                            max_penalty = 1e6,
                            delta = 1e-8,
                            accelerate = c("squarem", "none"),
                            gauge = c("pooled", "first", "innovation"),
                            latent_innovation = c("shared", "separate"),
                            seed = 1,
                            verbose = FALSE,
                            panel_id = NA_character_) {
  accelerate <- match.arg(accelerate)
  gauge <- match.arg(gauge)
  latent_innovation <- match.arg(latent_innovation)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  # ── Data preparation (shared with the older estimator) ───────────────────
  panel <- build_panel_matrices(
    edge_panel = edge_panel, time_col = time_col, row_col = row_col,
    col_col = col_col, value_col = value_col, row_prefix = row_prefix,
    col_prefix = col_prefix, row_pad = row_pad, col_pad = col_pad,
    transform = transform
  )

  X_row_list <- build_row_covariates(
    row_cov_df = row_cov_df, years = panel$years, row_ids = panel$row_ids,
    covar_names = row_covar_names, time_col = row_cov_time_col,
    row_id_col = row_cov_id_col, row_prefix = row_cov_prefix
  )
  X_col_list <- build_col_covariates(
    col_cov_df = col_cov_df, years = panel$years, col_ids = panel$col_ids,
    covar_names = col_covar_names, time_col = col_cov_time_col,
    col_id_col = col_cov_id_col, col_prefix = col_cov_prefix
  )
  X_dyad_list <- build_dyad_covariates(
    dyad_cov_df = dyad_cov_df, years = panel$years, row_ids = panel$row_ids,
    col_ids = panel$col_ids, covar_names = dyad_covar_names,
    time_col = dyad_cov_time_col, row_id_col = dyad_cov_row_id_col,
    col_id_col = dyad_cov_col_id_col, row_prefix = dyad_cov_row_prefix,
    col_prefix = dyad_cov_col_prefix
  )

  Tn <- length(panel$years)

  # Within-node demeaning removes the part of a node covariate that is collinear
  # with that node's own additive effect.
  if (demean_row_covariates && length(row_covar_names) > 0 && Tn > 1 &&
      !is.null(X_row_list[[1]])) {
    mu <- Reduce("+", X_row_list) / Tn
    X_row_list <- lapply(X_row_list, function(X) {
      out <- X - mu; dimnames(out) <- dimnames(X); out
    })
    if (verbose) cat("  [within-node demeaning applied to row covariates]\n")
  }
  if (demean_col_covariates && length(col_covar_names) > 0 && Tn > 1 &&
      !is.null(X_col_list[[1]])) {
    mu <- Reduce("+", X_col_list) / Tn
    X_col_list <- lapply(X_col_list, function(X) {
      out <- X - mu; dimnames(out) <- dimnames(X); out
    })
    if (verbose) cat("  [within-node demeaning applied to col covariates]\n")
  }

  if (all(vapply(X_row_list, is.null, logical(1))))  X_row_list <- NULL
  if (all(vapply(X_col_list, is.null, logical(1))))  X_col_list <- NULL
  if (all(vapply(X_dyad_list, is.null, logical(1)))) X_dyad_list <- NULL

  pdims <- .cov_dims(
    X_row  = if (is.null(X_row_list))  NULL else X_row_list[[1]],
    X_col  = if (is.null(X_col_list))  NULL else X_col_list[[1]],
    X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[1]]
  )
  P <- pdims$P

  # ── Step 1: initialisation and multistart ────────────────────────────────
  init <- .init_ame_trajectory(
    Y_list = panel$Y_list, K = K, P = P,
    n_starts = n_starts, perturb_sd = perturb_sd, seed = seed
  )

  blocks <- c("beta", "a", "b", "U", "V")
  one <- stats::setNames(rep(1, length(blocks)), blocks)   # spec's neutral start

  if (verbose) cat(sprintf("Multistart: %d candidate(s)\n", n_starts))

  cand_Q <- rep(NA_real_, n_starts)
  cand_params <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    res <- try(
      withCallingHandlers(
        .ame_inner_bcd(init$candidates[[s]], panel$Y_list, X_row_list,
                       X_col_list, X_dyad_list, one, one,
                       max_iter = inner_max_iter, eps_Q = eps_Q, delta = delta),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      if (verbose) cat(sprintf("  candidate %d/%d failed\n", s, n_starts))
      next
    }
    cand_Q[s] <- res$objective$Q
    cand_params[[s]] <- res$params
    if (verbose) {
      cat(sprintf("  candidate %d/%d  Q = %.6f (%d sweeps)\n",
                  s, n_starts, res$objective$Q, res$iterations))
    }
  }

  if (all(is.na(cand_Q))) {
    stop("All multistart candidates failed.", call. = FALSE)
  }
  best <- which.min(cand_Q)
  if (verbose) cat(sprintf("  selected candidate %d (Q = %.6f)\n", best, cand_Q[best]))

  # ── Steps 2-10: outer empirical-Bayes loop ───────────────────────────────
  if (verbose) cat("Outer empirical-Bayes loop:\n")
  outer <- .ame_outer_eb(
    params = cand_params[[best]], reference = init$reference,
    Y_list = panel$Y_list, X_row_list = X_row_list, X_col_list = X_col_list,
    X_dyad_list = X_dyad_list, lambda = one, gamma = one,
    outer_max_iter = outer_max_iter, inner_max_iter = inner_max_iter,
    eps_Q = eps_Q, eps_Omega = eps_Omega, eps_fit = eps_fit,
    eps_var = eps_var, max_penalty = max_penalty, delta = delta,
    accelerate = accelerate, gauge_anchor = gauge,
    tie_latent = (latent_innovation == "shared"), verbose = verbose
  )

  .build_dynamic_ame(
    outer = outer, init = init, panel = panel,
    X_row_list = X_row_list, X_col_list = X_col_list, X_dyad_list = X_dyad_list,
    pdims = pdims,
    row_covar_names = row_covar_names, col_covar_names = col_covar_names,
    dyad_covar_names = dyad_covar_names,
    cand_Q = cand_Q, best = best, K = K, panel_id = panel_id,
    settings = list(
      K = K, n_starts = n_starts, perturb_sd = perturb_sd,
      inner_max_iter = inner_max_iter, outer_max_iter = outer_max_iter,
      eps_Q = eps_Q, eps_Omega = eps_Omega, eps_fit = eps_fit,
      eps_var = eps_var, max_penalty = max_penalty, delta = delta,
      accelerate = accelerate, gauge = gauge,
      latent_innovation = latent_innovation, seed = seed,
      demean_row_covariates = demean_row_covariates,
      demean_col_covariates = demean_col_covariates,
      time_col = time_col, row_col = row_col, col_col = col_col,
      value_col = value_col, row_prefix = row_prefix, col_prefix = col_prefix
    )
  )
}


#' Assemble the returned `dynamic_ame` object
#'
#' Splits the fitted trajectories into per-period slices whose field names match
#' the older per-period fit objects, so the helpers in `post_estimation.R` can be
#' used unchanged, and builds the tidy summaries.
#'
#' @keywords internal
#' @noRd
.build_dynamic_ame <- function(outer, init, panel,
                               X_row_list, X_col_list, X_dyad_list, pdims,
                               row_covar_names, col_covar_names, dyad_covar_names,
                               cand_Q, best, K, panel_id, settings) {
  p <- outer$params
  Tn <- length(panel$years)
  years <- panel$years

  idx_row  <- if (pdims$Pr > 0) seq_len(pdims$Pr) else integer(0)
  idx_col  <- if (pdims$Pc > 0) pdims$Pr + seq_len(pdims$Pc) else integer(0)
  idx_dyad <- if (pdims$Pd > 0) pdims$Pr + pdims$Pc + seq_len(pdims$Pd) else integer(0)

  results <- vector("list", Tn)
  names(results) <- as.character(years)
  node_rows <- vector("list", Tn)
  coef_rows <- vector("list", Tn)

  for (t in seq_len(Tn)) {
    beta_t <- if (nrow(p$beta) > 0L) p$beta[, t] else numeric(0)

    results[[t]] <- list(
      period    = years[t],
      U         = p$U[[t]],
      V         = p$V[[t]],
      alpha     = p$a[, t],          # row (sender) effects
      beta      = p$b[, t],          # column (receiver) effects -- NOT coefficients
      coef_row  = beta_t[idx_row],
      coef_col  = beta_t[idx_col],
      coef_dyad = beta_t[idx_dyad]
    )

    U_df <- as.data.frame(p$U[[t]]); names(U_df) <- paste0("X", seq_len(K))
    V_df <- as.data.frame(p$V[[t]]); names(V_df) <- paste0("X", seq_len(K))

    node_rows[[t]] <- dplyr::bind_rows(
      dplyr::bind_cols(
        data.frame(panel_id = panel_id, period = years[t],
                   node = panel$row_ids, mode = "row",
                   alpha = p$a[, t], beta = NA_real_,
                   stringsAsFactors = FALSE), U_df),
      dplyr::bind_cols(
        data.frame(panel_id = panel_id, period = years[t],
                   node = panel$col_ids, mode = "col",
                   alpha = NA_real_, beta = p$b[, t],
                   stringsAsFactors = FALSE), V_df)
    )

    parts <- list()
    if (pdims$Pr > 0) parts <- c(parts, list(data.frame(
      panel_id = panel_id, period = years[t], covariate = row_covar_names,
      coefficient = beta_t[idx_row], type = "row", stringsAsFactors = FALSE)))
    if (pdims$Pc > 0) parts <- c(parts, list(data.frame(
      panel_id = panel_id, period = years[t], covariate = col_covar_names,
      coefficient = beta_t[idx_col], type = "col", stringsAsFactors = FALSE)))
    if (pdims$Pd > 0) parts <- c(parts, list(data.frame(
      panel_id = panel_id, period = years[t], covariate = dyad_covar_names,
      coefficient = beta_t[idx_dyad], type = "dyad", stringsAsFactors = FALSE)))
    coef_rows[[t]] <- if (length(parts)) dplyr::bind_rows(parts) else NULL
  }

  rownames(p$a) <- panel$row_ids
  rownames(p$b) <- panel$col_ids
  colnames(p$a) <- as.character(years)
  colnames(p$b) <- as.character(years)
  if (nrow(p$beta) > 0L) {
    rownames(p$beta) <- c(row_covar_names, col_covar_names, dyad_covar_names)
    colnames(p$beta) <- as.character(years)
  }

  structure(
    list(
      panel_id = panel_id,
      years = years,
      row_ids = panel$row_ids,
      col_ids = panel$col_ids,
      a = p$a, b = p$b, beta = p$beta, U = p$U, V = p$V,
      results = results,
      node_df = dplyr::bind_rows(node_rows),
      coef_df = if (any(!vapply(coef_rows, is.null, logical(1))))
        dplyr::bind_rows(coef_rows) else NULL,
      variance_components = list(
        Omega = outer$Omega, lambda = outer$lambda, gamma = outer$gamma
      ),
      sigma_eps2 = outer$sigma_eps2,
      n_obs = outer$n_obs,
      reference = init$reference,
      convergence = list(
        converged = outer$converged,
        outer_iterations = outer$iterations,
        inner_iterations = outer$inner_iterations,
        objective_trace = outer$objective_trace,
        Omega_trace = outer$Omega_trace,
        scale_trace = outer$scale_trace,
        # What the unconstrained split would have been at each pass. Kept even
        # when the latent innovation variances are tied: it is the direction
        # the data push a quantity they do not determine, and seeing it is how
        # a user learns the tie is doing work rather than sitting idle.
        latent_ratio_trace = outer$latent_ratio_trace,
        multistart_Q = cand_Q,
        multistart_selected = best,
        penalties_capped = outer$penalties_capped
      ),
      covariates = list(
        X_row_list = X_row_list, X_col_list = X_col_list,
        X_dyad_list = X_dyad_list,
        row_covar_names = row_covar_names, col_covar_names = col_covar_names,
        dyad_covar_names = dyad_covar_names, dims = pdims
      ),
      Y_list = panel$Y_list,
      settings = settings
    ),
    class = "dynamic_ame"
  )
}
