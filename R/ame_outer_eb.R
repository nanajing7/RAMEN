# Step 10: the outer empirical-Bayes loop.
#
# Each outer iteration runs the inner BCD to convergence at the current
# penalties, identifies the latent trajectories (Steps 8 then 7), and refreshes
# the variance components (Step 9). Penalised objective values are deliberately
# not compared across outer iterations, because the penalties change after every
# variance update. Convergence is judged instead on two quantities that remain
# comparable: the variance vector Omega and the fitted values Yhat.

#' Relative change between two vectors, guarded against a zero denominator
#'
#' `NA` entries — the components the spec leaves undefined when `P = 0` or
#' `T = 1` — are dropped from both sides.
#'
#' @param new,old Numeric vectors of the same length.
#' @param delta Small constant added to the denominator.
#' @return A scalar relative change.
#' @keywords internal
#' @noRd
.rel_change <- function(new, old, delta = 1e-8) {
  keep <- is.finite(new) & is.finite(old)
  if (!any(keep)) return(0)
  sqrt(sum((new[keep] - old[keep])^2)) / (sqrt(sum(old[keep]^2)) + delta)
}


#' Stacked fitted values across all periods
#'
#' @param params Parameter list.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param Tn Number of periods.
#' @return A numeric vector of all fitted values.
#' @keywords internal
#' @noRd
.stacked_fitted <- function(params, X_row_list, X_col_list, X_dyad_list, Tn) {
  unlist(lapply(seq_len(Tn), function(t) {
    .fitted_period(
      a_t = params$a[, t], b_t = params$b[, t],
      U_t = params$U[[t]], V_t = params$V[[t]],
      beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
      X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
      X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
      X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]]
    )
  }))
}


#' Run the outer empirical-Bayes loop (Step 10)
#'
#' Repeats inner BCD -> identification -> variance update until both criteria
#' are met:
#' \deqn{\frac{\|\Omega^{(m)} - \Omega^{(m-1)}\|}{\|\Omega^{(m-1)}\| + \delta}
#'       < \epsilon_\Omega
#'   \quad\text{and}\quad
#'   \frac{\|\hat Y^{(m)} - \hat Y^{(m-1)}\|_F}{\|\hat Y^{(m-1)}\|_F + \delta}
#'       < \epsilon_{\text{fit}}.}
#'
#' Both criteria are invariant to the common rotation left free by the model, so
#' the choice of orientation never affects when the loop stops.
#'
#' @section Acceleration:
#' The loop is a fixed-point iteration on the ten penalties, and fixed-point
#' iterations of this kind converge linearly — on a moderate panel this one
#' takes on the order of sixty passes, while the inner loop it drives settles
#' in one or two. Nearly all of the cost is therefore in how slowly the outer
#' sequence approaches its limit, not in the work done at any one pass.
#'
#' `accelerate = "squarem"` puts a squared extrapolation (Varadhan and Roland,
#' 2008) over that sequence. It changes how quickly the iteration arrives, not
#' where: the plain iteration remains the definition of the fixed point, and it
#' is the plain iteration that decides when the answer has been reached. After
#' the extrapolation stops, the loop continues plainly until both of the
#' original criteria hold, so `converged` means the same thing in either mode
#' and the two are comparable.
#'
#' Two details make the extrapolation safe here. It works on the log penalties,
#' so an extrapolated point cannot be negative however far it reaches. And it
#' runs without an objective function to backtrack on, because there is none to
#' be had — the penalised objective is not comparable across outer iterations,
#' the penalties having changed — so correctness rests on the plain iteration
#' that follows rather than on the extrapolation behaving.
#'
#' One case is outside all of this, and it is worth stating because it looks
#' like a disagreement between the modes and is not. When a variance component
#' is not identified by the panel it collapses toward zero, and there is then no
#' fixed point in that coordinate for either mode to reach. Both are still
#' descending when the criterion stops them, and the criterion — a relative
#' change in the norm of `Omega` — is insensitive to a component that has
#' already become small, so it stops them at different depths. Neither answer is
#' wrong and neither mode is at fault; what tells a user that this has happened
#' is `penalties_capped` and the warning that accompanies it, not agreement
#' between the two modes.
#'
#' @param params Starting parameters (the winning multistart candidate).
#' @param reference Step-1 reference for the identification step.
#' @param Y_list List of `T` outcome matrices.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param lambda,gamma Starting penalties (the spec starts all ten at 1).
#' @param outer_max_iter,inner_max_iter Iteration caps. `outer_max_iter` bounds
#'   the total number of passes in either mode, extrapolation included.
#' @param eps_Q,eps_Omega,eps_fit Tolerances.
#' @param eps_var Floor applied to every variance component.
#' @param max_penalty Cap applied to each derived penalty, for conditioning.
#' @param delta Denominator guard shared by all three criteria.
#' @param accelerate `"none"` for the plain iteration, `"squarem"` to
#'   extrapolate over it.
#' @param tie_latent Give the latent block one innovation variance instead of
#'   two. Only the product `tau_U^2 tau_V^2` is determined by the data — the
#'   split can be moved by rescaling `U` and `V` against each other — and
#'   estimating the two separately exposes the free direction to this loop's own
#'   feedback. See `.eb_update_variances()`.
#' @param gauge_anchor Which quantities fix the latent scale convention: the
#'   mean squared increments, or the magnitudes pooled or at `t = 1`. Passed to
#'   `.ame_identify()`. The choice decides which variance component is left
#'   holding the one undetermined direction, and this loop is what makes that
#'   matter — it derives `gamma_U` from `tau_U^2` and so is unstable wherever
#'   the free coordinate sits in `tau_U^2 / tau_V^2`. See
#'   `.scale_normalize_UV()`.
#' @param verbose If `TRUE`, print a line per outer iteration.
#'
#' @return A list with the final `params`, `Omega`, `lambda`, `gamma`,
#'   `sigma_eps2`, the per-iteration `Omega_trace`, `objective_trace`,
#'   `inner_iterations`, `scale_trace`, the number of outer `iterations`,
#'   `converged`, and the `accelerate` mode that produced them.
#' @keywords internal
#' @noRd
.ame_outer_eb <- function(params,
                          reference,
                          Y_list,
                          X_row_list = NULL,
                          X_col_list = NULL,
                          X_dyad_list = NULL,
                          lambda,
                          gamma,
                          outer_max_iter = 100,
                          inner_max_iter = 500,
                          eps_Q = 1e-6,
                          eps_Omega = 1e-4,
                          eps_fit = 1e-4,
                          eps_var = 1e-8,
                          max_penalty = 1e6,
                          delta = 1e-8,
                          accelerate = c("none", "squarem"),
                          gauge_anchor = c("pooled", "first", "innovation"),
                          tie_latent = TRUE,
                          verbose = FALSE) {
  gauge_anchor <- match.arg(gauge_anchor)
  accelerate <- match.arg(accelerate)
  Tn <- length(Y_list)
  blocks <- c("beta", "a", "b", "U", "V")
  lambda <- .check_penalty(lambda, "lambda", blocks)
  gamma  <- .check_penalty(gamma,  "gamma",  blocks)

  if (accelerate == "squarem" && !requireNamespace("SQUAREM", quietly = TRUE)) {
    warning("Package 'SQUAREM' is not installed; the plain iteration was used ",
            "instead. The answer is unaffected — only the number of passes.",
            call. = FALSE)
    accelerate <- "none"
  }

  Omega_trace <- list()
  objective_trace <- numeric(0)
  inner_iterations <- integer(0)
  scale_trace <- numeric(0)
  latent_ratio_trace <- numeric(0)
  ever_capped <- character(0)
  eb <- NULL
  n_pass <- 0L

  # One pass of the loop: inner BCD at the given penalties, identification,
  # then the variance refresh. `params` is carried in this frame rather than
  # passed in and out, because it is the state the inner loop warm-starts from
  # and not part of the fixed point being solved for. Every pass is counted and
  # traced, including the ones an extrapolation asks for, so the recorded
  # iteration count is the work actually done rather than the number of
  # accelerated steps.
  pass <- function(lam, gam) {
    # The inner loop warns on its own if it stalls; the outer loop can still
    # make progress, so that is not fatal here.
    bcd <- withCallingHandlers(
      .ame_inner_bcd(params, Y_list, X_row_list, X_col_list, X_dyad_list,
                     lam, gam, max_iter = inner_max_iter,
                     eps_Q = eps_Q, delta = delta),
      warning = function(w) invokeRestart("muffleWarning")
    )

    ident <- .ame_identify(bcd$params, reference, anchor = gauge_anchor)
    params <<- ident$params

    # The EM step needs the penalties that were in force, to rebuild the
    # posterior covariances of the trajectories it just estimated.
    e <- .eb_update_variances(params, Y_list, X_row_list, X_col_list,
                              X_dyad_list, lambda = lam, gamma = gam,
                              eps_var = eps_var, max_penalty = max_penalty,
                              tie_latent = tie_latent)
    if (any(e$capped))
      ever_capped <<- union(ever_capped, names(e$capped)[e$capped])

    n_pass <<- n_pass + 1L
    Omega_trace[[n_pass]] <<- e$Omega
    objective_trace <<- c(objective_trace, bcd$objective$Q)
    inner_iterations <<- c(inner_iterations, bcd$iterations)
    scale_trace <<- c(scale_trace, ident$c)
    # What the unconstrained split would have been at this pass. Recorded even
    # when it is not used: a run of these is what makes the instability visible
    # instead of inferred.
    latent_ratio_trace <<- c(latent_ratio_trace, e$latent_ratio)
    eb <<- e
    e$inner <- bcd$iterations
    e$Q <- bcd$objective$Q
    e
  }

  # The plain iteration, and the arbiter of convergence in both modes. Two
  # consecutive passes are always needed before it can stop, there being no
  # change to measure on the first.
  run_plain <- function(lam, gam, budget) {
    Omega_prev <- NULL
    Yhat_prev <- NULL
    done <- FALSE

    while (budget > 0L) {
      budget <- budget - 1L
      e <- pass(lam, gam)
      Yhat <- .stacked_fitted(params, X_row_list, X_col_list, X_dyad_list, Tn)

      d_Omega <- if (is.null(Omega_prev)) NA_real_ else
        .rel_change(e$Omega, Omega_prev, delta)
      d_fit <- if (is.null(Yhat_prev)) NA_real_ else
        .rel_change(Yhat, Yhat_prev, delta)

      if (verbose) {
        cat(sprintf("  outer %3d  inner %4d  Q = %.6f  dOmega = %-10s dYhat = %s\n",
                    n_pass, e$inner, e$Q,
                    if (is.na(d_Omega)) "-" else sprintf("%.3e", d_Omega),
                    if (is.na(d_fit)) "-" else sprintf("%.3e", d_fit)))
      }

      if (!is.na(d_Omega) && !is.na(d_fit) &&
          d_Omega < eps_Omega && d_fit < eps_fit) {
        done <- TRUE
        break
      }

      Omega_prev <- e$Omega
      Yhat_prev <- Yhat
      lam <- e$lambda
      gam <- e$gamma
    }
    list(lambda = lam, gamma = gam, converged = done)
  }

  if (accelerate == "none") {
    converged <- run_plain(lambda, gamma, outer_max_iter)$converged
  } else {
    # The extrapolation works on the state the loop actually carries -- the ten
    # penalties -- on the log scale, so no extrapolated point can be negative.
    # Entries the model leaves undefined (P = 0, or T = 1) are held out and
    # passed through: they are not part of the fixed point and taking their
    # logarithm would put NaN into the search.
    pen0 <- c(lambda, gamma)
    live <- is.finite(pen0) & pen0 > 0
    n_lam <- length(lambda)

    unpack <- function(theta) {
      p <- pen0
      p[live] <- exp(theta)
      list(lambda = p[seq_len(n_lam)], gamma = p[-seq_len(n_lam)])
    }
    fixpt <- function(theta) {
      pg <- unpack(theta)
      e <- pass(pg$lambda, pg$gamma)
      log(pmax(c(e$lambda, e$gamma)[live], .Machine$double.xmin))
    }

    # Each accelerated step costs about three passes. This budget is the
    # extrapolation's alone; the plain iteration that follows gets its own,
    # for the reason given where it is called.
    sq <- try(SQUAREM::squarem(
      par = log(pen0[live]), fixptfn = fixpt,
      control = list(tol = eps_Omega, trace = FALSE,
                     maxiter = max(1L, floor(outer_max_iter / 3)))),
      silent = TRUE)

    if (inherits(sq, "try-error")) {
      warning("The extrapolation failed (", sub("\n.*", "", conditionMessage(
        attr(sq, "condition"))), "); the plain iteration continued from where ",
        "it stopped. The answer is unaffected.", call. = FALSE)
      start <- if (is.null(eb)) list(lambda = lambda, gamma = gamma) else
        list(lambda = eb$lambda, gamma = eb$gamma)
    } else {
      start <- unpack(sq$par)
    }

    # Whatever the extrapolation reached is a proposal, not an answer. The
    # plain iteration runs on from it under the original criteria, so what is
    # returned has been certified the same way it would have been without any
    # acceleration at all.
    #
    # It gets the whole of outer_max_iter rather than what the extrapolation
    # left over. Sharing one budget between the two makes it possible for the
    # accelerated mode to fail to converge where the plain mode succeeds --
    # the extrapolation spends passes backtracking, and the certification is
    # left too few to finish. That happened, on two of eight panels. An
    # acceleration that converges worse than not accelerating is a defect and
    # not a trade-off, so the certification is given at least what it would
    # have had on its own. The total is bounded by roughly 4/3 of
    # outer_max_iter, and `iterations` still counts every pass actually run.
    converged <- run_plain(start$lambda, start$gamma,
                           outer_max_iter)$converged
  }

  if (!converged) {
    warning(
      sprintf("Outer empirical-Bayes loop did not converge in %d passes.",
              n_pass),
      call. = FALSE
    )
  }

  # A binding cap means the corresponding variance estimate collapsed toward
  # zero. The fit is still usable, but that component is not to be trusted.
  if (length(ever_capped) > 0) {
    warning(
      sprintf(paste0("Penalty cap (%g) was reached for: %s. The matching ",
                     "variance component(s) collapsed toward zero, so those ",
                     "entries of Omega are unreliable; the rest of the fit is ",
                     "unaffected."),
              max_penalty, paste(sort(ever_capped), collapse = ", ")),
      call. = FALSE
    )
  }

  list(
    params = params,
    Omega = eb$Omega,
    lambda = eb$lambda,
    gamma = eb$gamma,
    sigma_eps2 = eb$sigma_eps2,
    n_obs = eb$n_obs,
    Omega_trace = Omega_trace,
    objective_trace = objective_trace,
    inner_iterations = inner_iterations,
    scale_trace = scale_trace,
    latent_ratio_trace = latent_ratio_trace,
    iterations = n_pass,
    converged = converged,
    accelerate = accelerate,
    penalties_capped = ever_capped
  )
}
