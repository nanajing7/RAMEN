# Public entry point for inference on a dynamic AME fit.

#' Parametric bootstrap inference for a dynamic bipartite AME model
#'
#' Runs either or both of two parametric bootstrap designs. They answer
#' different questions and report different quantities.
#'
#' The CONDITIONAL design is the default, because it is the one that produces
#' standard errors and confidence intervals -- the quantities an analysis
#' reports. The full-model design measures how well the estimator recovers
#' known generating values; its dispersion is governed by how widely those
#' generating values vary and is NOT a standard error. It answers a question
#' about the estimator rather than about this dataset, so it is run only when
#' asked for.
#'
#' @section Design "conditional" — how precise are these estimates?:
#' The fitted systematic component is held fixed and only the observation errors
#' are regenerated, `eps ~ N(0, sigma_eps2)`, on the observed cells. Every
#' replicate therefore carries exactly the dependence structure the model
#' estimated — sender and receiver heterogeneity through the additive and latent
#' effects, temporal persistence through the random-walk priors — and the
#' dependence is preserved by construction rather than by the resampling scheme,
#' which is why the errors can be drawn independently across cells.
#'
#' Reports standard errors and percentile confidence intervals for the
#' coefficient trajectories, the additive effects, the variance components, and
#' the interaction matrices.
#'
#' @section Design "full" — how well does the estimator work?:
#' Every state trajectory is redrawn from the priors that the fitted variance
#' components parameterise and a complete new panel is generated, so each
#' replicate is a fresh realization of the model with a *known* set of
#' generating parameters.
#'
#' Because those parameters differ from replicate to replicate, the dispersion
#' of the estimates is not informative — it is governed by how widely the
#' generating values vary, and an estimator that recovered every replicate
#' exactly would show the same dispersion. What is reported instead is the
#' difference between each estimate and the value that generated it, summarised
#' by the mean error (systematic bias) and the root mean squared error.
#'
#' @section What is held fixed in both designs:
#' Nodes, periods, and covariates. The model specifies the conditional
#' distribution of the outcome given the covariates and does not model the
#' covariates themselves; the node sets and time window are conditioned on
#' rather than treated as draws from a larger population. Dyad-periods that are
#' unobserved in the data remain unobserved in every replicate, so each
#' replicate carries the same information as the original panel. The latent
#' dimension `K` is held at the value used for the original fit, so results do
#' not reflect uncertainty about its selection.
#'
#' @section Latent factors:
#' `U_t` and `V_t` are identified only up to a common orthogonal transformation
#' and a reciprocal rescaling, so their individual elements are not comparable
#' across replicates and are never summarised. Their product `U_t V_t'` is
#' invariant to both and is used throughout. For the same reason the full-model
#' design reports `sigma_U2`, `sigma_V2`, `tau_U2`, and `tau_V2` only through
#' the gauge-invariant products.
#'
#' @section Cost:
#' Every replicate re-estimates the model in full, including the empirical-Bayes
#' variance-component update, so the penalties are re-derived within each
#' replicate rather than held at their original values. The conditional design
#' warm-starts each replicate at the original estimate, which is both cheaper
#' and keeps replicates within one basin of attraction; the full-model design
#' cannot, since its generating values differ, and runs the multistart in full.
#' Expect the full-model design to be several times slower per replicate.
#' Replications are independent — set `n_cores` to use them.
#'
#' @param fit A `dynamic_ame` object from [fit_dynamic_ame()].
#' @param design Which designs to run: `"conditional"` (default), `"full"`, or
#'   `"both"`.
#' @param B Replications for the conditional design.
#' @param B_full Replications for the full-model design. Defaults to `B`; it can
#'   usually be smaller, since the full-model design reports means and root mean
#'   squared errors rather than tail quantiles.
#' @param conf_level Confidence level for the conditional design's intervals.
#' @param warm_start Conditional design only: initialise each replicate at the
#'   original estimate.
#' @param dispersion Conditional design only. A fitted trajectory is a posterior
#'   mean and so moves less than the parameter it estimates; generating panels
#'   from it alone asks the estimator to track a path that is already smooth,
#'   which understates the standard errors. `"posterior"` (the default) restores
#'   the missing dispersion for `beta`, `a` and `b` before each panel is
#'   generated, using the same `sigma^2` and `tau^2` the model reports.
#'   `"none"` generates from the fitted values as they stand. `U` and `V` are
#'   never perturbed: `tau_U^2` and `tau_V^2` can be traded between the two
#'   factors without changing the fit, so they describe the identification
#'   convention rather than the data, and the intervals stay conservative in the
#'   part of the error that comes from smoothing the latent factors.
#' @param store_latent_draws Conditional design only: keep every draw of
#'   `U_t V_t'`. Required for exact percentile intervals on the interaction
#'   matrices; otherwise a normal approximation is used and flagged as such.
#'   Memory grows as `B * T * N * M`.
#' @param n_cores Parallel workers. Each replicate is given its own seed up
#'   front, so results do not depend on the number of workers or on scheduling.
#' @param seed Base seed. The two designs are given separate derived seeds.
#' @param verbose Report progress.
#'
#' @return An object of class `bootstrap_ame` with elements `conditional` and
#'   `full` (either may be `NULL` if not requested). Each can be printed on its
#'   own; printing the whole object summarises both.
#'
#' @examples
#' \dontrun{
#' fit <- fit_dynamic_ame(edge_panel = my_panel, K = 2)
#'
#' # inference: standard errors and confidence intervals
#' bc <- bootstrap_ame(fit, B = 1000, n_cores = 4, store_latent_draws = TRUE)
#' bc$conditional$summaries$beta       # estimate, SE, CI per period
#' bc$conditional$latent$sign_prob     # Pr(U_t V_t' > 0) per dyad
#'
#' # both designs, when the estimator study is wanted as well
#' bt <- bootstrap_ame(fit, design = "both", B = 1000, B_full = 500,
#'                     n_cores = 4, store_latent_draws = TRUE)
#'
#' bt                                  # both designs side by side
#' bt$full$accuracy$Omega              # bias and RMSE for the variances
#' bt$full$latent_recovery             # recovery of U_t V_t' by period
#' }
#' @export
bootstrap_ame <- function(fit,
                          design = c("conditional", "full", "both"),
                          B = 1000,
                          B_full = B,
                          conf_level = 0.95,
                          warm_start = TRUE,
                          dispersion = c("posterior", "none"),
                          store_latent_draws = FALSE,
                          n_cores = 1,
                          seed = 1,
                          verbose = FALSE) {
  .validate_dynamic_ame(fit)
  design <- match.arg(design)
  dispersion <- match.arg(dispersion)

  # Separate derived seeds, so the two designs are not driven by the same
  # random numbers, while the whole call stays reproducible from `seed`.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  design_seeds <- sample.int(.Machine$integer.max, 2)

  out_conditional <- NULL
  out_full <- NULL

  if (design %in% c("both", "conditional")) {
    if (verbose) cat("== Conditional design ==\n")
    out_conditional <- .bootstrap_conditional(
      fit, B = B, conf_level = conf_level, warm_start = warm_start,
      dispersion = dispersion,
      store_latent_draws = store_latent_draws, n_cores = n_cores,
      seed = design_seeds[1], verbose = verbose
    )
  }

  if (design %in% c("both", "full")) {
    if (verbose) cat("== Full-model design ==\n")
    out_full <- .bootstrap_full(
      fit, B = B_full, n_cores = n_cores,
      seed = design_seeds[2], verbose = verbose
    )
  }

  structure(
    list(
      conditional = out_conditional,
      full = out_full,
      design = design,
      fit = fit,
      settings = list(B = B, B_full = B_full, conf_level = conf_level,
                      warm_start = warm_start, dispersion = dispersion,
                      store_latent_draws = store_latent_draws,
                      n_cores = n_cores, seed = seed)
    ),
    class = "bootstrap_ame"
  )
}


#' Print a combined bootstrap result
#'
#' @param x A `bootstrap_ame` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.bootstrap_ame <- function(x, ...) {
  cat("Parametric bootstrap inference for a dynamic bipartite AME model\n")
  cat(strrep("=", 66), "\n", sep = "")

  if (!is.null(x$conditional)) {
    cat("\n[1] Conditional design - precision of the estimates for this network\n\n")
    print(x$conditional)
  }

  if (!is.null(x$full)) {
    cat("\n", strrep("-", 66), "\n", sep = "")
    cat("\n[2] Full-model design - accuracy of the estimator under the model\n\n")
    print(x$full)
  }

  if (!is.null(x$conditional) && !is.null(x$full)) {
    cat("\n", strrep("-", 66), "\n", sep = "")
    cat("\nThe two designs report different quantities on purpose. The\n")
    cat("conditional design gives standard errors and intervals for the\n")
    cat("observed network; the full-model design gives bias and RMSE against\n")
    cat("known generating values. Its estimate dispersion is not reported,\n")
    cat("since it would reflect the variation of those generating values\n")
    cat("rather than the performance of the estimator.\n")
  }

  invisible(x)
}
