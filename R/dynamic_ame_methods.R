# Methods and accessors for objects returned by fit_dynamic_ame().

#' Validate a fitted dynamic AME object
#'
#' @param object An object expected to inherit from `dynamic_ame`.
#' @return Invisibly `TRUE`.
#' @keywords internal
#' @noRd
.validate_dynamic_ame <- function(object) {
  if (!inherits(object, "dynamic_ame")) {
    stop("object must be of class 'dynamic_ame'.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Dimensions of a fitted dynamic AME model
#'
#' @param object A `dynamic_ame` object.
#' @return A named integer vector with `N`, `M`, `K`, `T`, and `P`.
#' @keywords internal
#' @noRd
.dynamic_ame_dims <- function(object) {
  c(N = nrow(object$a), M = nrow(object$b), K = ncol(object$U[[1]]),
    T = length(object$years), P = nrow(object$beta))
}


#' Print a fitted dynamic AME model
#'
#' @param x A `dynamic_ame` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.dynamic_ame <- function(x, ...) {
  .validate_dynamic_ame(x)
  d <- .dynamic_ame_dims(x)

  cat("Dynamic bipartite AME model (joint-trajectory MAP + empirical Bayes)\n")
  if (!is.na(x$panel_id)) cat("Panel:", x$panel_id, "\n")
  cat(sprintf("%d senders x %d receivers, %d periods (%s-%s), K = %d\n",
              d[["N"]], d[["M"]], d[["T"]],
              as.character(x$years[1]), as.character(x$years[d[["T"]]]),
              d[["K"]]))
  cat(sprintf("Covariates: %d  |  observed cells: %d of %d (%.1f%%)\n",
              d[["P"]], x$n_obs, d[["N"]] * d[["M"]] * d[["T"]],
              100 * x$n_obs / (d[["N"]] * d[["M"]] * d[["T"]])))

  cv <- x$convergence
  cat(sprintf("\nConverged: %s after %d outer iteration%s (%d inner sweeps total)\n",
              if (isTRUE(cv$converged)) "yes" else "NO",
              cv$outer_iterations,
              if (cv$outer_iterations == 1L) "" else "s",
              sum(cv$inner_iterations)))
  if (!isTRUE(cv$converged)) {
    cat("  the outer loop hit its iteration cap; treat the estimates with care\n")
  }
  if (length(cv$penalties_capped) > 0) {
    cat(sprintf("  penalty cap reached for: %s\n",
                paste(sort(cv$penalties_capped), collapse = ", ")))
  }

  cat(sprintf("\nResidual variance (sigma_eps^2): %.6f\n", x$sigma_eps2))
  cat("Variance components:\n")
  print(round(x$variance_components$Omega, 6))

  invisible(x)
}


#' Summarise a fitted dynamic AME model
#'
#' Reports the variance components alongside the penalties they imply, the
#' convergence history, and how much of each variance estimate came from the
#' posterior-variance term of the EM step — a large share means the data leave
#' that component only weakly determined.
#'
#' @param object A `dynamic_ame` object.
#' @param ... Ignored.
#' @return An object of class `summary.dynamic_ame`.
#' @export
summary.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  d <- .dynamic_ame_dims(object)
  vc <- object$variance_components

  blocks <- c("beta", "a", "b", "U", "V")
  comp <- data.frame(
    block = blocks,
    sigma2 = unname(vc$Omega[paste0("sigma_", blocks, "2")]),
    tau2 = unname(vc$Omega[paste0("tau_", blocks, "2")]),
    lambda = unname(vc$lambda[blocks]),
    gamma = unname(vc$gamma[blocks]),
    stringsAsFactors = FALSE
  )

  coef_summary <- NULL
  if (d[["P"]] > 0L) {
    coef_summary <- data.frame(
      covariate = rownames(object$beta),
      mean = rowMeans(object$beta),
      first = object$beta[, 1],
      last = object$beta[, d[["T"]]],
      sd_over_time = apply(object$beta, 1, stats::sd),
      stringsAsFactors = FALSE
    )
    rownames(coef_summary) <- NULL
  }

  # The latent innovation variances are tied when the fit asked for it, so they
  # are one estimate reported twice rather than two. What the data would have
  # made of the split is carried alongside, because a tie that is never shown
  # to bind looks like an assumption rather than a finding.
  lrt <- object$convergence$latent_ratio_trace
  latent <- list(
    shared = isTRUE(object$settings$latent_innovation == "shared"),
    kappa = unname(vc$Omega[["tau_U2"]]),
    ratio_unconstrained = if (length(lrt)) unname(utils::tail(lrt, 1)) else
      NA_real_
  )

  structure(
    list(
      panel_id = object$panel_id,
      dims = d,
      years = object$years,
      n_obs = object$n_obs,
      sigma_eps2 = object$sigma_eps2,
      components = comp,
      latent = latent,
      coef_summary = coef_summary,
      convergence = object$convergence
    ),
    class = "summary.dynamic_ame"
  )
}


#' Print a dynamic AME summary
#'
#' @param x A `summary.dynamic_ame` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.summary.dynamic_ame <- function(x, digits = 5, ...) {
  d <- x$dims
  cat("Dynamic bipartite AME model\n")
  if (!is.na(x$panel_id)) cat("Panel:", x$panel_id, "\n")
  cat(sprintf("%d x %d nodes, %d periods, K = %d, P = %d\n",
              d[["N"]], d[["M"]], d[["T"]], d[["K"]], d[["P"]]))
  cat(sprintf("Observed cells: %d\n", x$n_obs))

  cat(sprintf("\nResidual variance: %s\n", format(x$sigma_eps2, digits = digits)))
  cat("\nVariance components and the penalties they imply:\n")
  print(format(x$components, digits = digits), row.names = FALSE)
  cat("\n  sigma2 = initial-state variance, tau2 = innovation variance\n")
  cat("  lambda = sigma_eps2/sigma2, gamma = sigma_eps2/tau2\n")
  cat("  NA means the component is undefined (no covariates, or T = 1)\n")

  lt <- x$latent
  if (!is.null(lt) && isTRUE(lt$shared)) {
    cat(sprintf(
      "\n  U and V share one innovation variance, kappa = %s.\n",
      format(lt$kappa, digits = digits)))
    cat("  Only the product tau_U2 tau_V2 is determined by the data: rescaling\n")
    cat("  U against V moves the split without touching any U_t V_t'. The two\n")
    cat("  rows above are one estimate reported twice, not two estimates.\n")
    if (is.finite(lt$ratio_unconstrained)) {
      cat(sprintf(
        "  Estimated separately the last pass would have put the ratio at %s;\n",
        format(lt$ratio_unconstrained, digits = 3)))
      cat("  that number is a property of where the iteration was stopped, not\n")
      cat("  of the data, and is shown only so the constraint is visible.\n")
    }
  }

  if (!is.null(x$coef_summary)) {
    cat("\nCoefficient trajectories:\n")
    print(format(x$coef_summary, digits = digits), row.names = FALSE)
  }

  cv <- x$convergence
  cat(sprintf("\nConvergence: %s, %d outer iterations\n",
              if (isTRUE(cv$converged)) "reached" else "NOT reached",
              cv$outer_iterations))
  cat(sprintf("  inner sweeps per outer iteration: min %d, median %d, max %d\n",
              min(cv$inner_iterations), round(stats::median(cv$inner_iterations)),
              max(cv$inner_iterations)))
  if (length(cv$multistart_Q) > 1L) {
    cat(sprintf("  multistart: %d candidates, kept #%d (Q = %s)\n",
                length(cv$multistart_Q), cv$multistart_selected,
                format(min(cv$multistart_Q, na.rm = TRUE), digits = digits)))
  }
  if (length(cv$penalties_capped) > 0) {
    cat(sprintf("  penalty cap reached for: %s\n",
                paste(sort(cv$penalties_capped), collapse = ", ")))
  }

  invisible(x)
}


#' Fitted values from a dynamic AME model
#'
#' @param object A `dynamic_ame` object.
#' @param ... Ignored.
#' @return A list of `T` matrices of fitted values, named by period.
#' @export
fitted.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  Tn <- length(object$years)
  cv <- object$covariates

  out <- lapply(seq_len(Tn), function(t) {
    f <- .fitted_period(
      a_t = object$a[, t], b_t = object$b[, t],
      U_t = object$U[[t]], V_t = object$V[[t]],
      beta_t = if (nrow(object$beta) > 0L) object$beta[, t] else numeric(0),
      X_row  = if (is.null(cv$X_row_list))  NULL else cv$X_row_list[[t]],
      X_col  = if (is.null(cv$X_col_list))  NULL else cv$X_col_list[[t]],
      X_dyad = if (is.null(cv$X_dyad_list)) NULL else cv$X_dyad_list[[t]]
    )
    dimnames(f) <- list(object$row_ids, object$col_ids)
    f
  })
  names(out) <- as.character(object$years)
  out
}


#' Residuals from a dynamic AME model
#'
#' Cells that were unobserved, or whose covariates were incomplete, stay `NA` —
#' they took no part in the fit.
#'
#' @param object A `dynamic_ame` object.
#' @param ... Ignored.
#' @return A list of `T` residual matrices, named by period.
#' @export
residuals.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  fit <- fitted(object)
  out <- lapply(seq_along(fit), function(t) object$Y_list[[t]] - fit[[t]])
  names(out) <- names(fit)
  out
}


#' Extract the estimated variance components
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return For `dynamic_ame`, a list with `Omega` (the eleven variance
#'   components), and the derived penalties `lambda` and `gamma`.
#' @export
variance_components <- function(object, ...) UseMethod("variance_components")

#' @rdname variance_components
#' @export
variance_components.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  object$variance_components
}


#' Extract the regression-coefficient trajectory
#'
#' @param object A fitted model object.
#' @param long If `TRUE`, return a long data frame instead of the `P x T`
#'   matrix.
#' @param ... Passed to methods.
#' @return A `P x T` matrix, or a long data frame with `covariate`, `period`,
#'   and `coefficient`. `NULL` when the model has no covariates.
#' @export
coef_trajectory <- function(object, long = FALSE, ...) UseMethod("coef_trajectory")

#' @rdname coef_trajectory
#' @export
coef_trajectory.dynamic_ame <- function(object, long = FALSE, ...) {
  .validate_dynamic_ame(object)
  if (nrow(object$beta) == 0L) return(NULL)
  if (!long) return(object$beta)
  object$coef_df
}


#' Extract the fitted multiplicative component
#'
#' Returns `U_t V_t'` for each period. This is the part of the latent structure
#' the model actually identifies: the individual factors are defined only up to
#' the gauge the fit fixes, but their product is invariant.
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return A list of `T` matrices, named by period.
#' @export
latent_trajectory <- function(object, ...) UseMethod("latent_trajectory")

#' @rdname latent_trajectory
#' @export
latent_trajectory.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  out <- lapply(seq_along(object$years), function(t) {
    m <- object$U[[t]] %*% t(object$V[[t]])
    dimnames(m) <- list(object$row_ids, object$col_ids)
    m
  })
  names(out) <- as.character(object$years)
  out
}


#' Decompose the fit into its additive, covariate, and latent parts
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return For `dynamic_ame`, a list of `T` decompositions, each with
#'   `additive`, `covariate`, `latent`, and `fitted` matrices that sum
#'   accordingly.
#' @export
decompose_fit <- function(object, ...) UseMethod("decompose_fit")

#' @rdname decompose_fit
#' @export
decompose_fit.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  Tn <- length(object$years)
  N <- nrow(object$a)
  M <- nrow(object$b)
  cv <- object$covariates

  out <- lapply(seq_len(Tn), function(t) {
    additive <- matrix(object$a[, t], N, M, byrow = FALSE) +
      matrix(object$b[, t], N, M, byrow = TRUE)
    covariate <- .cov_term(
      beta_t = if (nrow(object$beta) > 0L) object$beta[, t] else numeric(0),
      X_row  = if (is.null(cv$X_row_list))  NULL else cv$X_row_list[[t]],
      X_col  = if (is.null(cv$X_col_list))  NULL else cv$X_col_list[[t]],
      X_dyad = if (is.null(cv$X_dyad_list)) NULL else cv$X_dyad_list[[t]],
      N = N, M = M
    )
    latent <- object$U[[t]] %*% t(object$V[[t]])
    nm <- list(object$row_ids, object$col_ids)
    dimnames(additive) <- dimnames(covariate) <- dimnames(latent) <- nm
    fitted <- additive + covariate + latent
    list(additive = additive, covariate = covariate, latent = latent,
         fitted = fitted)
  })
  names(out) <- as.character(object$years)
  out
}


#' Convergence diagnostics from a dynamic AME fit
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return For `dynamic_ame`, a list with the convergence flag, iteration
#'   counts, the objective trace, the variance-component trace, the scale
#'   factors applied at each identification step, and the multistart objectives.
#' @export
convergence_info <- function(object, ...) UseMethod("convergence_info")

#' @rdname convergence_info
#' @export
convergence_info.dynamic_ame <- function(object, ...) {
  .validate_dynamic_ame(object)
  object$convergence
}
