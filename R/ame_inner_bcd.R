# The inner block coordinate descent loop: Steps 2-6 swept repeatedly at fixed
# penalties until the penalized MAP objective stabilizes.

#' Run the inner block coordinate descent at fixed penalties
#'
#' Sweeps the five trajectory blocks in the order given by the spec — `beta`
#' (Step 2), `a` (Step 3), `b` (Step 4), `U` (Step 5), `V` (Step 6) — and
#' repeats until the relative change in the objective falls below `eps_Q`:
#' \deqn{\frac{|Q^{(k)} - Q^{(k-1)}|}{|Q^{(k-1)}| + \delta} < \epsilon_Q.}
#'
#' Every block update is the exact minimiser of `Q` over that block with the
#' others held fixed, so `Q` is non-increasing along the whole sweep and the
#' sequence converges. `U` and `V` are updated alternately rather than jointly
#' because the multiplicative term is jointly non-convex in the pair but
#' conditionally quadratic in either one.
#'
#' The penalties are held fixed here; they are refreshed by the empirical-Bayes
#' step between outer iterations.
#'
#' @param params Starting parameter list (`a`, `b`, `beta`, `U`, `V`).
#' @param Y_list List of `T` outcome matrices, `NA` where unobserved.
#' @param X_row_list,X_col_list,X_dyad_list Covariate lists; may be `NULL`.
#' @param lambda,gamma Named penalty vectors over the blocks `beta`, `a`, `b`,
#'   `U`, `V`; a single scalar is recycled across all five.
#' @param max_iter Maximum number of sweeps.
#' @param eps_Q Relative-change tolerance on the objective.
#' @param delta Small constant guarding the denominator.
#' @param verbose If `TRUE`, report the objective at each sweep.
#'
#' @return A list with the updated `params`, the final objective breakdown
#'   `objective`, the vector of objective values `Q_trace` (one entry per sweep,
#'   preceded by the starting value), the number of sweeps `iterations`, and
#'   `converged`.
#' @keywords internal
#' @noRd
.ame_inner_bcd <- function(params,
                           Y_list,
                           X_row_list = NULL,
                           X_col_list = NULL,
                           X_dyad_list = NULL,
                           lambda,
                           gamma,
                           max_iter = 500,
                           eps_Q = 1e-6,
                           delta = 1e-8,
                           verbose = FALSE) {
  blocks <- c("beta", "a", "b", "U", "V")
  lambda <- .check_penalty(lambda, "lambda", blocks)
  gamma  <- .check_penalty(gamma,  "gamma",  blocks)

  obj <- .ame_objective(params, Y_list, X_row_list, X_col_list, X_dyad_list,
                        lambda, gamma)
  Q_trace <- obj$Q
  converged <- FALSE
  iter <- 0L

  for (k in seq_len(max_iter)) {
    iter <- k
    Q_prev <- obj$Q

    # Spec order: Step 2 -> 3 -> 4 -> 5 -> 6.
    params <- .update_beta_block(params, Y_list, X_row_list, X_col_list,
                                 X_dyad_list, lambda["beta"], gamma["beta"])
    params <- .update_a_block(params, Y_list, X_row_list, X_col_list,
                              X_dyad_list, lambda["a"], gamma["a"])
    params <- .update_b_block(params, Y_list, X_row_list, X_col_list,
                              X_dyad_list, lambda["b"], gamma["b"])
    params <- .update_U_block(params, Y_list, X_row_list, X_col_list,
                              X_dyad_list, lambda["U"], gamma["U"])
    params <- .update_V_block(params, Y_list, X_row_list, X_col_list,
                              X_dyad_list, lambda["V"], gamma["V"])

    obj <- .ame_objective(params, Y_list, X_row_list, X_col_list, X_dyad_list,
                          lambda, gamma)
    Q_trace <- c(Q_trace, obj$Q)

    if (verbose) {
      cat(sprintf("    inner %3d  Q = %.8f\n", k, obj$Q))
    }

    if (abs(Q_prev - obj$Q) / (abs(Q_prev) + delta) < eps_Q) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      sprintf("Inner BCD did not converge in %d sweeps (last relative change %.3g).",
              max_iter,
              abs(Q_trace[length(Q_trace) - 1] - obj$Q) /
                (abs(Q_trace[length(Q_trace) - 1]) + delta)),
      call. = FALSE
    )
  }

  list(
    params = params,
    objective = obj,
    Q_trace = Q_trace,
    iterations = iter,
    converged = converged
  )
}
