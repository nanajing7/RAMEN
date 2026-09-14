# The MAP objective minimised by the inner block coordinate descent of
# fit_dynamic_ame, at fixed penalties.

#' Penalised least-squares objective for a parameter trajectory
#'
#' Computes
#' \deqn{Q = \sum_t \sum_{(i,j) \in O_t} r_{ijt}^2
#'         + \sum_\bullet \lambda_\bullet \|\theta_{\bullet,1}\|^2
#'         + \sum_\bullet \gamma_\bullet \sum_{t=2}^T
#'           \|\theta_{\bullet,t} - \theta_{\bullet,t-1}\|^2 }
#' where `r` is the residual over observed cells and the five blocks are
#' `beta`, `a`, `b`, `U`, `V`. The first penalty group comes from the Gaussian
#' priors on the initial states, the second from the Gaussian random-walk
#' innovations.
#'
#' This is exactly the quantity the block updates minimise, so each block
#' update must leave `Q` non-increasing — the property the monotonicity tests
#' check. The variance components themselves are estimated separately in the
#' empirical-Bayes step; here `lambda` and `gamma` are held fixed.
#'
#' @param params List with `a` (`N x T`), `b` (`M x T`), `beta` (`P x T`), and
#'   `U`, `V` (lists of length `T` of `N x K` and `M x K` matrices).
#' @param Y_list List of length `T` of `N x M` outcome matrices, `NA` where
#'   unobserved.
#' @param X_row_list,X_col_list,X_dyad_list Lists of length `T` of covariate
#'   pieces; elements may be `NULL`.
#' @param lambda Named numeric vector of initial-state penalties, with names
#'   `beta`, `a`, `b`, `U`, `V`.
#' @param gamma Named numeric vector of random-walk penalties, same names.
#'
#' @return A list with the total `Q`, the residual sum of squares `ssr`, the
#'   named vectors `penalty_initial` and `penalty_rw`, and `n_obs` (the number
#'   of observed cells).
#' @keywords internal
#' @noRd
.ame_objective <- function(params,
                           Y_list,
                           X_row_list = NULL,
                           X_col_list = NULL,
                           X_dyad_list = NULL,
                           lambda,
                           gamma) {
  Tn <- length(Y_list)
  blocks <- c("beta", "a", "b", "U", "V")

  lambda <- .check_penalty(lambda, "lambda", blocks)
  gamma  <- .check_penalty(gamma,  "gamma",  blocks)

  # ── Residual sum of squares over observed cells ──────────────────────────
  ssr <- 0
  n_obs <- 0L

  for (t in seq_len(Tn)) {
    r <- .resid_period(
      Y_t    = Y_list[[t]],
      a_t    = params$a[, t],
      b_t    = params$b[, t],
      U_t    = params$U[[t]],
      V_t    = params$V[[t]],
      beta_t = if (nrow(params$beta) > 0L) params$beta[, t] else numeric(0),
      X_row  = if (is.null(X_row_list))  NULL else X_row_list[[t]],
      X_col  = if (is.null(X_col_list))  NULL else X_col_list[[t]],
      X_dyad = if (is.null(X_dyad_list)) NULL else X_dyad_list[[t]]
    )
    ssr <- ssr + sum(r^2, na.rm = TRUE)
    n_obs <- n_obs + sum(!is.na(r))
  }

  # ── Trajectories as a uniform list of per-period matrices ────────────────
  # Columns of a/b/beta and elements of U/V are all treated the same way:
  # squared Frobenius norm at t = 1, squared successive differences after.
  traj <- list(
    beta = .cols_to_list(params$beta),
    a    = .cols_to_list(params$a),
    b    = .cols_to_list(params$b),
    U    = params$U,
    V    = params$V
  )

  pen_init <- stats::setNames(numeric(length(blocks)), blocks)
  pen_rw   <- stats::setNames(numeric(length(blocks)), blocks)

  for (nm in blocks) {
    th <- traj[[nm]]

    pen_init[nm] <- lambda[nm] * sum(th[[1]]^2)

    if (Tn > 1L) {
      d <- 0
      for (t in 2:Tn) d <- d + sum((th[[t]] - th[[t - 1]])^2)
      pen_rw[nm] <- gamma[nm] * d
    }
  }

  list(
    Q               = ssr + sum(pen_init) + sum(pen_rw),
    ssr             = ssr,
    penalty_initial = pen_init,
    penalty_rw      = pen_rw,
    n_obs           = n_obs
  )
}


#' Split a `q x T` matrix into a list of `T` length-`q` vectors
#'
#' Lets `a`, `b`, and `beta` be handled with the same code path as the `U` and
#' `V` lists. A zero-row matrix yields zero-length vectors, so blocks with no
#' covariates contribute nothing.
#'
#' @param M A `q x T` matrix.
#' @return A list of length `T`.
#' @keywords internal
#' @noRd
.cols_to_list <- function(M) {
  lapply(seq_len(ncol(M)), function(t) M[, t])
}


#' Validate and normalise a penalty vector
#'
#' Accepts either a named vector covering all five blocks or a single scalar to
#' be recycled across them.
#'
#' @param x Penalty input.
#' @param nm Name used in error messages.
#' @param blocks Required block names.
#'
#' @return A named numeric vector indexed by `blocks`.
#' @keywords internal
#' @noRd
.check_penalty <- function(x, nm, blocks) {
  if (length(x) == 1L && is.null(names(x))) {
    x <- stats::setNames(rep(as.numeric(x), length(blocks)), blocks)
  }

  if (is.null(names(x))) {
    stop(sprintf("%s must be a named vector or a single scalar.", nm), call. = FALSE)
  }

  missing_nm <- setdiff(blocks, names(x))
  if (length(missing_nm) > 0L) {
    stop(
      sprintf("%s is missing entries for: %s", nm, paste(missing_nm, collapse = ", ")),
      call. = FALSE
    )
  }

  x <- x[blocks]
  if (any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("%s must be finite and non-negative.", nm), call. = FALSE)
  }

  x
}
