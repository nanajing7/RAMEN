
#' Fit one-period ALS with row/column/dyad covariates and temporal smoothing
#'
#' Model:
#'   Y_ij = alpha_i + beta_j + X_row_i * coef_row + X_col_j * coef_col
#'        + X_dyad_ij * coef_dyad + U_i'V_j + error
#'
#' @param Y Outcome matrix.
#' @param K Latent dimension.
#' @param lambda Ridge penalty on U and V.
#' @param gamma Temporal smoothing penalty.
#' @param U_prev,V_prev Previous-period latent factors.
#' @param alpha_prev,beta_prev Previous-period intercepts.
#' @param U_init,V_init Initial values for latent factors.
#' @param X_row Row covariate matrix.
#' @param X_col Column covariate matrix.
#' @param X_dyad Dyad covariate array.
#' @param coef_row_prev,coef_col_prev,coef_dyad_prev Previous-period covariate coefficients.
#' @param coef_row_init,coef_col_init,coef_dyad_init Initial covariate coefficients.
#' @param max_iter Maximum iterations.
#' @param tol Relative tolerance.
#'
#' @return A list containing parameter estimates and diagnostics.
als_factorize_joint_cov <- function(Y,
                                    K = 2,
                                    lambda = 0.1,
                                    gamma = 1,
                                    U_prev = NULL,
                                    V_prev = NULL,
                                    alpha_prev = NULL,
                                    beta_prev = NULL,
                                    U_init = NULL,
                                    V_init = NULL,
                                    X_row = NULL,
                                    X_col = NULL,
                                    X_dyad = NULL,
                                    coef_row_prev = NULL,
                                    coef_col_prev = NULL,
                                    coef_dyad_prev = NULL,
                                    coef_row_init = NULL,
                                    coef_col_init = NULL,
                                    coef_dyad_init = NULL,
                                    max_iter = 1000,
                                    tol = 1e-5) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  N <- nrow(Y)
  M <- ncol(Y)

  if (length(K) != 1 || K < 1 || K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (K > min(N, M)) {
    stop(sprintf("K = %d exceeds min(nrow(Y), ncol(Y)) = %d.", K, min(N, M)), call. = FALSE)
  }

  if (xor(is.null(U_prev), is.null(V_prev))) {
    stop("U_prev and V_prev must both be NULL or both be provided.", call. = FALSE)
  }
  if (xor(is.null(U_init), is.null(V_init))) {
    stop("U_init and V_init must both be NULL or both be provided.", call. = FALSE)
  }

  validate_matrix <- function(X, n_expected, k_expected, name) {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (!all(dim(X) == c(n_expected, k_expected))) {
      stop(sprintf("%s must have dimension %d x %d.", name, n_expected, k_expected), call. = FALSE)
    }
    storage.mode(X) <- "double"
    X
  }

  validate_vector <- function(x, n_expected, name) {
    x <- as.numeric(x)
    if (length(x) != n_expected) {
      stop(sprintf("%s must have length %d.", name, n_expected), call. = FALSE)
    }
    x
  }

  validate_array3 <- function(X, d1, d2, name) {
    if (length(dim(X)) != 3) {
      stop(sprintf("%s must be a 3D array.", name), call. = FALSE)
    }
    if (!all(dim(X)[1:2] == c(d1, d2))) {
      stop(sprintf("%s first two dimensions must be %d x %d.", name, d1, d2), call. = FALSE)
    }
    storage.mode(X) <- "double"
    X
  }

  if (!is.null(U_prev)) U_prev <- validate_matrix(U_prev, N, K, "U_prev")
  if (!is.null(V_prev)) V_prev <- validate_matrix(V_prev, M, K, "V_prev")
  if (!is.null(U_init)) U_init <- validate_matrix(U_init, N, K, "U_init")
  if (!is.null(V_init)) V_init <- validate_matrix(V_init, M, K, "V_init")
  if (!is.null(alpha_prev)) alpha_prev <- validate_vector(alpha_prev, N, "alpha_prev")
  if (!is.null(beta_prev)) beta_prev <- validate_vector(beta_prev, M, "beta_prev")

  Pr <- 0
  Pc <- 0
  Pd <- 0

  if (!is.null(X_row)) {
    X_row <- as.matrix(X_row)
    if (nrow(X_row) != N) stop("X_row must have N rows.", call. = FALSE)
    storage.mode(X_row) <- "double"
    X_row[is.na(X_row)] <- 0
    Pr <- ncol(X_row)
  }

  if (!is.null(X_col)) {
    X_col <- as.matrix(X_col)
    if (nrow(X_col) != M) stop("X_col must have M rows.", call. = FALSE)
    storage.mode(X_col) <- "double"
    X_col[is.na(X_col)] <- 0
    Pc <- ncol(X_col)
  }

  if (!is.null(X_dyad)) {
    X_dyad <- validate_array3(X_dyad, N, M, "X_dyad")
    X_dyad[is.na(X_dyad)] <- 0
    Pd <- dim(X_dyad)[3]
  }

  if (!is.null(coef_row_prev))  coef_row_prev  <- validate_vector(coef_row_prev,  Pr, "coef_row_prev")
  if (!is.null(coef_col_prev))  coef_col_prev  <- validate_vector(coef_col_prev,  Pc, "coef_col_prev")
  if (!is.null(coef_dyad_prev)) coef_dyad_prev <- validate_vector(coef_dyad_prev, Pd, "coef_dyad_prev")

  if (!is.null(U_init) && !is.null(V_init)) {
    U <- U_init
    V <- V_init
  } else if (!is.null(U_prev) && !is.null(V_prev)) {
    U <- U_prev
    V <- V_prev
  } else {
    init <- svd_init(Y, K)
    U <- init$U
    V <- init$V
  }

  alpha <- if (!is.null(alpha_prev)) alpha_prev else rowMeans(Y)
  beta  <- if (!is.null(beta_prev)) beta_prev else colMeans(Y) - mean(Y)

  coef_row <- if (Pr > 0) {
    if (!is.null(coef_row_init)) as.numeric(coef_row_init)
    else if (!is.null(coef_row_prev)) as.numeric(coef_row_prev)
    else rep(0, Pr)
  } else numeric(0)

  coef_col <- if (Pc > 0) {
    if (!is.null(coef_col_init)) as.numeric(coef_col_init)
    else if (!is.null(coef_col_prev)) as.numeric(coef_col_prev)
    else rep(0, Pc)
  } else numeric(0)

  coef_dyad <- if (Pd > 0) {
    if (!is.null(coef_dyad_init)) as.numeric(coef_dyad_init)
    else if (!is.null(coef_dyad_prev)) as.numeric(coef_dyad_prev)
    else rep(0, Pd)
  } else numeric(0)

  if (length(coef_row) != Pr) stop("coef_row_init has wrong length.", call. = FALSE)
  if (length(coef_col) != Pc) stop("coef_col_init has wrong length.", call. = FALSE)
  if (length(coef_dyad) != Pd) stop("coef_dyad_init has wrong length.", call. = FALSE)

  row_ids <- rownames(Y)
  col_ids <- colnames(Y)
  latent_ids <- paste0("k", seq_len(K))

  rownames(U) <- row_ids
  rownames(V) <- col_ids
  colnames(U) <- latent_ids
  colnames(V) <- latent_ids
  names(alpha) <- row_ids
  names(beta) <- col_ids

  gamma_U <- if (!is.null(U_prev)) gamma else 0
  gamma_V <- if (!is.null(V_prev)) gamma else 0
  gamma_alpha <- if (!is.null(alpha_prev)) gamma else 0
  gamma_beta  <- if (!is.null(beta_prev)) gamma else 0

  gamma_coef_row  <- if (!is.null(coef_row_prev)  && Pr > 0) gamma else 0
  gamma_coef_col  <- if (!is.null(coef_col_prev)  && Pc > 0) gamma else 0
  gamma_coef_dyad <- if (!is.null(coef_dyad_prev) && Pd > 0) gamma else 0

  obj_old <- NA_real_
  obj <- NA_real_

  build_cov_term <- function() {
    out <- matrix(0, N, M)

    if (Pr > 0) {
      row_term <- as.vector(X_row %*% coef_row)
      out <- out + matrix(row_term, N, M, byrow = FALSE)
    }

    if (Pc > 0) {
      col_term <- as.vector(X_col %*% coef_col)
      out <- out + matrix(col_term, N, M, byrow = TRUE)
    }

    if (Pd > 0) {
      dyad_term <- matrix(0, N, M)
      for (p in seq_len(Pd)) {
        dyad_term <- dyad_term + coef_dyad[p] * X_dyad[, , p]
      }
      out <- out + dyad_term
    }

    out
  }

  for (iter in seq_len(max_iter)) {
    UVt <- U %*% t(V)

    R_cov <- Y -
      matrix(alpha, N, M, byrow = FALSE) -
      matrix(beta, N, M, byrow = TRUE) -
      UVt

    Z_list <- list()
    prior_mean <- numeric(0)
    prior_weight <- numeric(0)

    if (Pr > 0) {
      Z_row <- X_row[rep(seq_len(N), times = M), , drop = FALSE]
      Z_list <- c(Z_list, list(Z_row))
      if (gamma_coef_row > 0) {
        prior_mean <- c(prior_mean, coef_row_prev)
        prior_weight <- c(prior_weight, rep(gamma_coef_row, Pr))
      } else {
        prior_mean <- c(prior_mean, rep(0, Pr))
        prior_weight <- c(prior_weight, rep(0, Pr))
      }
    }

    if (Pc > 0) {
      Z_col <- X_col[rep(seq_len(M), each = N), , drop = FALSE]
      Z_list <- c(Z_list, list(Z_col))
      if (gamma_coef_col > 0) {
        prior_mean <- c(prior_mean, coef_col_prev)
        prior_weight <- c(prior_weight, rep(gamma_coef_col, Pc))
      } else {
        prior_mean <- c(prior_mean, rep(0, Pc))
        prior_weight <- c(prior_weight, rep(0, Pc))
      }
    }

    if (Pd > 0) {
      Z_dyad <- sapply(seq_len(Pd), function(p) as.vector(X_dyad[, , p]))
      if (Pd == 1) Z_dyad <- matrix(Z_dyad, ncol = 1)
      Z_list <- c(Z_list, list(Z_dyad))
      if (gamma_coef_dyad > 0) {
        prior_mean <- c(prior_mean, coef_dyad_prev)
        prior_weight <- c(prior_weight, rep(gamma_coef_dyad, Pd))
      } else {
        prior_mean <- c(prior_mean, rep(0, Pd))
        prior_weight <- c(prior_weight, rep(0, Pd))
      }
    }

    if (length(Z_list) > 0) {
      Z <- do.call(cbind, Z_list)
      y_vec <- as.vector(R_cov)

      XtX <- crossprod(Z)
      Xty <- crossprod(Z, y_vec)

      if (any(prior_weight > 0)) {
        XtX <- XtX + diag(prior_weight, nrow = length(prior_weight))
        Xty <- Xty + prior_weight * prior_mean
      }

      coef_all <- qr.solve(XtX, Xty)

      idx <- 1
      if (Pr > 0) {
        coef_row <- coef_all[idx:(idx + Pr - 1)]
        idx <- idx + Pr
      }
      if (Pc > 0) {
        coef_col <- coef_all[idx:(idx + Pc - 1)]
        idx <- idx + Pc
      }
      if (Pd > 0) {
        coef_dyad <- coef_all[idx:(idx + Pd - 1)]
      }
    }

    cov_term <- build_cov_term()

    alpha_num <- rowSums(Y - cov_term - matrix(beta, N, M, byrow = TRUE) - UVt)
    if (gamma_alpha > 0) alpha_num <- alpha_num + gamma * alpha_prev
    alpha <- alpha_num / (M + gamma_alpha)

    beta_num <- colSums(Y - cov_term - matrix(alpha, N, M, byrow = FALSE) - UVt)
    if (gamma_beta > 0) beta_num <- beta_num + gamma * beta_prev
    beta <- beta_num / (N + gamma_beta)

    R <- Y - cov_term -
      matrix(alpha, N, M, byrow = FALSE) -
      matrix(beta, N, M, byrow = TRUE)

    A_U <- crossprod(V) + (lambda + gamma_U) * diag(K)
    rhs_U <- R %*% V
    if (gamma_U > 0) rhs_U <- rhs_U + gamma * U_prev
    U <- rhs_U %*% solve(A_U)

    A_V <- crossprod(U) + (lambda + gamma_V) * diag(K)
    rhs_V <- t(R) %*% U
    if (gamma_V > 0) rhs_V <- rhs_V + gamma * V_prev
    V <- rhs_V %*% solve(A_V)

    UVt <- U %*% t(V)
    cov_term <- build_cov_term()

    resid <- Y - cov_term -
      matrix(alpha, N, M, byrow = FALSE) -
      matrix(beta, N, M, byrow = TRUE) -
      UVt

    fit_term <- sum(resid^2)
    ridge_term <- lambda * (sum(U^2) + sum(V^2))

    temporal_term <- 0
    if (gamma_U > 0) temporal_term <- temporal_term + gamma * sum((U - U_prev)^2)
    if (gamma_V > 0) temporal_term <- temporal_term + gamma * sum((V - V_prev)^2)
    if (gamma_alpha > 0) temporal_term <- temporal_term + gamma * sum((alpha - alpha_prev)^2)
    if (gamma_beta > 0) temporal_term <- temporal_term + gamma * sum((beta - beta_prev)^2)
    if (gamma_coef_row > 0) temporal_term <- temporal_term + gamma_coef_row * sum((coef_row - coef_row_prev)^2)
    if (gamma_coef_col > 0) temporal_term <- temporal_term + gamma_coef_col * sum((coef_col - coef_col_prev)^2)
    if (gamma_coef_dyad > 0) temporal_term <- temporal_term + gamma_coef_dyad * sum((coef_dyad - coef_dyad_prev)^2)

    obj <- fit_term + ridge_term + temporal_term

    if (!is.na(obj_old)) {
      rel_change <- abs(obj_old - obj) / (abs(obj_old) + 1e-8)
      if (rel_change < tol) break
    }
    obj_old <- obj
  }

  rownames(U) <- row_ids
  rownames(V) <- col_ids
  colnames(U) <- latent_ids
  colnames(V) <- latent_ids
  names(alpha) <- row_ids
  names(beta) <- col_ids

  list(
    U = U,
    V = V,
    alpha = alpha,
    beta = beta,
    coef_row = coef_row,
    coef_col = coef_col,
    coef_dyad = coef_dyad,
    objective = obj,
    iterations = iter
  )
}

