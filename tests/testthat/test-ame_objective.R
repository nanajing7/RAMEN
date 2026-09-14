# The objective and the covariate helpers. Both are checked against
# element-by-element loops written straight from the model equation, so an error
# in the vectorised code cannot pass.

test_that("the covariate term matches a cell-by-cell computation", {
  d <- sim_ame_panel(seed = 5)
  p <- random_params(d)
  t <- 2L

  ct <- .cov_term(p$beta[, t], d$X_row[[t]], d$X_col[[t]], d$X_dyad[[t]],
                  d$N, d$M)
  brute <- matrix(0, d$N, d$M)
  for (i in seq_len(d$N)) {
    for (j in seq_len(d$M)) {
      x <- c(d$X_row[[t]][i, ], d$X_col[[t]][j, ], d$X_dyad[[t]][i, j, ])
      brute[i, j] <- sum(x * p$beta[, t])
    }
  }
  expect_equal(ct, brute, tolerance = 1e-12)
})


test_that("the design matrix uses the same stacking as the covariate term", {
  d <- sim_ame_panel(seed = 6)
  p <- random_params(d)
  t <- 1L
  Z <- .cov_design(d$X_row[[t]], d$X_col[[t]], d$X_dyad[[t]], d$N, d$M)
  expect_equal(dim(Z), c(d$N * d$M, d$P))
  expect_equal(as.vector(Z %*% p$beta[, t]),
               as.vector(.cov_term(p$beta[, t], d$X_row[[t]], d$X_col[[t]],
                                   d$X_dyad[[t]], d$N, d$M)),
               tolerance = 1e-12)
  expect_equal(.cov_dims(d$X_row[[t]], d$X_col[[t]], d$X_dyad[[t]]),
               list(Pr = 1L, Pc = 1L, Pd = 1L, P = 3L))
})


test_that("missing covariates propagate consistently in both helpers", {
  N <- 4L; M <- 5L
  X_row <- matrix(c(0.5, 1.5, NA, 2.0), N, 1)
  X_dyad <- array(stats::rnorm(N * M), c(N, M, 1))
  X_dyad[2, 3, 1] <- NA
  beta <- c(2, -1)

  ct <- .cov_term(beta, X_row, NULL, X_dyad, N, M)
  zt <- matrix(as.vector(.cov_design(X_row, NULL, X_dyad, N, M) %*% beta), N, M)

  expect_identical(is.na(ct), is.na(zt))
  expect_equal(ct[!is.na(ct)], zt[!is.na(zt)], tolerance = 1e-12)

  # a missing node covariate knocks out that node's whole row; a missing dyadic
  # covariate knocks out only its own cell
  expect_true(all(is.na(ct[3, ])))
  expect_true(is.na(ct[2, 3]))
  expect_true(all(!is.na(ct[2, -3])))
  expect_equal(sum(is.na(ct)), M + 1L)
})


test_that("Q matches an independent element-by-element computation", {
  d <- sim_ame_panel(N = 9, M = 6, Tn = 4, Pr = 2, Pc = 1, Pd = 2, seed = 7)
  p <- random_params(d)
  lam <- c(beta = .3, a = .4, b = .5, U = .6, V = .7)
  gam <- c(beta = 1.1, a = 1.2, b = 1.3, U = 1.4, V = 1.5)

  ssr <- 0
  nobs <- 0L
  for (t in seq_len(d$Tn)) {
    for (i in seq_len(d$N)) {
      for (j in seq_len(d$M)) {
        y <- d$Y[[t]][i, j]
        if (is.na(y)) next
        x <- c(d$X_row[[t]][i, ], d$X_col[[t]][j, ], d$X_dyad[[t]][i, j, ])
        f <- p$a[i, t] + p$b[j, t] + sum(x * p$beta[, t]) +
          sum(p$U[[t]][i, ] * p$V[[t]][j, ])
        ssr <- ssr + (y - f)^2
        nobs <- nobs + 1L
      }
    }
  }
  dsum <- function(f) sum(vapply(2:d$Tn, function(t) sum((f(t) - f(t - 1))^2),
                                 numeric(1)))
  pi_b <- c(beta = lam[["beta"]] * sum(p$beta[, 1]^2),
            a = lam[["a"]] * sum(p$a[, 1]^2),
            b = lam[["b"]] * sum(p$b[, 1]^2),
            U = lam[["U"]] * sum(p$U[[1]]^2),
            V = lam[["V"]] * sum(p$V[[1]]^2))
  pr_b <- c(beta = gam[["beta"]] * dsum(function(t) p$beta[, t]),
            a = gam[["a"]] * dsum(function(t) p$a[, t]),
            b = gam[["b"]] * dsum(function(t) p$b[, t]),
            U = gam[["U"]] * dsum(function(t) p$U[[t]]),
            V = gam[["V"]] * dsum(function(t) p$V[[t]]))

  res <- .ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, lam, gam)
  expect_equal(res$ssr, ssr, tolerance = 1e-9)
  expect_equal(res$n_obs, nobs)
  expect_equal(unname(res$penalty_initial), unname(pi_b), tolerance = 1e-9)
  expect_equal(unname(res$penalty_rw), unname(pr_b), tolerance = 1e-9)
  expect_equal(res$Q, ssr + sum(pi_b) + sum(pr_b), tolerance = 1e-9)
})


test_that("Q handles the degenerate configurations the spec allows", {
  d <- sim_ame_panel(seed = 8)
  p <- random_params(d)

  # no covariates
  p0 <- p
  p0$beta <- matrix(0, 0, d$Tn)
  r0 <- .ame_objective(p0, d$Y, NULL, NULL, NULL, 1, 1)
  expect_true(is.finite(r0$Q))
  expect_equal(unname(r0$penalty_initial[["beta"]]), 0)
  expect_equal(unname(r0$penalty_rw[["beta"]]), 0)

  # a single period leaves no random-walk term
  p1 <- list(a = p$a[, 1, drop = FALSE], b = p$b[, 1, drop = FALSE],
             beta = p$beta[, 1, drop = FALSE], U = p$U[1], V = p$V[1])
  r1 <- .ame_objective(p1, d$Y[1], d$X_row[1], d$X_col[1], d$X_dyad[1], 1, 1)
  expect_true(all(r1$penalty_rw == 0))

  # zero penalties reduce Q to the residual sum of squares
  rz <- .ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, 0, 0)
  expect_equal(rz$Q, rz$ssr, tolerance = 1e-12)

  # a constant trajectory has no innovation penalty
  pc <- p
  for (t in 2:d$Tn) {
    pc$a[, t] <- pc$a[, 1]; pc$b[, t] <- pc$b[, 1]
    pc$beta[, t] <- pc$beta[, 1]
    pc$U[[t]] <- pc$U[[1]]; pc$V[[t]] <- pc$V[[1]]
  }
  expect_equal(sum(.ame_objective(pc, d$Y, d$X_row, d$X_col, d$X_dyad, 1, 1)$penalty_rw),
               0, tolerance = 1e-12)
})


test_that("penalty arguments are validated and recycled", {
  d <- sim_ame_panel(seed = 9)
  p <- random_params(d)
  lam <- c(beta = .5, a = .5, b = .5, U = .5, V = .5)

  expect_equal(.ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, 0.5, 2)$Q,
               .ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, lam,
                              lam * 4)$Q,
               tolerance = 1e-12)
  # name order is irrelevant
  expect_equal(.ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                              lam[c(5, 4, 3, 2, 1)], lam)$Q,
               .ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, lam, lam)$Q,
               tolerance = 1e-12)
  expect_error(.ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, lam[1:3], lam),
               "missing entries")
  expect_error(.ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, -lam, lam),
               "non-negative")
})


test_that("fitted and residual helpers agree with each other", {
  d <- sim_ame_panel(seed = 10)
  p <- random_params(d)
  f <- .fitted_period(p$a[, 1], p$b[, 1], p$U[[1]], p$V[[1]], p$beta[, 1],
                      d$X_row[[1]], d$X_col[[1]], d$X_dyad[[1]])
  r <- .resid_period(d$Y[[1]], p$a[, 1], p$b[, 1], p$U[[1]], p$V[[1]],
                     p$beta[, 1], d$X_row[[1]], d$X_col[[1]], d$X_dyad[[1]])
  expect_equal(r, d$Y[[1]] - f, tolerance = 1e-12)
  expect_identical(is.na(r), is.na(d$Y[[1]]))
})
