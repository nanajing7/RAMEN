# Step 9, the EM M-step for the variance components.
#
# The plug-in version of this step (sums of squares of the point estimates, with
# no posterior-variance term) drives the innovation variances to zero: the
# estimates are already smoothed by the penalty, so their differences understate
# the innovation, which raises the penalty, which smooths them further. The
# regression test below pins that behaviour down.

test_that("the point-estimate sums of squares follow the spec's divisors", {
  # With the posterior term switched off by an overwhelming signal, the M-step
  # must reduce to the spec's plug-in formulas.
  set.seed(5)
  N <- 40L; M <- 30L; K <- 2L; Tn <- 5L
  s <- list(a = matrix(0, N, Tn), b = matrix(0, M, Tn), beta = matrix(0, 0, Tn),
            U = vector("list", Tn), V = vector("list", Tn))
  s$a[, 1] <- stats::rnorm(N); s$b[, 1] <- stats::rnorm(M)
  s$U[[1]] <- matrix(stats::rnorm(N * K), N, K)
  s$V[[1]] <- matrix(stats::rnorm(M * K), M, K)
  for (t in 2:Tn) {
    s$a[, t] <- s$a[, t - 1] + stats::rnorm(N, sd = .5)
    s$b[, t] <- s$b[, t - 1] + stats::rnorm(M, sd = .5)
    s$U[[t]] <- s$U[[t - 1]] + matrix(stats::rnorm(N * K, sd = .5), N, K)
    s$V[[t]] <- s$V[[t - 1]] + matrix(stats::rnorm(M * K, sd = .5), M, K)
  }
  Y <- lapply(seq_len(Tn), function(t)
    .fitted_period(s$a[, t], s$b[, t], s$U[[t]], s$V[[t]]) +
      matrix(stats::rnorm(N * M, sd = 0.02), N, M))

  eb <- .eb_update_variances(s, Y, NULL, NULL, NULL,
                             lambda = unit_penalty() * 1e-4,
                             gamma = unit_penalty() * 1e-4)

  # a tight posterior leaves the trace term negligible
  expect_lt(eb$trace_share$innovation[["a"]], 0.02)

  plug_tau_a <- sum(vapply(2:Tn, function(t) sum((s$a[, t] - s$a[, t - 1])^2),
                           numeric(1))) / (N * (Tn - 1))
  plug_sig_a <- sum(s$a[, 1]^2) / N
  plug_tau_U <- sum(vapply(2:Tn, function(t) sum((s$U[[t]] - s$U[[t - 1]])^2),
                           numeric(1))) / (N * K * (Tn - 1))
  expect_equal(eb$Omega[["tau_a2"]], plug_tau_a, tolerance = 0.02)
  expect_equal(eb$Omega[["sigma_a2"]], plug_sig_a, tolerance = 0.02)

  # The latent innovation variances are tied by default, so what is reported
  # for `tau_U2` is the geometric mean of the two plug-in values rather than
  # either one. The formula itself is checked with the tie switched off; that
  # the tie then holds is checked separately.
  eb_sep <- .eb_update_variances(s, Y, NULL, NULL, NULL,
                                 lambda = unit_penalty() * 1e-4,
                                 gamma = unit_penalty() * 1e-4,
                                 tie_latent = FALSE)
  expect_equal(eb_sep$Omega[["tau_U2"]], plug_tau_U, tolerance = 0.02)
  expect_equal(eb$Omega[["tau_U2"]],
               sqrt(eb_sep$Omega[["tau_U2"]] * eb_sep$Omega[["tau_V2"]]),
               tolerance = 1e-10)
})


test_that("sigma_eps2 is the mean squared residual over observed cells", {
  d <- sim_ame_panel(seed = 61)
  p <- random_params(d)
  eb <- .eb_update_variances(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                             lambda = unit_penalty(), gamma = unit_penalty())
  ssr <- 0; nobs <- 0L
  for (t in seq_len(d$Tn)) {
    r <- .resid_period(d$Y[[t]], p$a[, t], p$b[, t], p$U[[t]], p$V[[t]],
                       p$beta[, t], d$X_row[[t]], d$X_col[[t]], d$X_dyad[[t]])
    ssr <- ssr + sum(r^2, na.rm = TRUE)
    nobs <- nobs + sum(!is.na(r))
  }
  expect_equal(eb$Omega[["sigma_eps2"]], ssr / nobs, tolerance = 1e-12)
  expect_equal(eb$n_obs, nobs)
})


test_that("Omega carries the spec's eleven entries in order", {
  d <- sim_ame_panel(seed = 62)
  eb <- .eb_update_variances(random_params(d), d$Y, d$X_row, d$X_col, d$X_dyad,
                             lambda = unit_penalty(), gamma = unit_penalty())
  expect_named(eb$Omega,
               c("sigma_eps2", "sigma_beta2", "sigma_a2", "sigma_b2",
                 "sigma_U2", "sigma_V2", "tau_beta2", "tau_a2", "tau_b2",
                 "tau_U2", "tau_V2"))
  expect_equal(unname(eb$lambda),
               unname(eb$Omega[["sigma_eps2"]] /
                        eb$Omega[c("sigma_beta2", "sigma_a2", "sigma_b2",
                                   "sigma_U2", "sigma_V2")]),
               tolerance = 1e-12)
  expect_equal(unname(eb$gamma),
               unname(eb$Omega[["sigma_eps2"]] /
                        eb$Omega[c("tau_beta2", "tau_a2", "tau_b2",
                                   "tau_U2", "tau_V2")]),
               tolerance = 1e-12)
})


test_that("the posterior correction keeps the innovation variances alive", {
  # Regression test for the variance collapse. On this panel the plug-in update
  # drove tau_U2 from ~9e-2 down to the eps_var floor of 1e-8 over 25 outer
  # iterations, taking gamma_U to 4e7 and destroying the fit.
  d <- sim_ame_panel(N = 20, M = 14, Tn = 8, Pr = 0, Pc = 0, Pd = 0,
                     sd_eps = 0.55, na_frac = 0.12, seed = 1234)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = 0, n_starts = 1)

  p <- ini$candidates[[1]]
  lam <- unit_penalty(); gam <- unit_penalty()
  tau_U <- numeric(0)
  for (m in 1:12) {
    bc <- suppressWarnings(
      .ame_inner_bcd(p, d$Y, NULL, NULL, NULL, lam, gam,
                     eps_Q = 1e-6, max_iter = 300))
    p <- .ame_identify(bc$params, ini$reference)$params
    eb <- .eb_update_variances(p, d$Y, NULL, NULL, NULL,
                               lambda = lam, gamma = gam)
    tau_U <- c(tau_U, eb$Omega[["tau_U2"]])
    lam <- eb$lambda; gam <- eb$gamma
  }

  expect_gt(min(tau_U), 1e-3)             # nowhere near the 1e-8 floor
  expect_gt(tau_U[12], 0.3 * tau_U[1])    # no monotone collapse
  expect_lt(max(gam), 100)                # penalties stay in a sane range
  expect_false(any(eb$capped))

  # and the correction is a real, non-dominating share of the estimate
  sh <- eb$trace_share$innovation[c("a", "b", "U", "V")]
  expect_true(all(sh > 0.01))
  expect_true(all(sh < 0.95))
})


test_that("the latent innovation variances are tied, and only they are", {
  # Only the product `tau_U^2 tau_V^2` is determined by the data: rescaling U
  # against V moves the split without touching any `U_t V_t'`. The default
  # therefore reports the geometric mean for both and leaves the other
  # components alone.
  d <- sim_ame_panel(N = 14, M = 11, K = 2, Tn = 5, Pr = 0, Pc = 0, Pd = 1,
                     sd_eps = 0.4, na_frac = 0, seed = 808)
  p <- random_params(d, seed = 3)
  arg <- list(p, d$Y, NULL, NULL, d$X_dyad, lambda = unit_penalty(),
              gamma = unit_penalty())

  tied <- do.call(.eb_update_variances, arg)
  free <- do.call(.eb_update_variances, c(arg, list(tie_latent = FALSE)))

  expect_identical(tied$Omega[["tau_U2"]], tied$Omega[["tau_V2"]])
  expect_equal(tied$Omega[["tau_U2"]],
               sqrt(free$Omega[["tau_U2"]] * free$Omega[["tau_V2"]]),
               tolerance = 1e-12)
  # the product -- the part the data determine -- is what the tie preserves
  expect_equal(tied$Omega[["tau_U2"]] * tied$Omega[["tau_V2"]],
               free$Omega[["tau_U2"]] * free$Omega[["tau_V2"]],
               tolerance = 1e-12)

  # nothing else moves
  others <- setdiff(names(tied$Omega), c("tau_U2", "tau_V2"))
  expect_equal(tied$Omega[others], free$Omega[others], tolerance = 1e-12)

  # and the split the data would have produced is still reported, so the
  # direction the iteration is being pulled in stays visible
  expect_equal(tied$latent_ratio,
               free$Omega[["tau_U2"]] / free$Omega[["tau_V2"]],
               tolerance = 1e-12)
  expect_equal(free$latent_ratio, tied$latent_ratio, tolerance = 1e-12)
})


test_that("undefined components are NA and their penalties cannot matter", {
  d <- sim_ame_panel(Pr = 0, Pc = 0, Pd = 0, seed = 63)
  p <- random_params(d)
  eb <- .eb_update_variances(p, d$Y, NULL, NULL, NULL,
                             lambda = unit_penalty(), gamma = unit_penalty())
  expect_true(is.na(eb$Omega[["sigma_beta2"]]))
  expect_true(is.na(eb$Omega[["tau_beta2"]]))
  expect_equal(unname(eb$lambda[["beta"]]), 1)
  expect_equal(unname(eb$gamma[["beta"]]), 1)

  # the fallback provably has no effect: Q is unchanged by any beta penalty
  q1 <- .ame_objective(p, d$Y, NULL, NULL, NULL, eb$lambda, eb$gamma)$Q
  l2 <- eb$lambda; l2[["beta"]] <- 999
  g2 <- eb$gamma;  g2[["beta"]] <- 999
  expect_equal(.ame_objective(p, d$Y, NULL, NULL, NULL, l2, g2)$Q, q1,
               tolerance = 1e-12)

  # a single period leaves every innovation variance undefined
  p1 <- list(a = p$a[, 1, drop = FALSE], b = p$b[, 1, drop = FALSE],
             beta = p$beta[, 1, drop = FALSE], U = p$U[1], V = p$V[1])
  e1 <- .eb_update_variances(p1, d$Y[1], NULL, NULL, NULL,
                             lambda = unit_penalty(), gamma = unit_penalty())
  expect_true(all(is.na(e1$Omega[grep("^tau", names(e1$Omega))])))
  expect_true(all(e1$gamma == 1))
  expect_true(all(is.finite(e1$Omega[c("sigma_eps2", "sigma_a2", "sigma_b2",
                                       "sigma_U2", "sigma_V2")])))
})


test_that("variances are floored and penalties capped without NaN", {
  d <- sim_ame_panel(seed = 64)
  p <- random_params(d)
  pc <- p
  for (t in 2:d$Tn) {
    pc$a[, t] <- pc$a[, 1]; pc$b[, t] <- pc$b[, 1]
    pc$beta[, t] <- pc$beta[, 1]
    pc$U[[t]] <- pc$U[[1]]; pc$V[[t]] <- pc$V[[1]]
  }
  eb <- .eb_update_variances(pc, d$Y, d$X_row, d$X_col, d$X_dyad,
                             lambda = unit_penalty(), gamma = unit_penalty(),
                             eps_var = 1e-8, max_penalty = 1e6)
  expect_true(all(is.finite(eb$gamma)))
  expect_true(all(eb$gamma <= 1e6))
  expect_false(any(is.nan(eb$Omega)))

  # capping never distorts the reported variances
  expect_true(all(eb$Omega[!is.na(eb$Omega)] >= 1e-8))
})


test_that("the variance components respond to the U/V scale gauge as the spec says", {
  d <- sim_ame_panel(seed = 65)
  p <- random_params(d)
  base <- .eb_update_variances(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                               lambda = unit_penalty(), gamma = unit_penalty())
  cc <- 2.5
  ps <- p
  ps$U <- lapply(ps$U, function(U) cc * U)
  ps$V <- lapply(ps$V, function(V) V / cc)
  sc <- .eb_update_variances(ps, d$Y, d$X_row, d$X_col, d$X_dyad,
                             lambda = unit_penalty(), gamma = unit_penalty())

  # the fit is unchanged, so the observation variance must be too
  expect_equal(sc$Omega[["sigma_eps2"]], base$Omega[["sigma_eps2"]],
               tolerance = 1e-10)
  # this is exactly the drift Step 8 exists to prevent
  expect_gt(sc$Omega[["sigma_U2"]], base$Omega[["sigma_U2"]])
  expect_lt(sc$Omega[["sigma_V2"]], base$Omega[["sigma_V2"]])
})


test_that("an entirely unobserved panel is rejected", {
  d <- sim_ame_panel(seed = 66)
  Yna <- lapply(d$Y, function(Y) {
    Y[] <- NA_real_; Y
  })
  expect_error(
    .eb_update_variances(random_params(d), Yna, d$X_row, d$X_col, d$X_dyad,
                         lambda = unit_penalty(), gamma = unit_penalty()),
    "No observed cells"
  )
})
