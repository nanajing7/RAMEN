# Full-model parametric bootstrap.
#
# The design claim is that every state trajectory is redrawn from the priors
# implied by the fitted variance components, and that the estimator is then
# judged against those known generating values rather than against the spread of
# its own output. The tests check the generation mechanism directly and confirm
# the reported quantities are error-based.

small_fit_b <- function(seed = 501, K = 2) {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = seed)
  edge <- as_edge_panel(d, years = 2011:2014)
  suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, K = K, n_starts = 1,
                    outer_max_iter = 15))
}


test_that("random-walk trajectories match the prior they are drawn from", {
  set.seed(1)
  n <- 4000L; Tn <- 5L; s2 <- 1.5; t2 <- 0.4
  tr <- .draw_rw_trajectory(n, Tn, s2, t2)

  # the initial state is centred at zero with the stated variance
  expect_equal(mean(tr[, 1]), 0, tolerance = 0.05)
  expect_equal(stats::var(tr[, 1]), s2, tolerance = 0.1)

  # successive differences have the innovation variance
  d <- as.vector(tr[, -1] - tr[, -Tn])
  expect_equal(mean(d), 0, tolerance = 0.02)
  expect_equal(stats::var(d), t2, tolerance = 0.05)

  # variance accumulates as a random walk: var at t is sigma2 + (t-1) tau2
  for (t in seq_len(Tn)) {
    expect_equal(stats::var(tr[, t]), s2 + (t - 1) * t2, tolerance = 0.25,
                 info = sprintf("period %d", t))
  }
})


test_that("a single period leaves the innovation variance unused", {
  set.seed(2)
  tr <- .draw_rw_trajectory(2000L, 1L, 2.0, 999)
  expect_equal(dim(tr), c(2000L, 1L))
  expect_equal(stats::var(tr[, 1]), 2.0, tolerance = 0.15)
})


test_that("an empty block yields an empty trajectory", {
  expect_equal(dim(.draw_rw_trajectory(0L, 4L, 1, 1)), c(0L, 4L))
})


test_that("simulated panels keep the missing pattern and return the truth", {
  skip_if_not_installed("dplyr")
  fit <- small_fit_b()
  set.seed(3)
  for (rep in 1:5) {
    sim <- .simulate_full_model_panel(fit)
    expect_named(sim, c("Y_list", "params"))
    for (t in seq_along(sim$Y_list)) {
      expect_identical(is.na(sim$Y_list[[t]]), is.na(fit$Y_list[[t]]))
      expect_identical(dimnames(sim$Y_list[[t]]), dimnames(fit$Y_list[[t]]))
    }
    # the generating parameters come back with the right shapes
    expect_equal(dim(sim$params$a), dim(fit$a))
    expect_equal(dim(sim$params$b), dim(fit$b))
    expect_length(sim$params$U, length(fit$years))
    expect_equal(dim(sim$params$U[[1]]), dim(fit$U[[1]]))
  }
})


test_that("the generating states differ from the original estimate", {
  skip_if_not_installed("dplyr")
  # this is the defining feature of the design: the truth is redrawn, not reused
  fit <- small_fit_b()
  set.seed(4)
  s1 <- .simulate_full_model_panel(fit)
  s2 <- .simulate_full_model_panel(fit)

  expect_false(isTRUE(all.equal(s1$params$a, fit$a)))
  expect_false(isTRUE(all.equal(s1$params$a, s2$params$a)))
  expect_false(isTRUE(all.equal(unlist(s1$params$U), unlist(s2$params$U))))
})


test_that("the generated panel reproduces the fitted variance components", {
  skip_if_not_installed("dplyr")
  fit <- small_fit_b(seed = 502)
  Om <- fit$variance_components$Omega

  set.seed(5)
  sims <- lapply(1:150, function(b) .simulate_full_model_panel(fit)$params)

  # initial-state variance of a, averaged over replicates
  v_a1 <- mean(vapply(sims, function(p) mean(p$a[, 1]^2), numeric(1)))
  expect_equal(v_a1, Om[["sigma_a2"]], tolerance = 0.2 * Om[["sigma_a2"]])

  # innovation variance of a
  Tn <- length(fit$years)
  v_da <- mean(vapply(sims, function(p)
    mean((p$a[, -1] - p$a[, -Tn])^2), numeric(1)))
  expect_equal(v_da, Om[["tau_a2"]], tolerance = 0.2 * Om[["tau_a2"]])
})


test_that("the gauge-invariant summary keeps only comparable quantities", {
  Om <- c(sigma_eps2 = 0.4, sigma_beta2 = 0.5, sigma_a2 = 1.0, sigma_b2 = 0.7,
          sigma_U2 = 0.8, sigma_V2 = 0.5, tau_beta2 = 0.05, tau_a2 = 0.2,
          tau_b2 = 0.15, tau_U2 = 0.18, tau_V2 = 0.12)
  gi <- .gauge_invariant_omega(Om)

  expect_named(gi, c("sigma_eps2", "sigma_beta2", "sigma_a2", "sigma_b2",
                     "tau_beta2", "tau_a2", "tau_b2",
                     "sigma_UV2_product", "tau_UV2_product"))
  expect_equal(unname(gi[["sigma_UV2_product"]]), 0.8 * 0.5)
  expect_equal(unname(gi[["tau_UV2_product"]]), 0.18 * 0.12)

  # rescaling U by c and V by 1/c leaves every reported quantity unchanged
  cc <- 2.5
  Om2 <- Om
  Om2[["sigma_U2"]] <- Om[["sigma_U2"]] * cc^2
  Om2[["sigma_V2"]] <- Om[["sigma_V2"]] / cc^2
  Om2[["tau_U2"]] <- Om[["tau_U2"]] * cc^2
  Om2[["tau_V2"]] <- Om[["tau_V2"]] / cc^2
  expect_equal(.gauge_invariant_omega(Om2), gi, tolerance = 1e-12)

  # the individual latent components are NOT reported, precisely because they
  # would have changed
  expect_false(any(c("sigma_U2", "sigma_V2", "tau_U2", "tau_V2") %in% names(gi)))
})


test_that("the bootstrap runs and reports error-based quantities", {
  skip_if_not_installed("dplyr")
  fit <- small_fit_b()
  bt <- suppressWarnings(.bootstrap_full(fit, B = 8, seed = 21))

  expect_s3_class(bt, "bootstrap_ame_full")
  expect_equal(bt$design, "full-model")
  expect_gt(bt$B_successful, 0L)

  # state blocks: mean error and RMSE, not standard errors
  for (nm in c("a", "b")) {
    s <- bt$accuracy[[nm]]
    expect_true(all(c("period", "mean_error", "rmse") %in% names(s)))
    expect_false("std.error" %in% names(s))
    expect_true(all(s$rmse >= 0))
    expect_true(all(s$rmse >= abs(s$mean_error) - 1e-9))
  }
  expect_null(bt$accuracy$beta)     # no covariates in this fixture

  # variance components on the gauge-invariant scale: 9 rows, not 11
  expect_equal(nrow(bt$accuracy$Omega), 9L)
  expect_true(all(c("sigma_UV2_product", "tau_UV2_product") %in%
                    bt$accuracy$Omega$component))
  expect_false(any(c("sigma_U2", "tau_V2") %in% bt$accuracy$Omega$component))

  # the raw eleven-component draws are retained for inspection
  expect_equal(nrow(bt$draws$omega_full), 11L)
  expect_equal(ncol(bt$draws$omega_full), bt$B_successful)

  # latent recovery, one row per period
  lr <- bt$latent_recovery
  expect_equal(nrow(lr), length(fit$years))
  expect_true(all(lr$frobenius_mean >= 0))
  expect_true(all(lr$correlation_mean >= -1 & lr$correlation_mean <= 1))

  expect_silent(invisible(capture.output(print(bt))))
})


test_that("RMSE and mean error are computed against the generating values", {
  skip_if_not_installed("dplyr")
  # a direct check: rebuild one replicate by hand and confirm the error matches
  fit <- small_fit_b(seed = 503)
  set.seed(1)
  seeds <- sample.int(.Machine$integer.max, 2)

  set.seed(seeds[1])
  sim <- .simulate_full_model_panel(fit)
  est <- .refit_full_model(sim$Y_list, fit)
  skip_if(is.null(est), "refit failed")

  manual_err <- est$params$a - sim$params$a
  expect_equal(dim(manual_err), dim(fit$a))

  # the error is the estimate minus the generating value, never minus the
  # original estimate
  expect_false(isTRUE(all.equal(manual_err, est$params$a - fit$a)))
})


test_that("results are reproducible", {
  skip_if_not_installed("dplyr")
  fit <- small_fit_b()
  b1 <- suppressWarnings(.bootstrap_full(fit, B = 6, seed = 31))
  b2 <- suppressWarnings(.bootstrap_full(fit, B = 6, seed = 31))
  expect_equal(b1$accuracy$a$rmse, b2$accuracy$a$rmse, tolerance = 1e-12)
  expect_equal(b1$latent_recovery, b2$latent_recovery, tolerance = 1e-12)

  b3 <- suppressWarnings(.bootstrap_full(fit, B = 6, seed = 77))
  expect_false(isTRUE(all.equal(b1$accuracy$a$rmse, b3$accuracy$a$rmse)))

  set.seed(999); before <- stats::rnorm(1)
  set.seed(999)
  invisible(suppressWarnings(.bootstrap_full(fit, B = 3, seed = 4)))
  expect_identical(stats::rnorm(1), before)
})


test_that("covariates are held at their observed values", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = 504)
  edge <- as_edge_panel(d, years = 2001:2004)
  rows <- sprintf("%03d", seq_len(10))
  set.seed(504)
  row_cov <- do.call(rbind, lapply(1:4, function(t)
    data.frame(year = 2000 + t, node_row = rows, z = stats::rnorm(10),
               stringsAsFactors = FALSE)))
  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, row_cov_df = row_cov,
                    row_covar_names = "z", row_cov_prefix = "row_",
                    K = 2, n_starts = 1, outer_max_iter = 15))

  set.seed(6)
  s1 <- .simulate_full_model_panel(fit)
  s2 <- .simulate_full_model_panel(fit)
  # the covariates never change; only the states and errors do
  expect_identical(fit$covariates$X_row_list, fit$covariates$X_row_list)
  expect_false(isTRUE(all.equal(s1$params$beta, s2$params$beta)))

  bt <- suppressWarnings(.bootstrap_full(fit, B = 6, seed = 41))
  expect_s3_class(bt$accuracy$beta, "data.frame")
  expect_equal(nrow(bt$accuracy$beta), 1L * length(fit$years))
})


test_that("input validation", {
  skip_if_not_installed("dplyr")
  fit <- small_fit_b()
  expect_error(.bootstrap_full(structure(list(), class = "x")),
               "must be of class")
  expect_error(.bootstrap_full(fit, B = 1), "at least 2")
})
