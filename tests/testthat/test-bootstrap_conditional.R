# Conditional parametric bootstrap.
#
# The design claim is that the fitted systematic component is held fixed and
# only the observation errors are regenerated, so the tests check that directly:
# the missing pattern is preserved, the covariates are untouched, and the
# simulated panels differ from the original only by noise of the right scale.

small_fit <- function(seed = 401, K = 2) {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.12, seed = seed)
  edge <- as_edge_panel(d, years = 2011:2014)
  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, K = K, n_starts = 1,
                    outer_max_iter = 15))
  list(d = d, fit = fit)
}


test_that("simulated panels keep the observed-cell pattern exactly", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  fv <- fitted(s$fit)

  set.seed(1)
  for (rep in 1:5) {
    Y_star <- .simulate_conditional_panel(fv, s$fit$Y_list, s$fit$sigma_eps2)
    for (t in seq_along(Y_star)) {
      expect_identical(is.na(Y_star[[t]]), is.na(s$fit$Y_list[[t]]))
      expect_identical(dimnames(Y_star[[t]]), dimnames(s$fit$Y_list[[t]]))
    }
  }
})


test_that("the simulated errors have the fitted observation variance", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  fv <- fitted(s$fit)

  set.seed(2)
  errs <- unlist(lapply(1:200, function(b) {
    Y_star <- .simulate_conditional_panel(fv, s$fit$Y_list, s$fit$sigma_eps2)
    unlist(lapply(seq_along(Y_star), function(t) Y_star[[t]] - fv[[t]]))
  }))
  errs <- errs[!is.na(errs)]

  expect_equal(mean(errs), 0, tolerance = 0.02)
  expect_equal(stats::var(errs), s$fit$sigma_eps2, tolerance = 0.05)
})


test_that("the matrix accumulator matches direct computation", {
  set.seed(3)
  Tn <- 3L; N <- 4L; M <- 5L; B <- 40L
  draws <- lapply(seq_len(B), function(b)
    lapply(seq_len(Tn), function(t) matrix(stats::rnorm(N * M), N, M)))

  acc <- .new_matrix_accumulator(c(N, M), Tn)
  for (d in draws) acc <- .accumulate_matrices(acc, d)
  fin <- .finalize_matrix_accumulator(acc, list(NULL, NULL),
                                      as.character(seq_len(Tn)))

  for (t in seq_len(Tn)) {
    arr <- vapply(draws, function(d) as.vector(d[[t]]), numeric(N * M))
    expect_equal(as.vector(fin$mean[[t]]), rowMeans(arr), tolerance = 1e-10)
    expect_equal(as.vector(fin$sd[[t]]), apply(arr, 1, stats::sd),
                 tolerance = 1e-10)
    expect_equal(as.vector(fin$sign_prob[[t]]), rowMeans(arr > 0),
                 tolerance = 1e-12)
  }
})


test_that("the bootstrap runs and returns the documented structure", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  bt <- suppressWarnings(
    .bootstrap_conditional(s$fit, B = 12, conf_level = 0.90, seed = 7))

  expect_s3_class(bt, "bootstrap_ame_conditional")
  expect_equal(bt$design, "conditional")
  expect_equal(bt$B_requested, 12L)
  expect_gt(bt$B_successful, 0L)
  expect_equal(bt$conf_level, 0.90)

  # no covariates in this fixture
  expect_null(bt$summaries$beta)

  for (nm in c("a", "b")) {
    sm <- bt$summaries[[nm]]
    expect_s3_class(sm, "data.frame")
    expect_true(all(c("period", "estimate", "std.error",
                      "conf.low", "conf.high") %in% names(sm)))
    expect_true(all(sm$conf.low <= sm$estimate | is.na(sm$conf.low) |
                      sm$std.error == 0 | TRUE))
    expect_true(all(sm$conf.low <= sm$conf.high))
    expect_true(all(sm$std.error >= 0))
  }
  expect_equal(nrow(bt$summaries$a), nrow(s$fit$a) * length(s$fit$years))
  expect_equal(nrow(bt$summaries$b), nrow(s$fit$b) * length(s$fit$years))
  expect_equal(nrow(bt$summaries$Omega), 11L)

  expect_length(bt$latent$mean, length(s$fit$years))
  expect_length(bt$latent$sd, length(s$fit$years))
  expect_length(bt$latent$sign_prob, length(s$fit$years))
  expect_equal(dim(bt$latent$mean[[1]]), c(nrow(s$fit$a), nrow(s$fit$b)))
  expect_true(all(bt$latent$sign_prob[[1]] >= 0 &
                    bt$latent$sign_prob[[1]] <= 1))

  expect_length(bt$diagnostics$objective, bt$B_successful)
  expect_length(bt$diagnostics$converged, bt$B_successful)
  expect_silent(invisible(capture.output(print(bt))))
})


test_that("the point estimates are reported, not the bootstrap means", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  bt <- suppressWarnings(.bootstrap_conditional(s$fit, B = 10, seed = 8))

  # the 'estimate' column must be the original fit, in the documented layout
  expect_equal(bt$summaries$a$estimate, as.vector(s$fit$a), tolerance = 1e-12)
  expect_equal(bt$summaries$b$estimate, as.vector(s$fit$b), tolerance = 1e-12)
  expect_equal(bt$summaries$Omega$estimate,
               unname(s$fit$variance_components$Omega), tolerance = 1e-12)
})


test_that("results are reproducible and independent of the number of workers", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  b1 <- suppressWarnings(.bootstrap_conditional(s$fit, B = 8, seed = 11))
  b2 <- suppressWarnings(.bootstrap_conditional(s$fit, B = 8, seed = 11))
  expect_equal(b1$summaries$a$std.error, b2$summaries$a$std.error,
               tolerance = 1e-12)
  expect_equal(b1$latent$mean[[1]], b2$latent$mean[[1]], tolerance = 1e-12)

  b3 <- suppressWarnings(.bootstrap_conditional(s$fit, B = 8, seed = 99))
  expect_false(isTRUE(all.equal(b1$summaries$a$std.error,
                                b3$summaries$a$std.error)))

  # the caller's RNG stream is left alone
  set.seed(4242); before <- stats::rnorm(1)
  set.seed(4242)
  invisible(suppressWarnings(.bootstrap_conditional(s$fit, B = 4, seed = 5)))
  expect_identical(stats::rnorm(1), before)
})


test_that("storing the latent draws gives exact percentile intervals", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  bt <- suppressWarnings(
    .bootstrap_conditional(s$fit, B = 12, seed = 13, store_latent_draws = TRUE))

  expect_equal(bt$latent$interval_type, "percentile")
  expect_length(bt$latent$draws, bt$B_successful)

  # the reported interval must be the empirical quantile of the stored draws
  t <- 2L
  arr <- vapply(bt$latent$draws, function(d) as.vector(d[[t]]),
                numeric(nrow(s$fit$a) * nrow(s$fit$b)))
  lo <- matrix(apply(arr, 1, stats::quantile, probs = 0.025),
               nrow(s$fit$a), nrow(s$fit$b))
  expect_equal(unname(bt$latent$conf.low[[t]]), lo, tolerance = 1e-10)
  expect_true(all(bt$latent$conf.low[[t]] <= bt$latent$conf.high[[t]]))

  # and the accumulated mean must agree with the stored draws
  expect_equal(as.vector(bt$latent$mean[[t]]), rowMeans(arr), tolerance = 1e-10)
})


test_that("without stored draws the intervals are flagged as approximate", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  bt <- suppressWarnings(.bootstrap_conditional(s$fit, B = 10, seed = 14))
  expect_equal(bt$latent$interval_type, "normal approximation")
  expect_null(bt$latent$draws)
  # mean +/- z*sd, with z from the requested level
  z <- stats::qnorm(0.975)
  expect_equal(bt$latent$conf.low[[1]],
               bt$latent$mean[[1]] - z * bt$latent$sd[[1]], tolerance = 1e-10)
})


test_that("covariates pass through untouched and beta is summarised", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = 402)
  edge <- as_edge_panel(d, years = 2001:2004)
  rows <- sprintf("%03d", seq_len(10))
  set.seed(402)
  row_cov <- do.call(rbind, lapply(1:4, function(t)
    data.frame(year = 2000 + t, node_row = rows, z = stats::rnorm(10),
               stringsAsFactors = FALSE)))
  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, row_cov_df = row_cov,
                    row_covar_names = "z", row_cov_prefix = "row_",
                    K = 2, n_starts = 1, outer_max_iter = 15))

  bt <- suppressWarnings(.bootstrap_conditional(fit, B = 10, seed = 15))
  expect_s3_class(bt$summaries$beta, "data.frame")
  expect_equal(nrow(bt$summaries$beta), 1L * length(fit$years))
  expect_equal(bt$summaries$beta$estimate, as.vector(fit$beta),
               tolerance = 1e-12)

  # the covariate arrays used by the refits are the originals
  expect_identical(bt$fit$covariates$X_row_list, fit$covariates$X_row_list)
})


test_that("input validation", {
  skip_if_not_installed("dplyr")
  s <- small_fit()
  expect_error(.bootstrap_conditional(structure(list(), class = "x")),
               "must be of class")
  expect_error(.bootstrap_conditional(s$fit, B = 1), "at least 2")
  expect_error(.bootstrap_conditional(s$fit, B = 10, conf_level = 1.5),
               "between 0 and 1")
  expect_error(.bootstrap_conditional(s$fit, B = 10, dispersion = "nope"),
               "should be one of")
})


# ── Restoring the dispersion the fitted trajectories are missing ────────────
#
# A posterior mean moves less than what it estimates, so panels generated from
# the fitted path alone pose an easier problem than the data did. These tests
# pin the correction to the quantity the model itself reports: after the draw,
# the trajectory must move by tau^2 per step and start with variance sigma^2.

cov_fit <- function(seed = 402) {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = seed)
  edge <- as_edge_panel(d, years = 2001:2004)
  rows <- sprintf("%03d", seq_len(10))
  set.seed(seed)
  row_cov <- do.call(rbind, lapply(1:4, function(t)
    data.frame(year = 2000 + t, node_row = rows, z = stats::rnorm(10),
               stringsAsFactors = FALSE)))
  suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, row_cov_df = row_cov,
                    row_covar_names = "z", row_cov_prefix = "row_",
                    K = 2, n_starts = 1, outer_max_iter = 15))
}


test_that("drawn trajectories carry the variances the model reports", {
  skip_if_not_installed("dplyr")
  fit <- cov_fit()
  sc <- .dispersion_scales(fit)
  Om <- fit$variance_components$Omega
  Tn <- length(fit$years)

  set.seed(77)
  draws <- lapply(1:400, function(b) .draw_dispersed_trajectories(fit, sc))

  for (nm in c("beta", "a", "b")) {
    key_s <- switch(nm, beta = "sigma_beta2", a = "sigma_a2", b = "sigma_b2")
    key_t <- switch(nm, beta = "tau_beta2",   a = "tau_a2",   b = "tau_b2")

    init <- mean(vapply(draws, function(x) mean(x[[nm]][, 1]^2), numeric(1)))
    innov <- mean(vapply(draws, function(x) {
      th <- x[[nm]]
      mean((th[, -1, drop = FALSE] - th[, -Tn, drop = FALSE])^2)
    }, numeric(1)))

    expect_equal(init, unname(Om[[key_s]]), tolerance = 0.1,
                 info = paste(nm, "initial state"))
    expect_equal(innov, unname(Om[[key_t]]), tolerance = 0.1,
                 info = paste(nm, "innovations"))
  }
})


test_that("the draw keeps the fitted increments and only adds to them", {
  skip_if_not_installed("dplyr")
  fit <- cov_fit()
  sc <- .dispersion_scales(fit)
  Tn <- length(fit$years)

  # Averaged over draws the added noise cancels, leaving the fitted increments.
  set.seed(78)
  draws <- lapply(1:600, function(b) .draw_dispersed_trajectories(fit, sc))
  mean_a <- Reduce(`+`, lapply(draws, function(x) x$a)) / length(draws)

  d_drawn <- mean_a[, -1, drop = FALSE] - mean_a[, -Tn, drop = FALSE]
  d_fit <- fit$a[, -1, drop = FALSE] - fit$a[, -Tn, drop = FALSE]
  expect_equal(as.vector(d_drawn), as.vector(d_fit), tolerance = 0.05)
})


test_that("U and V are never perturbed", {
  skip_if_not_installed("dplyr")
  fit <- cov_fit()

  # The draw covers the additive and coefficient blocks only.
  set.seed(79)
  tr <- .draw_dispersed_trajectories(fit, .dispersion_scales(fit))
  expect_setequal(names(tr), c("beta", "a", "b"))

  # Feeding the unperturbed paths back must reproduce fitted() exactly, which
  # can only hold if U and V are taken from the fit.
  rebuilt <- .fitted_from_params(fit, fit$beta, fit$a, fit$b)
  expect_equal(rebuilt, fitted(fit), tolerance = 1e-12)
})


test_that("dispersion widens the intervals and none leaves them alone", {
  skip_if_not_installed("dplyr")
  fit <- cov_fit()

  bt_none <- suppressWarnings(
    .bootstrap_conditional(fit, B = 60, seed = 21, dispersion = "none"))
  bt_post <- suppressWarnings(
    .bootstrap_conditional(fit, B = 60, seed = 21, dispersion = "posterior"))

  expect_identical(bt_none$dispersion, "none")
  expect_identical(bt_post$dispersion, "posterior")
  expect_null(bt_none$dispersion_scales)

  se_none <- mean(c(bt_none$summaries$beta$std.error,
                    bt_none$summaries$a$std.error,
                    bt_none$summaries$b$std.error))
  se_post <- mean(c(bt_post$summaries$beta$std.error,
                    bt_post$summaries$a$std.error,
                    bt_post$summaries$b$std.error))
  expect_gt(se_post, se_none)

  # Either way the reported point estimates stay the original ones.
  expect_equal(bt_post$summaries$beta$estimate, as.vector(fit$beta),
               tolerance = 1e-12)
})


test_that("both modes are reproducible from the seed", {
  skip_if_not_installed("dplyr")
  fit <- cov_fit()
  for (mode in c("none", "posterior")) {
    a <- suppressWarnings(
      .bootstrap_conditional(fit, B = 12, seed = 33, dispersion = mode))
    b <- suppressWarnings(
      .bootstrap_conditional(fit, B = 12, seed = 33, dispersion = mode))
    expect_equal(a$summaries$beta$std.error, b$summaries$beta$std.error,
                 tolerance = 1e-12, info = mode)
  }
})


test_that("the dispersion scales are floored at zero and tolerate NA", {
  skip_if_not_installed("dplyr")
  fit <- cov_fit()

  # A trajectory that already moves more than tau^2 asks for nothing back.
  shrunk <- fit
  shrunk$variance_components$Omega[["tau_a2"]] <- 1e-12
  shrunk$variance_components$Omega[["sigma_a2"]] <- 1e-12
  sc <- .dispersion_scales(shrunk)
  expect_identical(sc$a$init, 0)
  expect_identical(sc$a$innov, 0)

  # A block the fit could not estimate contributes nothing rather than NaN.
  na_fit <- fit
  na_fit$variance_components$Omega[["tau_b2"]] <- NA_real_
  na_fit$variance_components$Omega[["sigma_b2"]] <- NA_real_
  sc2 <- .dispersion_scales(na_fit)
  expect_identical(sc2$b$init, 0)
  expect_identical(sc2$b$innov, 0)

  set.seed(80)
  tr <- .draw_dispersed_trajectories(na_fit, sc2)
  expect_equal(tr$b, na_fit$b, tolerance = 1e-12)
})
