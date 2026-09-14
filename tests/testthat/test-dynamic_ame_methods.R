# Methods and accessors for the dynamic_ame class.

fitted_model <- function(seed = 301) {
  d <- sim_ame_panel(N = 12, M = 9, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.08, seed = seed)
  edge <- as_edge_panel(d, years = 2011:2014)
  rows <- sprintf("%03d", seq_len(d$N))
  set.seed(seed)
  row_cov <- do.call(rbind, lapply(seq_len(d$Tn), function(t)
    data.frame(year = 2010 + t, node_row = rows, z = stats::rnorm(d$N),
               stringsAsFactors = FALSE)))
  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, row_cov_df = row_cov,
                    row_covar_names = "z", row_cov_prefix = "row_",
                    K = 2, n_starts = 1, outer_max_iter = 25))
  list(d = d, fit = fit)
}


test_that("print and summary report without error and cover the essentials", {
  skip_if_not_installed("dplyr")
  fm <- fitted_model()
  out <- capture.output(print(fm$fit))
  expect_true(any(grepl("Dynamic bipartite AME", out)))
  expect_true(any(grepl("senders", out)))
  expect_true(any(grepl("Variance components", out)))

  s <- summary(fm$fit)
  expect_s3_class(s, "summary.dynamic_ame")
  expect_equal(nrow(s$components), 5L)
  expect_named(s$components, c("block", "sigma2", "tau2", "lambda", "gamma"))
  expect_equal(s$components$block, c("beta", "a", "b", "U", "V"))
  expect_equal(nrow(s$coef_summary), 1L)

  sout <- capture.output(print(s))
  expect_true(any(grepl("Variance components and the penalties", sout)))
  expect_true(any(grepl("Convergence", sout)))
})


test_that("summary carries no coefficient block when there are no covariates", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 10, M = 8, Tn = 4, Pr = 0, Pc = 0, Pd = 0, seed = 302)
  fit <- suppressWarnings(
    fit_dynamic_ame(as_edge_panel(d), K = 2, n_starts = 1, outer_max_iter = 20))
  s <- summary(fit)
  expect_null(s$coef_summary)
  expect_true(is.na(s$components$sigma2[s$components$block == "beta"]))
  expect_silent(invisible(capture.output(print(s))))
})


test_that("fitted() reproduces the model equation cell by cell", {
  skip_if_not_installed("dplyr")
  fm <- fitted_model(303)
  fit <- fm$fit
  fv <- fitted(fit)

  expect_length(fv, length(fit$years))
  expect_named(fv, as.character(fit$years))
  expect_equal(dimnames(fv[[1]]), list(fit$row_ids, fit$col_ids))

  t <- 2L
  Xr <- fit$covariates$X_row_list[[t]]
  manual <- matrix(0, nrow(fit$a), nrow(fit$b))
  for (i in seq_len(nrow(fit$a))) {
    for (j in seq_len(nrow(fit$b))) {
      manual[i, j] <- fit$a[i, t] + fit$b[j, t] +
        sum(Xr[i, ] * fit$beta[, t]) +
        sum(fit$U[[t]][i, ] * fit$V[[t]][j, ])
    }
  }
  expect_equal(unname(fv[[t]]), manual, tolerance = 1e-10)
})


test_that("residuals() equals Y minus fitted and keeps the missing mask", {
  skip_if_not_installed("dplyr")
  fm <- fitted_model(304)
  fit <- fm$fit
  r <- residuals(fit)
  fv <- fitted(fit)
  for (t in seq_along(r)) {
    expect_equal(r[[t]], fit$Y_list[[t]] - fv[[t]], tolerance = 1e-12)
    expect_identical(is.na(r[[t]]), is.na(fit$Y_list[[t]]))
  }
  # the residual variance the fit reports is the mean square of these
  ssr <- sum(vapply(r, function(x) sum(x^2, na.rm = TRUE), numeric(1)))
  nobs <- sum(vapply(r, function(x) sum(!is.na(x)), integer(1)))
  expect_equal(ssr / nobs, fit$sigma_eps2, tolerance = 1e-8)
})


test_that("the accessors return what they promise", {
  skip_if_not_installed("dplyr")
  fm <- fitted_model(305)
  fit <- fm$fit

  vc <- variance_components(fit)
  expect_named(vc, c("Omega", "lambda", "gamma"))
  expect_length(vc$Omega, 11L)
  expect_identical(vc, fit$variance_components)

  ct <- coef_trajectory(fit)
  expect_equal(dim(ct), c(1L, length(fit$years)))
  expect_equal(rownames(ct), "z")
  ctl <- coef_trajectory(fit, long = TRUE)
  expect_s3_class(ctl, "data.frame")
  expect_equal(ctl$coefficient, unname(ct["z", ]), tolerance = 1e-12)

  lt <- latent_trajectory(fit)
  expect_length(lt, length(fit$years))
  expect_equal(unname(lt[[3]]), unname(fit$U[[3]] %*% t(fit$V[[3]])),
               tolerance = 1e-12)
  expect_equal(dimnames(lt[[1]]), list(fit$row_ids, fit$col_ids))

  ci <- convergence_info(fit)
  expect_true(is.logical(ci$converged))
  expect_length(ci$objective_trace, ci$outer_iterations)
})


test_that("coef_trajectory is NULL when the model has no covariates", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 10, M = 8, Tn = 4, Pr = 0, Pc = 0, Pd = 0, seed = 306)
  fit <- suppressWarnings(
    fit_dynamic_ame(as_edge_panel(d), K = 2, n_starts = 1, outer_max_iter = 20))
  expect_null(coef_trajectory(fit))
  expect_null(coef_trajectory(fit, long = TRUE))
})


test_that("decompose_fit splits the fit into parts that sum back", {
  skip_if_not_installed("dplyr")
  fm <- fitted_model(307)
  fit <- fm$fit
  dec <- decompose_fit(fit)
  fv <- fitted(fit)

  expect_length(dec, length(fit$years))
  for (t in seq_along(dec)) {
    expect_named(dec[[t]], c("additive", "covariate", "latent", "fitted"))
    expect_equal(dec[[t]]$additive + dec[[t]]$covariate + dec[[t]]$latent,
                 dec[[t]]$fitted, tolerance = 1e-12)
    expect_equal(dec[[t]]$fitted, fv[[t]], tolerance = 1e-10)
    expect_equal(dec[[t]]$latent, latent_trajectory(fit)[[t]], tolerance = 1e-12)
  }
})


test_that("every method rejects objects of the wrong class", {
  bad <- structure(list(), class = "not_a_fit")
  expect_error(print.dynamic_ame(bad), "must be of class")
  expect_error(summary.dynamic_ame(bad), "must be of class")
  expect_error(fitted.dynamic_ame(bad), "must be of class")
  expect_error(residuals.dynamic_ame(bad), "must be of class")
  expect_error(variance_components.dynamic_ame(bad), "must be of class")
  expect_error(coef_trajectory.dynamic_ame(bad), "must be of class")
  expect_error(latent_trajectory.dynamic_ame(bad), "must be of class")
  expect_error(decompose_fit.dynamic_ame(bad), "must be of class")
  expect_error(convergence_info.dynamic_ame(bad), "must be of class")
})
