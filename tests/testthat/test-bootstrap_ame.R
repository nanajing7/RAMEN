# The public bootstrap entry point.
#
# The two designs are tested in detail in their own files; here the concern is
# that the wrapper dispatches correctly, keeps the two results separate, and
# stays reproducible.

wrapper_fit <- function(seed = 601) {
  d <- sim_ame_panel(N = 9, M = 7, Tn = 3, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = seed)
  edge <- as_edge_panel(d, years = 2011:2013)
  suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, K = 2, n_starts = 1,
                    outer_max_iter = 10))
}


test_that("both designs run and are returned separately", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  bt <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 5, B_full = 4, seed = 1))

  expect_s3_class(bt, "bootstrap_ame")
  expect_equal(bt$design, "both")
  expect_s3_class(bt$conditional, "bootstrap_ame_conditional")
  expect_s3_class(bt$full, "bootstrap_ame_full")

  # each design used its own B
  expect_equal(bt$conditional$B_requested, 5L)
  expect_equal(bt$full$B_requested, 4L)

  # the sub-objects report their own quantities
  expect_true("std.error" %in% names(bt$conditional$summaries$a))
  expect_true(all(c("mean_error", "rmse") %in% names(bt$full$accuracy$a)))
  expect_false("std.error" %in% names(bt$full$accuracy$a))
})


test_that("a single design can be requested", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()

  bc <- suppressWarnings(bootstrap_ame(fit, design = "conditional", B = 4))
  expect_s3_class(bc$conditional, "bootstrap_ame_conditional")
  expect_null(bc$full)

  bf <- suppressWarnings(bootstrap_ame(fit, design = "full", B_full = 4))
  expect_s3_class(bf$full, "bootstrap_ame_full")
  expect_null(bf$conditional)

  expect_error(bootstrap_ame(fit, design = "nonsense"), "should be one of")
})


test_that("B_full defaults to B but can be set independently", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  bt <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 4, seed = 2))
  expect_equal(bt$full$B_requested, 4L)
  expect_equal(bt$settings$B_full, 4)
})


test_that("the two designs are driven by different random numbers", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  bt <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 4, B_full = 4, seed = 3))
  # derived seeds differ, so the designs are not sharing a stream
  expect_false(identical(bt$conditional$settings$seed,
                         bt$full$settings$seed))
})


test_that("the whole call is reproducible from the base seed", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  b1 <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 4, B_full = 3, seed = 11))
  b2 <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 4, B_full = 3, seed = 11))
  expect_equal(b1$conditional$summaries$a$std.error,
               b2$conditional$summaries$a$std.error, tolerance = 1e-12)
  expect_equal(b1$full$accuracy$a$rmse, b2$full$accuracy$a$rmse,
               tolerance = 1e-12)

  b3 <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 4, B_full = 3, seed = 22))
  expect_false(isTRUE(all.equal(b1$conditional$summaries$a$std.error,
                                b3$conditional$summaries$a$std.error)))

  set.seed(7); before <- stats::rnorm(1)
  set.seed(7)
  invisible(suppressWarnings(bootstrap_ame(fit, design = "both", B = 3, B_full = 3, seed = 5)))
  expect_identical(stats::rnorm(1), before)
})


test_that("printing works for the whole object and for each design", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  bt <- suppressWarnings(bootstrap_ame(fit, design = "both", B = 4, B_full = 3, seed = 4))

  out <- capture.output(print(bt))
  expect_true(any(grepl("Conditional design", out)))
  expect_true(any(grepl("Full-model design", out)))
  expect_true(any(grepl("different quantities on purpose", out)))

  expect_silent(invisible(capture.output(print(bt$conditional))))
  expect_silent(invisible(capture.output(print(bt$full))))

  # with only one design the closing note is omitted
  bc <- suppressWarnings(bootstrap_ame(fit, design = "conditional", B = 4))
  out1 <- capture.output(print(bc))
  expect_false(any(grepl("different quantities on purpose", out1)))
})


test_that("arguments are passed through to the underlying designs", {
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  bt <- suppressWarnings(
    bootstrap_ame(fit, design = "conditional", B = 5, conf_level = 0.80,
                  store_latent_draws = TRUE, warm_start = FALSE, seed = 9))

  expect_equal(bt$conditional$conf_level, 0.80)
  expect_false(bt$conditional$warm_start)
  expect_equal(bt$conditional$latent$interval_type, "percentile")
  expect_length(bt$conditional$latent$draws, bt$conditional$B_successful)
})


test_that("a non-fit object is rejected", {
  expect_error(bootstrap_ame(structure(list(), class = "x")),
               "must be of class")
})

test_that("the default design is conditional only", {
  # The default is what an analysis gets when it does not think about it, and
  # what it should get is the design that produces standard errors and
  # confidence intervals. The full-model design reports bias and RMSE against
  # known generating values -- a question about the estimator, not about this
  # dataset -- so it must be asked for.
  skip_if_not_installed("dplyr")
  fit <- wrapper_fit()
  bd <- suppressWarnings(bootstrap_ame(fit, B = 4, seed = 7))
  expect_equal(bd$design, "conditional")
  expect_false(is.null(bd$conditional))
  expect_null(bd$full)
})
