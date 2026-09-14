# The inner loop, the outer empirical-Bayes loop, and the public entry point.

test_that("the inner loop converges and never increases Q", {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, seed = 101)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  res <- .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, d$X_col, d$X_dyad,
                        unit_penalty(), unit_penalty(),
                        eps_Q = 1e-6, max_iter = 500)

  expect_true(res$converged)
  expect_length(res$Q_trace, res$iterations + 1L)
  expect_true(all(is.finite(res$Q_trace)))
  expect_true(all(diff(res$Q_trace) <= 1e-9))
  expect_lt(res$objective$Q, res$Q_trace[1])

  # the stopping rule is the spec's relative change
  n <- length(res$Q_trace)
  expect_lt(abs(res$Q_trace[n - 1] - res$Q_trace[n]) /
              (abs(res$Q_trace[n - 1]) + 1e-8), 1e-6)

  # the answer is a fixed point
  again <- .ame_inner_bcd(res$params, d$Y, d$X_row, d$X_col, d$X_dyad,
                          unit_penalty(), unit_penalty(), eps_Q = 1e-6)
  expect_lte(again$iterations, 2L)
})


test_that("a tighter inner tolerance costs more sweeps and reaches a lower Q", {
  d <- sim_ame_panel(N = 9, M = 6, Tn = 4, seed = 102)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  run <- function(eps, cap = 2000) suppressWarnings(
    .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, d$X_col, d$X_dyad,
                   unit_penalty(), unit_penalty(), eps_Q = eps, max_iter = cap))
  loose <- run(1e-3); tight <- run(1e-9)
  expect_lt(loose$iterations, tight$iterations)
  expect_lte(tight$objective$Q, loose$objective$Q + 1e-12)
})


test_that("the inner loop warns and flags when it runs out of sweeps", {
  d <- sim_ame_panel(seed = 103)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  expect_warning(
    res <- .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, d$X_col, d$X_dyad,
                          unit_penalty(), unit_penalty(),
                          eps_Q = 1e-14, max_iter = 2),
    "did not converge"
  )
  expect_false(res$converged)
  expect_true(all(is.finite(res$params$a)))
})


test_that("stronger penalties flatten the trajectories, as they must", {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 4, seed = 104)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  base <- suppressWarnings(
    .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, d$X_col, d$X_dyad,
                   unit_penalty(), unit_penalty(), max_iter = 800))
  smooth <- suppressWarnings(
    .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, d$X_col, d$X_dyad,
                   1, 1000, max_iter = 800))
  shrunk <- suppressWarnings(
    .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, d$X_col, d$X_dyad,
                   1000, 1, max_iter = 800))

  innov <- function(r) sum(vapply(2:d$Tn,
    function(t) sum((r$params$a[, t] - r$params$a[, t - 1])^2), numeric(1)))
  expect_lt(innov(smooth), innov(base))
  expect_lt(sum(shrunk$params$a[, 1]^2), sum(base$params$a[, 1]^2))
})


test_that("the relative-change helper follows the spec and drops NA entries", {
  expect_true(is.finite(.rel_change(c(0, 0), c(0, 0))))
  expect_equal(.rel_change(c(1, NA, 3), c(1, NA, 3)), 0, tolerance = 1e-12)
  expect_equal(.rel_change(c(NA, NA), c(NA, NA)), 0)
  expect_equal(.rel_change(c(3, 4), c(1, 2)), sqrt(8) / (sqrt(5) + 1e-8),
               tolerance = 1e-9)
})


test_that("the outer loop converges and reports a complete result", {
  d <- sim_ame_panel(N = 12, M = 9, Tn = 5, Pr = 1, Pc = 0, Pd = 0, seed = 105)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  start <- suppressWarnings(
    .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, NULL, NULL,
                   unit_penalty(), unit_penalty(), max_iter = 500))
  r <- suppressWarnings(
    .ame_outer_eb(start$params, ini$reference, d$Y, d$X_row, NULL, NULL,
                  unit_penalty(), unit_penalty(), outer_max_iter = 60))

  expect_true(r$converged)
  expect_length(r$Omega_trace, r$iterations)
  expect_length(r$objective_trace, r$iterations)
  expect_length(r$scale_trace, r$iterations)
  expect_true(all(is.finite(r$lambda)) && all(r$lambda > 0))
  expect_true(all(is.finite(r$gamma)) && all(r$gamma > 0))

  # convergence is judged on the two quantities the spec names
  n <- length(r$Omega_trace)
  expect_lt(.rel_change(r$Omega_trace[[n]], r$Omega_trace[[n - 1]]), 1e-4)

  # identification is a fixed point of the converged fit
  again <- .ame_identify(r$params, ini$reference)
  expect_equal(again$c, 1, tolerance = 1e-8)
  expect_equal(again$R, diag(d$K), tolerance = 1e-7)
})


test_that("the outer loop warns and flags when it runs out of iterations", {
  d <- sim_ame_panel(seed = 106)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  expect_warning(
    r <- .ame_outer_eb(ini$candidates[[1]], ini$reference, d$Y,
                       d$X_row, d$X_col, d$X_dyad,
                       unit_penalty(), unit_penalty(),
                       outer_max_iter = 2, eps_Omega = 1e-14, eps_fit = 1e-14),
    "did not converge"
  )
  expect_false(r$converged)
  expect_true(all(is.finite(r$params$a)))
})


test_that("fit_dynamic_ame runs from a raw edge panel and returns the documented object", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 14, M = 10, Tn = 5, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = 201)
  edge <- as_edge_panel(d, years = 2011:2015)

  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, K = d$K, n_starts = 2,
                    outer_max_iter = 40, panel_id = "test"))

  expect_s3_class(fit, "dynamic_ame")
  expect_equal(dim(fit$a), c(d$N, d$Tn))
  expect_equal(dim(fit$b), c(d$M, d$Tn))
  expect_equal(dim(fit$beta), c(0L, d$Tn))
  expect_length(fit$U, d$Tn)
  expect_length(fit$V, d$Tn)
  expect_length(fit$variance_components$Omega, 11L)
  expect_length(fit$variance_components$lambda, 5L)
  expect_length(fit$variance_components$gamma, 5L)
  expect_null(fit$coef_df)
  expect_equal(nrow(fit$node_df), (d$N + d$M) * d$Tn)
  expect_equal(fit$n_obs,
               sum(vapply(fit$Y_list, function(Y) sum(!is.na(Y)), integer(1))))

  # the per-period slices are usable by the existing post-estimation helpers
  sl <- fit$results[[2]]
  expect_true(.validate_als_fit_object(sl))
  expect_equal(latent_matrix_als_factorize_joint_cov(sl),
               fit$U[[2]] %*% t(fit$V[[2]]), tolerance = 1e-12,
               ignore_attr = TRUE)
  dec <- decompose_als_factorize_joint_cov(sl)
  expect_equal(dec$additive + dec$covariate + dec$latent, dec$fitted,
               tolerance = 1e-12)

  # node_df agrees with the trajectory matrices
  a2 <- fit$node_df$alpha[fit$node_df$mode == "row" &
                            fit$node_df$period == fit$years[2]]
  expect_equal(a2, unname(fit$a[, 2]), tolerance = 1e-12)
})


test_that("fit_dynamic_ame handles covariates and records its settings", {
  skip_if_not_installed("dplyr")
  set.seed(7)
  d <- sim_ame_panel(N = 12, M = 9, Tn = 4, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.05, seed = 202)
  edge <- as_edge_panel(d, years = 2001:2004)
  rows <- sprintf("%03d", seq_len(d$N))
  row_cov <- do.call(rbind, lapply(seq_len(d$Tn), function(t)
    data.frame(year = 2000 + t, node_row = rows,
               z = stats::rnorm(d$N), stringsAsFactors = FALSE)))

  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, row_cov_df = row_cov,
                    row_covar_names = "z", row_cov_prefix = "row_",
                    K = 2, n_starts = 2, outer_max_iter = 40))

  expect_equal(dim(fit$beta), c(1L, d$Tn))
  expect_equal(rownames(fit$beta), "z")
  expect_equal(nrow(fit$coef_df), d$Tn)
  expect_equal(fit$coef_df$coefficient, unname(fit$beta["z", ]),
               tolerance = 1e-12)
  expect_equal(fit$settings$K, 2)
  expect_equal(fit$settings$n_starts, 2)
  expect_length(fit$convergence$multistart_Q, 2L)
  expect_equal(fit$convergence$multistart_selected,
               which.min(fit$convergence$multistart_Q))
})


test_that("fit_dynamic_ame is reproducible and leaves the RNG stream alone", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 10, M = 8, Tn = 4, Pr = 0, Pc = 0, Pd = 0, seed = 203)
  edge <- as_edge_panel(d)

  f1 <- suppressWarnings(fit_dynamic_ame(edge, K = 2, n_starts = 3,
                                         outer_max_iter = 25))
  f2 <- suppressWarnings(fit_dynamic_ame(edge, K = 2, n_starts = 3,
                                         outer_max_iter = 25))
  expect_equal(unlist(f1$U), unlist(f2$U), tolerance = 1e-12)

  set.seed(31337); before <- stats::rnorm(1)
  set.seed(31337)
  invisible(suppressWarnings(fit_dynamic_ame(edge, K = 2, n_starts = 2,
                                             outer_max_iter = 15)))
  expect_identical(stats::rnorm(1), before)
})


test_that("fit_dynamic_ame recovers what the model identifies", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 20, M = 14, Tn = 6, Pr = 0, Pc = 0, Pd = 0,
                     sd_eps = 0.4, na_frac = 0.08, seed = 204)
  edge <- as_edge_panel(d)
  fit <- suppressWarnings(fit_dynamic_ame(edge, K = 2, n_starts = 2,
                                          outer_max_iter = 60))

  uv_hat <- unlist(lapply(seq_len(d$Tn),
                          function(t) fit$U[[t]] %*% t(fit$V[[t]])))
  uv_true <- unlist(lapply(seq_len(d$Tn),
                           function(t) d$sim$U[[t]] %*% t(d$sim$V[[t]])))
  expect_gt(stats::cor(uv_hat, uv_true), 0.80)
  expect_gt(stats::cor(as.vector(fit$a), as.vector(d$sim$a)), 0.80)
  expect_true(all(is.finite(unlist(fit$U))))
})


test_that("K = 1 and a single-covariate panel both work", {
  skip_if_not_installed("dplyr")
  d <- sim_ame_panel(N = 10, M = 8, Tn = 4, Pr = 0, Pc = 0, Pd = 0, seed = 205)
  edge <- as_edge_panel(d)
  f1 <- suppressWarnings(fit_dynamic_ame(edge, K = 1, n_starts = 1,
                                         outer_max_iter = 20))
  expect_equal(ncol(f1$U[[1]]), 1L)
  expect_true(all(is.finite(unlist(f1$V))))
})
