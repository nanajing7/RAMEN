# The acceleration is meant to arrive at the same place sooner. Everything
# here is about the first half of that sentence: the plain iteration defines
# the estimate, so if the two ever disagree it is the acceleration that is
# wrong, and no speed-up would make up for it.

# The panel is deliberately large enough, and quiet enough, for every variance
# component to be identified. On a smaller one some component collapses toward
# zero, and then there is no fixed point for the two modes to agree on: both
# are still descending when the convergence criterion -- a relative change in
# the norm of Omega, which a small component cannot move -- halts them, and
# they halt at different depths. That behaviour is real and is documented with
# the accelerator, but it is not what these tests are about.
fixture <- function(seed = 1, Tn = 6) {
  d <- sim_ame_panel(N = 25, M = 20, K = 2, Tn = Tn, Pr = 0, Pc = 0, Pd = 2,
                     sd_eps = 0.3, na_frac = 0.05, seed = seed)
  init <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1L, seed = 7)
  one <- unit_penalty()
  start <- suppressWarnings(
    .ame_inner_bcd(init$candidates[[1]], d$Y, NULL, NULL, d$X_dyad,
                   one, one, max_iter = 200, eps_Q = 1e-5)$params
  )
  list(d = d, reference = init$reference, start = start, one = one)
}

run_outer <- function(fx, mode) {
  suppressWarnings(.ame_outer_eb(
    params = fx$start, reference = fx$reference, Y_list = fx$d$Y,
    X_row_list = NULL, X_col_list = NULL, X_dyad_list = fx$d$X_dyad,
    lambda = fx$one, gamma = fx$one,
    outer_max_iter = 300, inner_max_iter = 500,
    eps_Q = 1e-6, eps_Omega = 1e-4, eps_fit = 1e-4,
    accelerate = mode
  ))
}


test_that("acceleration reaches the same fixed point", {
  skip_if_not_installed("SQUAREM")
  fx <- fixture()
  plain <- run_outer(fx, "none")
  sq    <- run_outer(fx, "squarem")

  # Both stop on the same criterion, a relative change of 1e-4 in Omega and in
  # the fitted values, so agreement is expected at that order and not beyond.
  expect_equal(sq$Omega, plain$Omega, tolerance = 1e-3)
  expect_equal(sq$lambda, plain$lambda, tolerance = 1e-3)
  expect_equal(sq$gamma, plain$gamma, tolerance = 1e-3)
  expect_equal(sq$sigma_eps2, plain$sigma_eps2, tolerance = 1e-3)
})


test_that("acceleration reaches the same parameters, not merely the same variances", {
  skip_if_not_installed("SQUAREM")
  fx <- fixture()
  plain <- run_outer(fx, "none")
  sq    <- run_outer(fx, "squarem")

  # Omega agreeing while the trajectories differ would mean the two runs found
  # different points that happen to imply the same variances, which is exactly
  # the failure a variance-only check would miss.
  expect_equal(sq$params$beta, plain$params$beta, tolerance = 1e-2)
  expect_equal(sq$params$a, plain$params$a, tolerance = 1e-2)
  expect_equal(sq$params$b, plain$params$b, tolerance = 1e-2)

  # U and V are identified only up to a rotation, so they are compared through
  # the product, which is not.
  for (t in seq_along(sq$params$U)) {
    expect_equal(sq$params$U[[t]] %*% t(sq$params$V[[t]]),
                 plain$params$U[[t]] %*% t(plain$params$V[[t]]),
                 tolerance = 1e-2)
  }
})


test_that("acceleration certifies convergence by the original criterion", {
  skip_if_not_installed("SQUAREM")
  fx <- fixture()
  sq <- run_outer(fx, "squarem")

  # The extrapolation stops on its own tolerance; what makes `converged` mean
  # the same thing in both modes is that the plain iteration runs on from
  # there until both original criteria hold.
  expect_true(sq$converged)
  expect_identical(sq$accelerate, "squarem")
  expect_identical(run_outer(fx, "none")$accelerate, "none")
})


test_that("the recorded pass count is the work actually done", {
  skip_if_not_installed("SQUAREM")
  fx <- fixture()
  sq <- run_outer(fx, "squarem")

  # An accelerated step costs about three passes, and reporting steps rather
  # than passes would inflate the speed-up threefold. Every trace grows once
  # per pass, so their lengths are the check.
  expect_length(sq$objective_trace, sq$iterations)
  expect_length(sq$inner_iterations, sq$iterations)
  expect_length(sq$Omega_trace, sq$iterations)
})


test_that("acceleration never converges where the plain iteration would not", {
  skip_if_not_installed("SQUAREM")
  fx <- fixture()

  # Regression test. The certification phase used to be given whatever the
  # extrapolation left of a shared budget, and the extrapolation spends passes
  # backtracking; on two of eight real panels that left too few to finish, and
  # the accelerated fit reported itself unconverged where the plain one had
  # converged. A budget only a little larger than the plain iteration needs
  # reproduces that: under the old split the certification would get nothing.
  need <- run_outer(fx, "none")$iterations
  tight <- function(mode) suppressWarnings(.ame_outer_eb(
    params = fx$start, reference = fx$reference, Y_list = fx$d$Y,
    X_row_list = NULL, X_col_list = NULL, X_dyad_list = fx$d$X_dyad,
    lambda = fx$one, gamma = fx$one,
    outer_max_iter = need + 2L, inner_max_iter = 500,
    eps_Q = 1e-6, eps_Omega = 1e-4, eps_fit = 1e-4, accelerate = mode))

  plain <- tight("none")
  sq <- tight("squarem")
  expect_true(plain$converged)
  expect_true(sq$converged)
})


test_that("a collapsing component has no fixed point for either mode to reach", {
  skip_if_not_installed("SQUAREM")
  # Kept as a record of a real limitation rather than to check a behaviour: on
  # a panel too small to identify the innovation variance of `a`, tightening
  # the criterion keeps moving the plain iteration's estimate, which is what a
  # run that has not settled looks like. Neither mode is at fault, and the two
  # will not agree on that component; what protects a user is that the fit
  # reports its capped penalties, not that the modes agree.
  d <- sim_ame_panel(N = 12, M = 9, K = 2, Tn = 4, Pr = 0, Pc = 0, Pd = 2,
                     sd_eps = 0.6, na_frac = 0.05, seed = 1)
  init <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1L, seed = 7)
  one <- unit_penalty()
  start <- suppressWarnings(
    .ame_inner_bcd(init$candidates[[1]], d$Y, NULL, NULL, d$X_dyad,
                   one, one, max_iter = 200, eps_Q = 1e-5)$params)
  go <- function(eps) suppressWarnings(.ame_outer_eb(
    params = start, reference = init$reference, Y_list = d$Y,
    X_dyad_list = d$X_dyad, lambda = one, gamma = one,
    outer_max_iter = 400, inner_max_iter = 500, eps_Q = 1e-6,
    eps_Omega = eps, eps_fit = eps, accelerate = "none"))

  loose <- go(1e-4)
  tight <- go(1e-5)
  expect_true(loose$converged)
  # Both runs called themselves converged, and the estimate still fell by a
  # factor of three between them: the criterion stopped the first, it did not
  # certify it.
  expect_lt(tight$Omega[["tau_a2"]], loose$Omega[["tau_a2"]] / 2)
  expect_lt(tight$Omega[["tau_a2"]], 0.01)   # true value is 0.1225
})


test_that("acceleration takes no more passes than the plain iteration", {
  skip_if_not_installed("SQUAREM")
  fx <- fixture(seed = 3, Tn = 8)
  plain <- run_outer(fx, "none")
  sq    <- run_outer(fx, "squarem")

  # On a panel the plain loop already clears in a few passes there is nothing
  # to accelerate, and the extrapolation's own overhead can cost more than it
  # saves. That is a property of the problem, not a defect, so the comparison
  # is only made where it means something.
  skip_if(plain$iterations < 20,
          sprintf("the plain iteration converged in %d passes here; too few for the comparison to say anything",
                  plain$iterations))
  expect_lte(sq$iterations, plain$iterations)
})


test_that("fit_dynamic_ame carries the mode into the fit and the bootstrap", {
  d <- sim_ame_panel(N = 10, M = 8, K = 2, Tn = 3, Pr = 0, Pc = 0, Pd = 1,
                     sd_eps = 0.5, na_frac = 0, seed = 5)
  ep <- as_edge_panel(d)
  dyad <- do.call(rbind, lapply(seq_len(d$Tn), function(t) {
    g <- expand.grid(i = seq_len(d$N), j = seq_len(d$M))
    data.frame(year = t,
               node_row = sprintf("%03d", g$i), node_col = sprintf("%03d", g$j),
               x1 = d$X_dyad[[t]][cbind(g$i, g$j, 1)],
               stringsAsFactors = FALSE)
  }))

  fit <- suppressWarnings(fit_dynamic_ame(
    edge_panel = ep, dyad_cov_df = dyad, dyad_covar_names = "x1",
    row_prefix = "", col_prefix = "",
    dyad_cov_row_prefix = "", dyad_cov_col_prefix = "",
    K = 2, n_starts = 1, accelerate = "none"))

  # The bootstrap re-estimates every replicate the way the point estimate was
  # estimated, so the mode has to travel with the fit rather than being asked
  # for again at bootstrap time.
  expect_identical(fit$settings$accelerate, "none")

  # A fit made before this argument existed carries no such setting, and must
  # still refit rather than failing on a NULL.
  old <- fit
  old$settings$accelerate <- NULL
  expect_no_error(
    suppressWarnings(bootstrap_ame(old, design = "conditional", B = 2,
                                   n_cores = 1, seed = 1))
  )
})
