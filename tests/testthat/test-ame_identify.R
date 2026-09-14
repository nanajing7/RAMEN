# Steps 7 and 8. These fix the two gauge freedoms the observation model leaves
# in the latent factors; they must not move the fit, and applying them twice
# must change nothing.
#
# Note the implementation applies the scale before the rotation, the reverse of
# the order written in the paper. Both orders give identical fitted values and
# identical variance components; only scale-then-rotate is idempotent, because
# rescaling changes the relative weight of U and V inside the Procrustes
# objective.

setup_fit <- function(seed = 55) {
  d <- sim_ame_panel(N = 11, M = 7, Tn = 4, Pr = 2, Pc = 0, Pd = 0, seed = seed)
  ini <- .init_ame_trajectory(d$Y, K = d$K, P = d$P, n_starts = 1)
  lam <- c(beta = .4, a = .5, b = .6, U = .7, V = .8)
  gam <- c(beta = 1.1, a = 1.2, b = 1.3, U = 1.4, V = 1.5)
  fit <- suppressWarnings(
    .ame_inner_bcd(ini$candidates[[1]], d$Y, d$X_row, NULL, NULL, lam, gam,
                   max_iter = 200))
  list(d = d, ini = ini, p = fit$params, lam = lam, gam = gam)
}

stacked_fit <- function(p, d) {
  unlist(lapply(seq_len(d$Tn), function(t)
    .fitted_period(p$a[, t], p$b[, t], p$U[[t]], p$V[[t]], p$beta[, t],
                   d$X_row[[t]])))
}


test_that("identification leaves the fitted values exactly unchanged", {
  s <- setup_fit()
  id <- .ame_identify(s$p, s$ini$reference)

  for (t in seq_len(s$d$Tn)) {
    expect_equal(id$params$U[[t]] %*% t(id$params$V[[t]]),
                 s$p$U[[t]] %*% t(s$p$V[[t]]), tolerance = 1e-10)
  }
  expect_equal(stacked_fit(id$params, s$d), stacked_fit(s$p, s$d),
               tolerance = 1e-10)
  expect_identical(id$params$a, s$p$a)
  expect_identical(id$params$b, s$p$b)
  expect_identical(id$params$beta, s$p$beta)
})


test_that("the rotation alone leaves Q and every penalty term intact", {
  # the spec claims Q-invariance for the rotation; the rescaling moves Q by
  # design, and runs between inner loops where Q is not compared
  s <- setup_fit()
  sn <- .scale_normalize_UV(s$p$U, s$p$V)
  p8 <- s$p; p8$U <- sn$U; p8$V <- sn$V
  R <- .orthogonal_procrustes(c(p8$U, p8$V),
                              c(s$ini$reference$U, s$ini$reference$V))
  p87 <- p8
  p87$U <- lapply(p8$U, function(U) U %*% R)
  p87$V <- lapply(p8$V, function(V) V %*% R)

  o8  <- .ame_objective(p8,  s$d$Y, s$d$X_row, NULL, NULL, s$lam, s$gam)
  o87 <- .ame_objective(p87, s$d$Y, s$d$X_row, NULL, NULL, s$lam, s$gam)
  expect_equal(o87$Q, o8$Q, tolerance = 1e-9)
  expect_equal(o87$penalty_initial, o8$penalty_initial, tolerance = 1e-10)
  expect_equal(o87$penalty_rw, o8$penalty_rw, tolerance = 1e-10)
})


test_that("the result satisfies both gauges at once, with N != M", {
  s <- setup_fit()
  K <- s$d$K

  # The rotation is orthogonal, and the scale convention holds wherever it was
  # anchored. Each is checked where it makes its claim and nowhere else:
  # equalising the increments says nothing about the magnitudes, and vice
  # versa, which is what makes them different conventions rather than the same
  # one written twice.
  Tn <- s$d$Tn
  inc <- function(L, n) sum(vapply(2:Tn, function(t)
    sum((L[[t]] - L[[t - 1L]])^2), numeric(1))) / (n * K * (Tn - 1L))

  idi <- .ame_identify(s$p, s$ini$reference, anchor = "innovation")
  expect_equal(crossprod(idi$R), diag(K), tolerance = 1e-12)
  expect_equal(inc(idi$params$U, s$d$N), inc(idi$params$V, s$d$M),
               tolerance = 1e-10)

  idf <- .ame_identify(s$p, s$ini$reference, anchor = "first")
  expect_equal(sum(idf$params$U[[1]]^2) / (s$d$N * K),
               sum(idf$params$V[[1]]^2) / (s$d$M * K), tolerance = 1e-10)

  idp <- .ame_identify(s$p, s$ini$reference, anchor = "pooled")
  expect_equal(crossprod(idp$R), diag(K), tolerance = 1e-12)
  expect_equal(
    sum(vapply(idp$params$U, function(U) sum(U^2), numeric(1))) /
      (s$d$N * K * s$d$Tn),
    sum(vapply(idp$params$V, function(V) sum(V^2), numeric(1))) /
      (s$d$M * K * s$d$Tn),
    tolerance = 1e-10)
})


test_that("identification is idempotent", {
  s <- setup_fit()
  id <- .ame_identify(s$p, s$ini$reference)
  id2 <- .ame_identify(id$params, s$ini$reference)

  expect_equal(id2$R, diag(s$d$K), tolerance = 1e-9)
  expect_equal(id2$c, 1, tolerance = 1e-10)
  expect_equal(unlist(id2$params$U), unlist(id$params$U), tolerance = 1e-10)

  id3 <- .ame_identify(id2$params, s$ini$reference)
  expect_equal(unlist(id3$params$U), unlist(id2$params$U), tolerance = 1e-10)
})


test_that("a scrambled solution is restored exactly", {
  s <- setup_fit()
  id <- .ame_identify(s$p, s$ini$reference)
  K <- s$d$K

  set.seed(2)
  Rt <- qr.Q(qr(matrix(stats::rnorm(K * K), K, K)))
  scr <- id$params
  scr$U <- lapply(scr$U, function(U) U %*% Rt)
  scr$V <- lapply(scr$V, function(V) V %*% Rt)
  back <- .ame_identify(scr, s$ini$reference)
  expect_equal(unlist(back$params$U), unlist(id$params$U), tolerance = 1e-9)
  expect_equal(back$R, t(Rt), tolerance = 1e-9)

  bad <- id$params
  bad$U <- lapply(bad$U, function(U) 40 * U)
  bad$V <- lapply(bad$V, function(V) V / 40)
  fixed <- .ame_identify(bad, s$ini$reference)
  expect_equal(fixed$c, 1 / 40, tolerance = 1e-10)
  expect_equal(unlist(fixed$params$U), unlist(id$params$U), tolerance = 1e-9)
})


test_that("identification stops the variance drift it exists to prevent", {
  s <- setup_fit()
  id <- .ame_identify(s$p, s$ini$reference)

  no_fix <- yes_fix <- id$params
  drift_no <- drift_yes <- numeric(5)
  for (m in 1:5) {
    no_fix$U <- lapply(no_fix$U, function(U) 1.6 * U)
    no_fix$V <- lapply(no_fix$V, function(V) V / 1.6)
    drift_no[m] <- .eb_update_variances(no_fix, s$d$Y, s$d$X_row, NULL, NULL,
                                        lambda = s$lam,
                                        gamma = s$gam)$Omega[["sigma_U2"]]

    yes_fix$U <- lapply(yes_fix$U, function(U) 1.6 * U)
    yes_fix$V <- lapply(yes_fix$V, function(V) V / 1.6)
    yes_fix <- .ame_identify(yes_fix, s$ini$reference)$params
    drift_yes[m] <- .eb_update_variances(yes_fix, s$d$Y, s$d$X_row, NULL, NULL,
                                         lambda = s$lam,
                                         gamma = s$gam)$Omega[["sigma_U2"]]
  }
  expect_gt(drift_no[5] / drift_no[1], 20)
  expect_equal(max(abs(drift_yes - drift_yes[1])) / drift_yes[1], 0,
               tolerance = 1e-8)
})


test_that("Omega is invariant to the choice of rotation", {
  s <- setup_fit()
  id <- .ame_identify(s$p, s$ini$reference)
  K <- s$d$K
  set.seed(3)
  Rr <- qr.Q(qr(matrix(stats::rnorm(K * K), K, K)))
  rot <- id$params
  rot$U <- lapply(rot$U, function(U) U %*% Rr)
  rot$V <- lapply(rot$V, function(V) V %*% Rr)

  o1 <- .eb_update_variances(id$params, s$d$Y, s$d$X_row, NULL, NULL,
                             lambda = s$lam, gamma = s$gam)$Omega
  o2 <- .eb_update_variances(rot, s$d$Y, s$d$X_row, NULL, NULL,
                             lambda = s$lam, gamma = s$gam)$Omega
  expect_equal(o1, o2, tolerance = 1e-9)
})


test_that("a mismatched reference is rejected", {
  s <- setup_fit()
  expect_error(
    .ame_identify(s$p, list(U = s$ini$reference$U[1], V = s$ini$reference$V[1])),
    "same number of periods"
  )
})
