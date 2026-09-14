# Gibbs sampler for the dynamic bipartite AME model.
#
# The load-bearing claim is that a sweep draws from the same conditional the
# block updates take the mean of. Two tests carry it: the sampling primitive is
# checked against a densely constructed Gaussian, and the posterior mean of
# beta is checked against the optimiser's point estimate on the same data. If
# the sampler targeted a different distribution, the second would fail even
# with the first passing.

# --- the sampling primitive --------------------------------------------

test_that("block draws have the mean and covariance of the dense Gaussian", {
  set.seed(11)
  Tn <- 4L; q <- 3L; off <- -0.7; scale <- 2.3

  # a symmetric positive-definite block-tridiagonal system
  D <- lapply(seq_len(Tn), function(t) {
    A <- matrix(stats::rnorm(q * q), q, q)
    crossprod(A) + (q + 2) * diag(q)
  })
  r <- lapply(seq_len(Tn), function(t) stats::rnorm(q))

  # the same matrix, densely
  H <- matrix(0, q * Tn, q * Tn)
  for (t in seq_len(Tn)) {
    ix <- ((t - 1L) * q + 1L):(t * q)
    H[ix, ix] <- D[[t]]
    if (t < Tn) {
      jx <- (t * q + 1L):((t + 1L) * q)
      H[ix, jx] <- off * diag(q)
      H[jx, ix] <- off * diag(q)
    }
  }
  Hinv <- solve(H)
  mu_ref <- as.vector(Hinv %*% unlist(r))

  B <- 20000L
  draws <- vapply(seq_len(B),
                  function(b) as.vector(.rtridiag_block(D, off, r, scale)),
                  numeric(q * Tn))

  # the empirical mean is itself a random quantity, so it is judged against
  # its own Monte Carlo error rather than a relative tolerance -- at this B a
  # fixed tolerance would be either vacuous or flaky depending on the scale
  se <- sqrt(diag(scale * Hinv) / B)
  expect_lt(max(abs(rowMeans(draws) - mu_ref) / se), 4)
  expect_equal(stats::cov(t(draws)), scale * Hinv,
               tolerance = 0.06, ignore_attr = TRUE)
})


test_that("the scalar path agrees with the block path at q = 1", {
  set.seed(12)
  Tn <- 5L; n <- 4L; off <- -0.45; scale <- 1.7

  D <- matrix(stats::runif(n * Tn, 2, 4), n, Tn)
  r <- matrix(stats::rnorm(n * Tn), n, Tn)

  # means must agree exactly: both call the same solver underneath
  set.seed(1); s_draw <- .rtridiag_scalar_vec(D, off, r, scale)
  expect_equal(dim(s_draw), c(n, Tn))

  mu_scalar <- .solve_scalar_tridiagonal_vec(D, off, r)
  for (i in seq_len(n)) {
    Di <- lapply(seq_len(Tn), function(t) matrix(D[i, t], 1, 1))
    ri <- lapply(seq_len(Tn), function(t) r[i, t])
    mu_block <- unlist(.solve_block_tridiagonal(Di, off, ri))
    expect_equal(mu_block, mu_scalar[i, ], tolerance = 1e-10,
                 ignore_attr = TRUE)
  }

  # and the dispersion must match the dense reference for one system
  i <- 2L
  Di <- lapply(seq_len(Tn), function(t) matrix(D[i, t], 1, 1))
  H <- diag(D[i, ])
  for (t in seq_len(Tn - 1L)) { H[t, t + 1L] <- off; H[t + 1L, t] <- off }

  B <- 20000L
  set.seed(99)
  dr <- vapply(seq_len(B), function(b)
    .rtridiag_scalar_vec(D, off, r, scale)[i, ], numeric(Tn))
  expect_equal(stats::cov(t(dr)), scale * solve(H),
               tolerance = 0.06, ignore_attr = TRUE)
})


test_that("a zero-variance scale collapses the draw onto the solve", {
  set.seed(13)
  Tn <- 3L; q <- 2L; off <- -0.3
  D <- lapply(seq_len(Tn), function(t) crossprod(matrix(stats::rnorm(4), 2, 2)) +
                3 * diag(q))
  r <- lapply(seq_len(Tn), function(t) stats::rnorm(q))
  expect_equal(.rtridiag_block(D, off, r, 0),
               matrix(unlist(.solve_block_tridiagonal(D, off, r)), q, Tn),
               tolerance = 1e-12)
})


# --- the sampler as a whole ---------------------------------------------

gibbs_fixture <- function(seed = 501, K = 2, Tn = 4L, N = 12L, M = 9L) {
  d <- sim_ame_panel(N = N, M = M, Tn = Tn, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0.1, seed = seed)
  edge <- as_edge_panel(d, years = 2001:(2000 + Tn))
  fit <- suppressWarnings(
    fit_dynamic_ame(edge_panel = edge, K = K, n_starts = 1,
                    outer_max_iter = 30))
  list(d = d, fit = fit)
}


test_that("the sampler runs and returns the documented structure", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()

  g <- gibbs_from_fit(s$fit, n_iter = 40, burn = 10, seed = 3)

  expect_s3_class(g, "gibbs_dynamic_ame")
  expect_equal(nrow(g$omega), 40L)
  expect_equal(ncol(g$omega), 11L)
  expect_equal(colnames(g$omega)[1], "sigma_eps2")
  expect_true(all(is.finite(g$omega[, "sigma_eps2"])))
  expect_true(all(g$omega[, "sigma_eps2"] > 0))
  expect_true(g$seconds > 0)
  expect_null(g$latent)

  # the fixture has no covariates, so there is no beta to store
  expect_null(g$beta)

  # the last state is a full parameter list, so a chain can be restarted
  expect_setequal(names(g$params), c("beta", "a", "b", "U", "V"))
  expect_equal(dim(g$params$a), c(nrow(s$fit$a), length(s$fit$years)))
  expect_equal(length(g$params$U), length(s$fit$years))
})


test_that("thinning and burn-in keep the right number of sweeps", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  g <- gibbs_from_fit(s$fit, n_iter = 30, burn = 5, thin = 3, seed = 4)
  expect_equal(nrow(g$omega), 10L)
  expect_equal(g$settings$n_kept, 10L)
})


test_that("results are reproducible from the seed", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  a <- gibbs_from_fit(s$fit, n_iter = 20, burn = 5, seed = 7)
  b <- gibbs_from_fit(s$fit, n_iter = 20, burn = 5, seed = 7)
  c2 <- gibbs_from_fit(s$fit, n_iter = 20, burn = 5, seed = 8)
  expect_equal(a$omega, b$omega, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(a$omega, c2$omega)))
})


test_that("unobserved cells never enter the sampler", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  # the fixture has missing cells; a draw that touched them would produce
  # NA residuals and hence NA variances
  expect_true(any(vapply(s$fit$Y_list, function(Y) any(is.na(Y)), logical(1))))
  g <- gibbs_from_fit(s$fit, n_iter = 25, burn = 5, seed = 9)
  expect_true(all(is.finite(g$omega[, c("sigma_eps2", "sigma_a2", "sigma_b2")])))
})


test_that("dimnames survive a sweep", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  g <- gibbs_from_fit(s$fit, n_iter = 6, burn = 2, seed = 10)
  expect_equal(rownames(g$params$a), s$fit$row_ids)
  expect_equal(rownames(g$params$b), s$fit$col_ids)
  expect_equal(rownames(g$params$U[[1]]), s$fit$row_ids)
  expect_equal(rownames(g$params$V[[1]]), s$fit$col_ids)
})


test_that("storing the latent products gives one panel per kept sweep", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  g <- gibbs_from_fit(s$fit, n_iter = 8, burn = 2, seed = 11,
                      store_latent = TRUE)
  expect_length(g$latent, 8L)
  expect_equal(dim(g$latent[[1]][[1]]),
               c(length(s$fit$row_ids), length(s$fit$col_ids)))
})


# --- the sampler and the optimiser target the same model ----------------

test_that("the posterior mean of beta lands near the point estimate", {
  skip_if_not_installed("dplyr")

  # a panel with one dyadic covariate, large enough that the posterior is
  # tight and close to Gaussian, so its mean and its mode should nearly
  # coincide. A sampler aimed at a different target would miss here even
  # though every primitive test above passed.
  Tn <- 4L; N <- 25L; M <- 20L
  d <- sim_ame_panel(N = N, M = M, Tn = Tn, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = 0, seed = 777)
  edge <- as_edge_panel(d, years = 2001:(2000 + Tn))

  set.seed(777)
  rows <- sort(unique(edge$node_row)); cols <- sort(unique(edge$node_col))
  dyad <- expand.grid(node_row = rows, node_col = cols,
                      year = 2001:(2000 + Tn), stringsAsFactors = FALSE)
  dyad$z <- stats::rnorm(nrow(dyad))
  # give the covariate a real coefficient so the comparison is not about zero
  edge <- merge(edge, dyad, by = c("node_row", "node_col", "year"))
  edge$value <- edge$value + 0.8 * edge$z

  fit <- suppressWarnings(fit_dynamic_ame(
    edge_panel = edge[, c("year", "node_row", "node_col", "value")],
    dyad_cov_df = edge[, c("year", "node_row", "node_col", "z")],
    dyad_covar_names = "z", row_prefix = "", col_prefix = "",
    dyad_cov_row_prefix = "", dyad_cov_col_prefix = "",
    K = 2, n_starts = 2))

  g <- gibbs_from_fit(fit, n_iter = 1500, burn = 500, seed = 21)

  post_mean <- colMeans(g$beta)
  map <- as.vector(fit$beta)
  post_sd <- apply(g$beta, 2, stats::sd)

  # within half a posterior standard deviation of the mode, period by period
  expect_lt(max(abs(post_mean - map) / post_sd), 0.5)

  # and both near the value the data were generated with
  expect_lt(abs(mean(post_mean) - 0.8), 0.15)
})


test_that("summary reports both effective sample sizes", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("posterior")
  s <- gibbs_fixture()
  g <- gibbs_from_fit(s$fit, n_iter = 200, burn = 50, seed = 12)
  sm <- summary(g, pars = "omega")
  expect_true(all(c("ess_bulk", "ess_tail") %in% names(sm)))
  expect_true(all(sm$ess_bulk[!is.na(sm$ess_bulk)] > 0))
  expect_true(all(sm$conf.low <= sm$conf.high))
})


test_that("R-hat needs several chains and rejects one", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("posterior")
  s <- gibbs_fixture()
  g1 <- gibbs_from_fit(s$fit, n_iter = 120, burn = 30, seed = 31)
  expect_error(gibbs_rhat(list(g1), pars = "omega"), "at least two chains")

  g2 <- gibbs_from_fit(s$fit, n_iter = 120, burn = 30, seed = 32)
  g3 <- gibbs_from_fit(s$fit, n_iter = 120, burn = 30, seed = 33)
  rh <- gibbs_rhat(list(g1, g2, g3), pars = "omega")
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(rh)))
  expect_true(all(is.finite(rh$rhat)))
})


test_that("R-hat works on unnamed draws and never returns nothing", {
  skip_if_not_installed("posterior")

  # A caller that takes max() over an empty result reads -Inf as a perfect
  # R-hat, so a check that computed nothing would look like convergence.
  # Unnamed columns used to produce exactly that.
  set.seed(41)
  mk <- function(mu) {
    Mx <- cbind(stats::rnorm(400, mu, 0.02), stats::rnorm(400, -0.6, 0.02))
    Mx                                       # deliberately no column names
  }
  rh <- gibbs_rhat(list(list(beta = mk(0.8)), list(beta = mk(0.8))),
                   pars = "beta")
  expect_equal(nrow(rh), 2L)
  expect_true(all(is.finite(rh$rhat)))
  expect_equal(rh$parameter, c("par1", "par2"))

  # all-NA draws are undefined, not converged
  na_ch <- list(list(beta = matrix(NA_real_, 100, 2)),
                list(beta = matrix(NA_real_, 100, 2)))
  expect_error(gibbs_rhat(na_ch, pars = "beta"), "would look like convergence")

  # chains of different width are a mismatch, not something to silently pool
  expect_error(
    gibbs_rhat(list(list(beta = mk(0.8)),
                    list(beta = cbind(mk(0.8), stats::rnorm(400)))),
               pars = "beta"),
    "same parameters")
})


test_that("R-hat separates chains that agree from one that does not", {
  skip_if_not_installed("posterior")

  # three chains at one value and one elsewhere: the pathology the timing
  # study's calibration has to detect, since a posterior with two modes
  # never brings the full set below the target however long it is run
  set.seed(42)
  mk <- function(mu) {
    Mx <- cbind(stats::rnorm(3000, mu, 0.02))
    colnames(Mx) <- "x1:1"; Mx
  }
  ch <- list(list(beta = mk(0.7614)), list(beta = mk(0.7623)),
             list(beta = mk(0.8047)), list(beta = mk(0.7624)))

  all_four <- gibbs_rhat(ch, pars = "beta")
  agreeing <- gibbs_rhat(ch[c(1, 2, 4)], pars = "beta")

  expect_gt(max(all_four$rhat), 1.05)
  expect_lt(max(agreeing$rhat), 1.01)
})


# ── the reciprocal scale gauge ───────────────────────────────────────────
#
# The likelihood sees only U_t V_t', so the pair can drift along
# (c U, V / c) forever with the product fixed. Nothing identified suffers,
# but sigma_U^2 and tau_U^2 are reported quantities, and left alone they come
# out describing where the chain drifted. The sampler therefore applies the
# same normalisation the estimator does, and these check that it changes what
# it should and nothing else.

test_that("normalising the scale leaves the latent product untouched", {
  set.seed(51)
  Tn <- 4L; N <- 7L; M <- 9L; K <- 2L
  U <- lapply(seq_len(Tn), function(t) matrix(stats::rnorm(N * K), N, K))
  V <- lapply(seq_len(Tn), function(t) matrix(stats::rnorm(M * K), M, K))
  # push the pair far off the convention, as a drifting chain would
  U <- lapply(U, function(x) 6 * x)
  V <- lapply(V, function(x) x / 6)

  before <- lapply(seq_len(Tn), function(t) U[[t]] %*% t(V[[t]]))
  # The convention the sampler pins itself to. The estimator's default has
  # since moved to the initial states; the sampler has not followed, so this
  # asks for what the sampler actually applies rather than for whatever the
  # default happens to be.
  sn <- .scale_normalize_UV(U, V, anchor = "pooled")
  after <- lapply(seq_len(Tn), function(t) sn$U[[t]] %*% t(sn$V[[t]]))

  expect_equal(after, before, tolerance = 1e-12)

  # and afterwards the two sides carry equal per-coordinate magnitude
  S_U <- sum(vapply(sn$U, function(x) sum(x^2), numeric(1))) / (N * K * Tn)
  S_V <- sum(vapply(sn$V, function(x) sum(x^2), numeric(1))) / (M * K * Tn)
  expect_equal(S_U, S_V, tolerance = 1e-10)
})


test_that("the sampler returns states in the normalised gauge", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  g <- gibbs_from_fit(s$fit, n_iter = 30, burn = 10, seed = 61)

  U <- g$params$U; V <- g$params$V
  N <- nrow(U[[1]]); M <- nrow(V[[1]]); K <- ncol(U[[1]]); Tn <- length(U)
  S_U <- sum(vapply(U, function(x) sum(x^2), numeric(1))) / (N * K * Tn)
  S_V <- sum(vapply(V, function(x) sum(x^2), numeric(1))) / (M * K * Tn)
  expect_equal(S_U, S_V, tolerance = 1e-8)

  expect_equal(rownames(U[[1]]), s$fit$row_ids)
  expect_equal(rownames(V[[1]]), s$fit$col_ids)
})


test_that("chains agree on the latent variances once the gauge is fixed", {
  skip_if_not_installed("dplyr")

  # This is the failure the normalisation exists for: four chains from
  # different seeds used to end at sigma_U^2 an order of magnitude apart
  # while their products agreed, because nothing pinned the split between
  # the two factors.
  s <- gibbs_fixture(seed = 505, Tn = 5L)
  ch <- lapply(1:3, function(k)
    gibbs_from_fit(s$fit, n_iter = 400, burn = 200, seed = 700 + k))

  keep <- 201:400
  su <- vapply(ch, function(g) mean(g$omega[keep, "sigma_U2"]), numeric(1))
  sv <- vapply(ch, function(g) mean(g$omega[keep, "sigma_V2"]), numeric(1))
  prod <- vapply(ch, function(g)
    mean(g$omega[keep, "sigma_U2"] * g$omega[keep, "sigma_V2"]), numeric(1))

  # the components now vary no more wildly than the product they multiply to,
  # which is the identified quantity and the natural yardstick
  spread <- function(x) diff(range(x)) / mean(x)
  expect_lt(spread(su), 3 * spread(prod) + 0.25)
  expect_lt(spread(sv), 3 * spread(prod) + 0.25)
})


test_that("input validation", {
  expect_error(gibbs_dynamic_ame(list()), "non-empty list")
  expect_error(gibbs_dynamic_ame(list(matrix(1, 2, 2)), n_iter = 0),
               "n_iter must be positive")
})


test_that("the latent innovation variance is shared, and drawn from pooled statistics", {
  d <- sim_ame_panel(N = 14, M = 10, K = 2, Tn = 5, Pr = 0, Pc = 0, Pd = 1,
                     sd_eps = 0.5, na_frac = 0, seed = 404)
  p <- random_params(d)
  pr <- list(a0 = 1e-3, b0 = 1e-3)

  set.seed(9)
  tied <- .gibbs_omega(p, d$Y, NULL, NULL, d$X_dyad, pr, tie_latent = TRUE)
  expect_identical(tied[["tau_U2"]], tied[["tau_V2"]])

  set.seed(9)
  sep <- .gibbs_omega(p, d$Y, NULL, NULL, d$X_dyad, pr, tie_latent = FALSE)
  expect_false(isTRUE(all.equal(sep[["tau_U2"]], sep[["tau_V2"]])))

  # Not an average of the two separate draws. Under a shared kappa the full
  # conditional is one inverse-gamma on the pooled sufficient statistics, and
  # averaging two draws from the separate conditionals would be a different
  # distribution -- the chain would no longer be sampling the posterior it
  # claims to. The two agree in order of magnitude and must not agree exactly.
  expect_false(isTRUE(all.equal(
    tied[["tau_U2"]], sqrt(sep[["tau_U2"]] * sep[["tau_V2"]]),
    tolerance = 1e-8)))
  expect_lt(abs(log(tied[["tau_U2"]] / sqrt(sep[["tau_U2"]] * sep[["tau_V2"]]))),
            log(3))

  # everything outside the latent block is untouched by the tie
  for (nm in c("sigma_eps2", "sigma_a2", "sigma_b2", "tau_a2", "tau_b2"))
    expect_identical(tied[[nm]], sep[[nm]])
})


test_that("the sampler's shared kappa reaches the sampler's own entry point", {
  skip_if_not_installed("dplyr")
  s <- gibbs_fixture()
  g <- gibbs_from_fit(s$fit, n_iter = 20, burn = 5, seed = 77)
  om <- g$omega
  # Every kept draw, not just the last: the tie is applied inside the sweep, so
  # a single draw agreeing could be coincidence.
  expect_true(all(om[, "tau_U2"] == om[, "tau_V2"]))
})
