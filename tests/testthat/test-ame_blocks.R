# The five trajectory block updates (Steps 2-6).
#
# Each update claims to be the exact minimiser of Q over its block. That is
# checked three independent ways: against a dense system assembled by hand from
# the spec's first-order conditions, by confirming random perturbations cannot
# lower Q, and by a numerical gradient.

pen_lam <- c(beta = .35, a = .45, b = .55, U = .65, V = .75)
pen_gam <- c(beta = 1.15, a = 1.20, b = 1.30, U = 1.45, V = 1.55)

Q_of <- function(p, d) {
  .ame_objective(p, d$Y, d$X_row, d$X_col, d$X_dyad, pen_lam, pen_gam)$Q
}

expect_exact_minimiser <- function(p_new, d, perturb, coords, n_draw = 100) {
  Q1 <- Q_of(p_new, d)

  # no random perturbation may improve on it
  beaten <- 0L
  for (s in seq_len(n_draw)) {
    if (Q_of(perturb(p_new), d) < Q1 - 1e-9) beaten <- beaten + 1L
  }
  expect_equal(beaten, 0L)

  # and the numerical gradient must vanish
  h <- 1e-6
  g <- vapply(coords, function(set) {
    up <- set(p_new, h); dn <- set(p_new, -h)
    (Q_of(up, d) - Q_of(dn, d)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(g)), 1e-4)
}


test_that("the additive diagonal follows the spec's boundary pattern", {
  n <- matrix(3, 2, 5)
  D <- .additive_diagonal(n, lambda = 10, gamma = 100)
  expect_equal(D[1, 1], 3 + 10 + 100)      # t = 1:  n + lambda + gamma
  expect_true(all(D[1, 2:4] == 3 + 200))   # interior: n + 2 gamma
  expect_equal(D[1, 5], 3 + 100)           # t = T:  n + gamma
  # a single period has no temporal neighbour at all
  expect_equal(.additive_diagonal(matrix(3, 1, 1), 10, 100)[1, 1], 13)
})


test_that("the a-update solves the spec's tridiagonal system exactly", {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 5, seed = 11)
  p <- random_params(d)
  p1 <- .update_a_block(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                        pen_lam[["a"]], pen_gam[["a"]])

  # rebuild node 1's system straight from the first-order conditions
  i <- 1L
  S <- numeric(d$Tn); nv <- numeric(d$Tn)
  for (t in seq_len(d$Tn)) {
    R <- d$Y[[t]] - matrix(p$b[, t], d$N, d$M, byrow = TRUE) -
      .cov_term(p$beta[, t], d$X_row[[t]], d$X_col[[t]], d$X_dyad[[t]], d$N, d$M) -
      p$U[[t]] %*% t(p$V[[t]])
    S[t] <- sum(R[i, ], na.rm = TRUE)
    nv[t] <- sum(!is.na(R[i, ]))
  }
  A <- matrix(0, d$Tn, d$Tn)
  for (t in seq_len(d$Tn)) {
    nb <- (t > 1) + (t < d$Tn)
    A[t, t] <- nv[t] + pen_gam[["a"]] * nb + if (t == 1) pen_lam[["a"]] else 0
    if (t < d$Tn) {
      A[t, t + 1] <- -pen_gam[["a"]]
      A[t + 1, t] <- -pen_gam[["a"]]
    }
  }
  expect_equal(as.numeric(solve(A, S)), p1$a[i, ], tolerance = 1e-10)

  expect_lt(Q_of(p1, d), Q_of(p, d))
  expect_exact_minimiser(
    p1, d,
    perturb = function(x) { x$a <- x$a + matrix(stats::rnorm(d$N * d$Tn, sd = .05),
                                                d$N, d$Tn); x },
    coords = list(function(x, h) { x$a[1, 1] <- x$a[1, 1] + h; x },
                  function(x, h) { x$a[5, 3] <- x$a[5, 3] + h; x },
                  function(x, h) { x$a[d$N, d$Tn] <- x$a[d$N, d$Tn] + h; x })
  )
})


test_that("the b-update is the exact minimiser and mirrors the a-update", {
  d <- sim_ame_panel(N = 10, M = 7, Tn = 5, seed = 12)
  p <- random_params(d)
  p2 <- .update_b_block(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                        pen_lam[["b"]], pen_gam[["b"]])
  expect_lt(Q_of(p2, d), Q_of(p, d))
  expect_exact_minimiser(
    p2, d,
    perturb = function(x) { x$b <- x$b + matrix(stats::rnorm(d$M * d$Tn, sd = .05),
                                                d$M, d$Tn); x },
    coords = list(function(x, h) { x$b[1, 1] <- x$b[1, 1] + h; x },
                  function(x, h) { x$b[d$M, d$Tn] <- x$b[d$M, d$Tn] + h; x })
  )
})


test_that("a single period reduces the additive update to a ridge mean", {
  d <- sim_ame_panel(Tn = 4, seed = 13)
  p <- random_params(d)
  p1 <- list(a = p$a[, 1, drop = FALSE], b = p$b[, 1, drop = FALSE],
             beta = p$beta[, 1, drop = FALSE], U = p$U[1], V = p$V[1])
  r <- .update_a_block(p1, d$Y[1], d$X_row[1], d$X_col[1], d$X_dyad[1],
                       pen_lam[["a"]], pen_gam[["a"]])
  R <- d$Y[[1]] - matrix(p1$b[, 1], d$N, d$M, byrow = TRUE) -
    .cov_term(p1$beta[, 1], d$X_row[[1]], d$X_col[[1]], d$X_dyad[[1]], d$N, d$M) -
    p1$U[[1]] %*% t(p1$V[[1]])
  expect_equal(unname(r$a[, 1]),
               unname(rowSums(R, na.rm = TRUE) /
                        (rowSums(!is.na(R)) + pen_lam[["a"]])),
               tolerance = 1e-12)
})


test_that("the U-update solves every node's block system exactly", {
  d <- sim_ame_panel(N = 9, M = 6, Tn = 4, seed = 21)
  p <- random_params(d)
  pU <- .update_U_block(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                        pen_lam[["U"]], pen_gam[["U"]])
  K <- d$K

  partial <- function(t) {
    d$Y[[t]] - matrix(p$a[, t], d$N, d$M) -
      matrix(p$b[, t], d$N, d$M, byrow = TRUE) -
      .cov_term(p$beta[, t], d$X_row[[t]], d$X_col[[t]], d$X_dyad[[t]], d$N, d$M)
  }

  for (i in seq_len(d$N)) {
    G <- vector("list", d$Tn); cc <- vector("list", d$Tn)
    for (t in seq_len(d$Tn)) {
      R <- partial(t)
      Gt <- matrix(0, K, K); ct <- numeric(K)
      for (j in seq_len(d$M)) {
        if (is.na(R[i, j])) next
        v <- p$V[[t]][j, ]
        Gt <- Gt + outer(v, v)
        ct <- ct + v * R[i, j]
      }
      G[[t]] <- Gt; cc[[t]] <- ct
    }
    A <- matrix(0, d$Tn * K, d$Tn * K); rr <- numeric(d$Tn * K)
    for (t in seq_len(d$Tn)) {
      idx <- ((t - 1) * K + 1):(t * K)
      nb <- (t > 1) + (t < d$Tn)
      A[idx, idx] <- G[[t]] +
        (pen_gam[["U"]] * nb + if (t == 1) pen_lam[["U"]] else 0) * diag(K)
      rr[idx] <- cc[[t]]
      if (t < d$Tn) {
        nx <- (t * K + 1):((t + 1) * K)
        A[idx, nx] <- -pen_gam[["U"]] * diag(K)
        A[nx, idx] <- -pen_gam[["U"]] * diag(K)
      }
    }
    # the dense solution is stacked period-major
    got <- as.vector(vapply(seq_len(d$Tn), function(t) pU$U[[t]][i, ],
                            numeric(K)))
    expect_equal(as.numeric(solve(A, rr)), got, tolerance = 1e-9)
  }

  expect_lt(Q_of(pU, d), Q_of(p, d))
  expect_exact_minimiser(
    pU, d,
    perturb = function(x) { x$U <- lapply(x$U, function(U)
      U + matrix(stats::rnorm(d$N * K, sd = .05), d$N, K)); x },
    coords = list(
      function(x, h) { x$U[[1]][1, 1] <- x$U[[1]][1, 1] + h; x },
      function(x, h) { x$U[[d$Tn]][d$N, K] <- x$U[[d$Tn]][d$N, K] + h; x })
  )
})


test_that("the V-update is the exact minimiser and equals U on the transpose", {
  d <- sim_ame_panel(N = 9, M = 6, Tn = 4, Pr = 0, Pc = 0, Pd = 0, seed = 22)
  p <- random_params(d)
  pV <- .update_V_block(p, d$Y, NULL, NULL, NULL, pen_lam[["V"]], pen_gam[["V"]])
  expect_lt(Q_of(pV, d), Q_of(p, d))

  # updating V on Y is updating U on t(Y) with the two node sets swapped
  pt <- list(a = p$b, b = p$a, beta = p$beta, U = p$V, V = p$U)
  uu <- .update_U_block(pt, lapply(d$Y, t), NULL, NULL, NULL,
                        pen_lam[["V"]], pen_gam[["V"]])$U
  expect_equal(unlist(pV$V), unlist(uu), tolerance = 1e-10)
})


test_that("the beta-update solves the single stacked system exactly", {
  d <- sim_ame_panel(N = 8, M = 6, Tn = 4, Pr = 2, Pc = 1, Pd = 2, seed = 31)
  p <- random_params(d)
  pB <- .update_beta_block(p, d$Y, d$X_row, d$X_col, d$X_dyad,
                           pen_lam[["beta"]], pen_gam[["beta"]])
  P <- d$P

  S <- vector("list", d$Tn); rr <- vector("list", d$Tn)
  for (t in seq_len(d$Tn)) {
    St <- matrix(0, P, P); rt <- numeric(P)
    for (i in seq_len(d$N)) {
      for (j in seq_len(d$M)) {
        y <- d$Y[[t]][i, j]
        if (is.na(y)) next
        x <- c(d$X_row[[t]][i, ], d$X_col[[t]][j, ], d$X_dyad[[t]][i, j, ])
        resid <- y - p$a[i, t] - p$b[j, t] -
          sum(p$U[[t]][i, ] * p$V[[t]][j, ])
        St <- St + outer(x, x)
        rt <- rt + x * resid
      }
    }
    S[[t]] <- St; rr[[t]] <- rt
  }
  A <- matrix(0, d$Tn * P, d$Tn * P); bb <- numeric(d$Tn * P)
  for (t in seq_len(d$Tn)) {
    idx <- ((t - 1) * P + 1):(t * P)
    nb <- (t > 1) + (t < d$Tn)
    A[idx, idx] <- S[[t]] +
      (pen_gam[["beta"]] * nb + if (t == 1) pen_lam[["beta"]] else 0) * diag(P)
    bb[idx] <- rr[[t]]
    if (t < d$Tn) {
      nx <- (t * P + 1):((t + 1) * P)
      A[idx, nx] <- -pen_gam[["beta"]] * diag(P)
      A[nx, idx] <- -pen_gam[["beta"]] * diag(P)
    }
  }
  expect_equal(as.numeric(solve(A, bb)), as.vector(pB$beta), tolerance = 1e-9)

  expect_lt(Q_of(pB, d), Q_of(p, d))
  expect_exact_minimiser(
    pB, d,
    perturb = function(x) { x$beta <- x$beta +
      matrix(stats::rnorm(P * d$Tn, sd = .05), P, d$Tn); x },
    coords = list(function(x, h) { x$beta[1, 1] <- x$beta[1, 1] + h; x },
                  function(x, h) { x$beta[P, d$Tn] <- x$beta[P, d$Tn] + h; x })
  )
})


test_that("the beta-update is a no-op when there are no covariates", {
  d <- sim_ame_panel(Pr = 0, Pc = 0, Pd = 0, seed = 32)
  p <- random_params(d)
  expect_identical(
    .update_beta_block(p, d$Y, NULL, NULL, NULL, pen_lam[["beta"]],
                       pen_gam[["beta"]]),
    p
  )
})


test_that("a full sweep in the spec's order never increases Q", {
  d <- sim_ame_panel(N = 9, M = 6, Tn = 4, seed = 41)
  p <- random_params(d)
  q <- Q_of(p, d)
  for (it in 1:6) {
    for (upd in list(
      function(x) .update_beta_block(x, d$Y, d$X_row, d$X_col, d$X_dyad,
                                     pen_lam[["beta"]], pen_gam[["beta"]]),
      function(x) .update_a_block(x, d$Y, d$X_row, d$X_col, d$X_dyad,
                                  pen_lam[["a"]], pen_gam[["a"]]),
      function(x) .update_b_block(x, d$Y, d$X_row, d$X_col, d$X_dyad,
                                  pen_lam[["b"]], pen_gam[["b"]]),
      function(x) .update_U_block(x, d$Y, d$X_row, d$X_col, d$X_dyad,
                                  pen_lam[["U"]], pen_gam[["U"]]),
      function(x) .update_V_block(x, d$Y, d$X_row, d$X_col, d$X_dyad,
                                  pen_lam[["V"]], pen_gam[["V"]]))) {
      p <- upd(p)
      q_new <- Q_of(p, d)
      expect_lte(q_new, q + 1e-9)
      q <- q_new
    }
  }
})


test_that("block updates survive missing covariates and empty periods", {
  d <- sim_ame_panel(N = 9, M = 7, Tn = 4, seed = 51)
  p <- random_params(d)

  Xr <- d$X_row
  Xr[[2]][3, 1] <- NA                      # node 3 loses a covariate in period 2
  expect_true(all(is.finite(
    .update_beta_block(p, d$Y, Xr, d$X_col, d$X_dyad,
                       pen_lam[["beta"]], pen_gam[["beta"]])$beta)))

  Y2 <- d$Y
  Y2[[2]][4, ] <- NA                       # node 4 unobserved throughout period 2
  expect_true(all(is.finite(
    .update_a_block(p, Y2, d$X_row, d$X_col, d$X_dyad,
                    pen_lam[["a"]], pen_gam[["a"]])$a)))
  expect_true(all(is.finite(unlist(
    .update_U_block(p, Y2, d$X_row, d$X_col, d$X_dyad,
                    pen_lam[["U"]], pen_gam[["U"]])$U))))
})
