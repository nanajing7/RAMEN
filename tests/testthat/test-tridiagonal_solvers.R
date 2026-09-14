# Numerical primitives. Every routine here is checked against a dense
# computation of the same quantity, so a regression cannot hide behind an
# internally consistent but wrong recursion.

dense_H <- function(D, off) {
  Tn <- length(D)
  q <- if (is.matrix(D[[1]])) nrow(D[[1]]) else 1L
  H <- matrix(0, Tn * q, Tn * q)
  for (t in seq_len(Tn)) {
    idx <- ((t - 1) * q + 1):(t * q)
    H[idx, idx] <- if (is.matrix(D[[t]])) D[[t]] else matrix(D[[t]], 1, 1)
    if (t < Tn) {
      nx <- (t * q + 1):((t + 1) * q)
      H[idx, nx] <- off * diag(q)
      H[nx, idx] <- off * diag(q)
    }
  }
  H
}

mk_D <- function(q, off, Tn) {
  lapply(seq_len(Tn), function(t) {
    A <- crossprod(matrix(stats::rnorm(q * q), q, q))
    A + (2 * abs(off) + q + 1) * diag(q)
  })
}


test_that("block-tridiagonal solve equals a dense solve", {
  set.seed(42)
  for (q in c(1L, 2L, 3L, 5L)) {
    for (Tn in c(1L, 2L, 3L, 6L)) {
      off <- -1.7
      D <- mk_D(q, off, Tn)
      r <- lapply(seq_len(Tn), function(t) stats::rnorm(q))
      x <- .solve_block_tridiagonal(D, off, r)
      expect_equal(unlist(x), as.vector(solve(dense_H(D, off), unlist(r))),
                   tolerance = 1e-9,
                   info = sprintf("q = %d, T = %d", q, Tn))
    }
  }
})


test_that("the solver accepts bare scalars and decouples when off = 0", {
  off <- -0.5
  D <- list(3, 4, 5)
  r <- list(1, 2, 3)
  expect_equal(unlist(.solve_block_tridiagonal(D, off, r)),
               as.vector(solve(dense_H(D, off), unlist(r))), tolerance = 1e-10)

  set.seed(1)
  D <- mk_D(3L, 0, 4L)
  r <- lapply(1:4, function(t) stats::rnorm(3))
  expect_equal(unlist(.solve_block_tridiagonal(D, 0, r)),
               unlist(lapply(1:4, function(t) as.numeric(solve(D[[t]], r[[t]])))),
               tolerance = 1e-12)
})


test_that("the solver validates its inputs", {
  expect_error(.solve_block_tridiagonal(list(), 1, list()), "at least one")
  expect_error(.solve_block_tridiagonal(list(1, 2), 1, list(1)), "same length")
  expect_error(.solve_block_tridiagonal(list(1), c(1, 2), list(1)), "single finite")
})


test_that("the vectorised scalar solver matches the block solver", {
  set.seed(7)
  for (trial in 1:10) {
    n <- 6L
    Tn <- sample(1:7, 1)
    off <- -stats::runif(1, 0.1, 3)
    D <- matrix(stats::runif(n * Tn, 1, 5), n, Tn) + 2 * abs(off)
    r <- matrix(stats::rnorm(n * Tn), n, Tn)
    xv <- .solve_scalar_tridiagonal_vec(D, off, r)
    for (i in seq_len(n)) {
      xb <- .solve_block_tridiagonal(as.list(D[i, ]), off, as.list(r[i, ]))
      expect_equal(unlist(xb), xv[i, ], tolerance = 1e-10)
    }
  }
})


test_that("a singular additive block is reported clearly", {
  expect_error(
    .solve_scalar_tridiagonal_vec(matrix(0, 2, 3), 0, matrix(1, 2, 3)),
    "Singular additive block"
  )
})


test_that("selected inversion equals the explicit inverse", {
  set.seed(17)
  for (q in c(1L, 2L, 3L, 4L)) {
    for (Tn in c(1L, 2L, 3L, 5L, 8L)) {
      off <- -stats::runif(1, 0.5, 4)
      D <- mk_D(q, off, Tn)
      si <- .block_tridiagonal_selected_inverse(D, off)
      Hinv <- solve(dense_H(D, off))

      for (t in seq_len(Tn)) {
        idx <- ((t - 1) * q + 1):(t * q)
        expect_equal(si$diag[[t]], Hinv[idx, idx, drop = FALSE],
                     tolerance = 1e-9, ignore_attr = TRUE,
                     info = sprintf("diag q = %d, T = %d, t = %d", q, Tn, t))
      }
      if (Tn > 1L) {
        for (t in seq_len(Tn - 1L)) {
          idx <- ((t - 1) * q + 1):(t * q)
          nx <- (t * q + 1):((t + 1) * q)
          expect_equal(si$offdiag[[t]], Hinv[nx, idx, drop = FALSE],
                       tolerance = 1e-9, ignore_attr = TRUE)
        }
      }
    }
  }
})


test_that("selected inversion matches the 2x2 closed form", {
  d1 <- 3; d2 <- 5; e <- -1.5
  si <- .block_tridiagonal_selected_inverse(list(d1, d2), e)
  det <- d1 * d2 - e^2
  expect_equal(si$diag[[1]][1, 1], d2 / det, tolerance = 1e-12)
  expect_equal(si$diag[[2]][1, 1], d1 / det, tolerance = 1e-12)
  expect_equal(si$offdiag[[1]][1, 1], -e / det, tolerance = 1e-12)
})


test_that("the vectorised scalar selected inverse matches the block version", {
  set.seed(23)
  n <- 5L; Tn <- 6L; off <- -1.4
  D <- matrix(stats::runif(n * Tn, 2, 6), n, Tn)
  sv <- .scalar_tridiagonal_selected_inverse_vec(D, off)
  for (i in seq_len(n)) {
    sb <- .block_tridiagonal_selected_inverse(as.list(D[i, ]), off)
    expect_equal(sv$diag[i, ], vapply(sb$diag, function(M) M[1, 1], numeric(1)),
                 tolerance = 1e-10)
    expect_equal(sv$offdiag[i, ],
                 vapply(sb$offdiag, function(M) M[1, 1], numeric(1)),
                 tolerance = 1e-10)
  }
})


test_that("the trace combination the EM step needs is correct", {
  # tr(Var(theta_t - theta_{t-1})) built from the selected inverse must equal
  # the same quantity formed directly from the dense inverse.
  set.seed(31)
  q <- 2L; Tn <- 5L; off <- -1.3
  D <- mk_D(q, off, Tn)
  si <- .block_tridiagonal_selected_inverse(D, off)
  Hinv <- solve(dense_H(D, off))

  for (t in 2:Tn) {
    from_sel <- sum(diag(si$diag[[t]])) + sum(diag(si$diag[[t - 1]])) -
      2 * sum(diag(si$offdiag[[t - 1]]))
    A <- matrix(0, q, Tn * q)
    A[, ((t - 1) * q + 1):(t * q)] <- diag(q)
    A[, ((t - 2) * q + 1):((t - 1) * q)] <- -diag(q)
    expect_equal(from_sel, sum(diag(A %*% Hinv %*% t(A))), tolerance = 1e-9)
    expect_gt(from_sel, 0)
  }
})


test_that("Procrustes recovers a known orthogonal transform, reflections included", {
  set.seed(11)
  K <- 3L
  Rtrue <- qr.Q(qr(matrix(stats::rnorm(K * K), K, K)))
  Blist <- lapply(1:4, function(t) matrix(stats::rnorm(10 * K), 10, K))
  Alist <- lapply(Blist, function(B) B %*% t(Rtrue))

  Rhat <- .orthogonal_procrustes(Alist, Blist)
  expect_equal(Rhat, Rtrue, tolerance = 1e-9)
  expect_equal(crossprod(Rhat), diag(K), tolerance = 1e-10)

  # a reflection must be recoverable too -- U V' is invariant under any
  # orthogonal transform, not only proper rotations
  Rref <- Rtrue
  Rref[, K] <- -Rref[, K]
  Alist2 <- lapply(Blist, function(B) B %*% t(Rref))
  Rhat2 <- .orthogonal_procrustes(Alist2, Blist)
  expect_equal(Rhat2, Rref, tolerance = 1e-9)
  expect_lt(det(Rhat2), 0)
})


test_that("pooled normalisation follows the spec and preserves U V'", {
  set.seed(3)
  N <- 25L; M <- 7L; K <- 3L; Tn <- 4L      # N != M on purpose
  Ul <- lapply(1:Tn, function(t) matrix(stats::rnorm(N * K), N, K) * 4)
  Vl <- lapply(1:Tn, function(t) matrix(stats::rnorm(M * K), M, K) * 0.3)

  SU <- sum(vapply(Ul, function(U) sum(U^2), numeric(1))) / (N * K * Tn)
  SV <- sum(vapply(Vl, function(V) sum(V^2), numeric(1))) / (M * K * Tn)

  sn <- .scale_normalize_UV(Ul, Vl, anchor = "pooled")
  expect_equal(sn$S_U, SU, tolerance = 1e-12)
  expect_equal(sn$S_V, SV, tolerance = 1e-12)
  expect_equal(sn$c, (SV / SU)^(1 / 4), tolerance = 1e-12)

  # the convention the spec states
  lhs <- sum(vapply(sn$U, function(U) sum(U^2), numeric(1))) / (N * K * Tn)
  rhs <- sum(vapply(sn$V, function(V) sum(V^2), numeric(1))) / (M * K * Tn)
  expect_equal(lhs, rhs, tolerance = 1e-10)

  # the fitted interaction is untouched
  for (t in seq_len(Tn)) {
    expect_equal(sn$U[[t]] %*% t(sn$V[[t]]), Ul[[t]] %*% t(Vl[[t]]),
                 tolerance = 1e-10)
  }

  # idempotent
  expect_equal(.scale_normalize_UV(sn$U, sn$V, anchor = "pooled")$c, 1,
               tolerance = 1e-12)
})


test_that("first-period normalisation equalises the initial states and preserves U V'", {
  set.seed(3)
  N <- 25L; M <- 7L; K <- 3L; Tn <- 4L
  Ul <- lapply(1:Tn, function(t) matrix(stats::rnorm(N * K), N, K) * 4)
  Vl <- lapply(1:Tn, function(t) matrix(stats::rnorm(M * K), M, K) * 0.3)

  sn <- .scale_normalize_UV(Ul, Vl, anchor = "first")
  expect_equal(sn$S_U, sum(Ul[[1]]^2) / (N * K), tolerance = 1e-12)
  expect_equal(sn$S_V, sum(Vl[[1]]^2) / (M * K), tolerance = 1e-12)

  # the convention: the two latent spaces are the same size at t = 1, and
  # nothing is claimed about any later period
  expect_equal(sum(sn$U[[1]]^2) / (N * K), sum(sn$V[[1]]^2) / (M * K),
               tolerance = 1e-10)

  for (t in seq_len(Tn)) {
    expect_equal(sn$U[[t]] %*% t(sn$V[[t]]), Ul[[t]] %*% t(Vl[[t]]),
                 tolerance = 1e-10)
  }
  expect_equal(.scale_normalize_UV(sn$U, sn$V, anchor = "first")$c, 1,
               tolerance = 1e-12)
})


test_that("the innovation convention equalises the increments, not the magnitudes", {
  # `anchor = "innovation"` makes the two sides agree on mean squared
  # increment. The magnitudes are then NOT equal, and that is the point:
  # imposing both would use one degree of freedom twice.
  #
  # It is not the default. Equalising the increments requires rescaling the
  # trajectories, and the factor does not return to one between outer passes --
  # the correction compounds, and after forty passes the trajectories carry a
  # factor of 1e7 and the solve goes singular. The same equality is imposed
  # instead by `tie_latent` in `.eb_update_variances()`, which constrains the
  # reported variance and never touches the trajectories.
  set.seed(41)
  N <- 30L; M <- 18L; K <- 2L; Tn <- 6L
  Ul <- Vl <- vector("list", Tn)
  Ul[[1]] <- matrix(stats::rnorm(N * K), N, K) * 3
  Vl[[1]] <- matrix(stats::rnorm(M * K), M, K) * 0.4
  for (t in 2:Tn) {
    Ul[[t]] <- Ul[[t-1]] + matrix(stats::rnorm(N * K, sd = 0.30), N, K)
    Vl[[t]] <- Vl[[t-1]] + matrix(stats::rnorm(M * K, sd = 0.10), M, K)
  }

  sn <- .scale_normalize_UV(Ul, Vl, anchor = "innovation")
  inc <- function(L, n) sum(vapply(2:Tn, function(t) sum((L[[t]] - L[[t-1]])^2),
                                   numeric(1))) / (n * K * (Tn - 1L))
  expect_equal(inc(sn$U, N), inc(sn$V, M), tolerance = 1e-10)

  # the magnitudes are left unequal, and the fit is untouched
  expect_gt(abs(sum(sn$U[[1]]^2) / (N * K) - sum(sn$V[[1]]^2) / (M * K)), 1e-6)
  for (t in seq_len(Tn))
    expect_equal(sn$U[[t]] %*% t(sn$V[[t]]), Ul[[t]] %*% t(Vl[[t]]),
                 tolerance = 1e-10)

  expect_equal(.scale_normalize_UV(sn$U, sn$V, anchor = "innovation")$c, 1,
               tolerance = 1e-12)
})


test_that("a single period has no increments and falls back to the initial states", {
  set.seed(9)
  U1 <- list(matrix(stats::rnorm(20), 10, 2))
  V1 <- list(matrix(stats::rnorm(12), 6, 2))
  expect_equal(.scale_normalize_UV(U1, V1)$c,
               .scale_normalize_UV(U1, V1, anchor = "first")$c,
               tolerance = 1e-12)
})


test_that("the innovation ratio survives a rescaling of the whole trajectory", {
  # Under U -> aU, V -> V/a the factor moves to c/a and the two changes cancel,
  # so the ratio of the normalised innovation variances is untouched. This
  # holds for BOTH anchors, which is worth pinning down: it means a pure scale
  # drift is not what sends the ratio off, and that choosing between the
  # anchors is a question about bias rather than about invariance.
  set.seed(11)
  N <- 30L; M <- 18L; K <- 2L; Tn <- 6L
  Ul <- Vl <- vector("list", Tn)
  Ul[[1]] <- matrix(stats::rnorm(N * K), N, K)
  Vl[[1]] <- matrix(stats::rnorm(M * K), M, K)
  for (t in 2:Tn) {
    Ul[[t]] <- Ul[[t-1]] + matrix(stats::rnorm(N * K, sd = 0.30), N, K)
    Vl[[t]] <- Vl[[t-1]] + matrix(stats::rnorm(M * K, sd = 0.10), M, K)
  }

  ratio_of <- function(U, V, anchor) {
    sn <- .scale_normalize_UV(U, V, anchor = anchor)
    tu <- sum(vapply(2:Tn, function(t) sum((sn$U[[t]] - sn$U[[t-1]])^2),
                     numeric(1))) / (N * K * (Tn - 1))
    tv <- sum(vapply(2:Tn, function(t) sum((sn$V[[t]] - sn$V[[t-1]])^2),
                     numeric(1))) / (M * K * (Tn - 1))
    tu / tv
  }

  a <- 7.3
  Ua <- lapply(Ul, function(U) a * U)
  Va <- lapply(Vl, function(V) V / a)
  for (anc in c("first", "pooled")) {
    expect_equal(ratio_of(Ua, Va, anc), ratio_of(Ul, Vl, anc),
                 tolerance = 1e-9)
  }
})


test_that("pooling the anchor pulls a real asymmetry toward symmetry", {
  # The two trajectories start at the same per-coordinate scale and U then
  # moves nine times as fast. Anchored at t = 1 the convention leaves that
  # alone; pooled, it does not — U spreads further as time passes, the pooled
  # magnitude picks the spread up, and the factor absorbs part of the very
  # asymmetry the innovation variances exist to report.
  set.seed(202)
  N <- 40L; M <- 40L; K <- 2L; Tn <- 8L
  Ul <- Vl <- vector("list", Tn)
  Ul[[1]] <- matrix(stats::rnorm(N * K), N, K)
  Vl[[1]] <- matrix(stats::rnorm(M * K), M, K)
  for (t in 2:Tn) {
    Ul[[t]] <- Ul[[t-1]] + matrix(stats::rnorm(N * K, sd = 0.30), N, K)
    Vl[[t]] <- Vl[[t-1]] + matrix(stats::rnorm(M * K, sd = 0.10), M, K)
  }

  ratio_of <- function(anchor) {
    sn <- .scale_normalize_UV(Ul, Vl, anchor = anchor)
    tu <- sum(vapply(2:Tn, function(t) sum((sn$U[[t]] - sn$U[[t-1]])^2),
                     numeric(1))) / (N * K * (Tn - 1))
    tv <- sum(vapply(2:Tn, function(t) sum((sn$V[[t]] - sn$V[[t-1]])^2),
                     numeric(1))) / (M * K * (Tn - 1))
    tu / tv
  }

  r_first <- ratio_of("first")
  r_pooled <- ratio_of("pooled")

  expect_gt(r_first, r_pooled)          # pooling shrinks it
  expect_gt(r_first, 4)                 # and the true ratio here is about 9
  expect_lt(r_pooled, r_first * 0.8)    # by a margin that is not rounding
})


test_that("scale normalisation warns instead of producing NaN when degenerate", {
  Uz <- lapply(1:2, function(t) matrix(0, 5, 2))
  Vz <- lapply(1:2, function(t) matrix(stats::rnorm(6), 3, 2))
  expect_warning(res <- .scale_normalize_UV(Uz, Vz), "degenerate")
  expect_equal(res$c, 1)
  expect_true(all(is.finite(unlist(res$U))))
  expect_true(all(is.finite(unlist(res$V))))
})
