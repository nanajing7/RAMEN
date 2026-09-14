# Shared fixtures for the dynamic AME tests.

#' Simulate a panel from the model itself
#'
#' Returns the true parameters alongside the data, so tests can check recovery
#' rather than only internal consistency.
sim_ame_panel <- function(N = 9, M = 7, K = 2, Tn = 4,
                          Pr = 1, Pc = 1, Pd = 1,
                          sd_eps = 0.5, na_frac = 0.1, seed = 1) {
  set.seed(seed)
  P <- Pr + Pc + Pd

  s <- list(a = matrix(0, N, Tn), b = matrix(0, M, Tn),
            beta = matrix(0, P, Tn),
            U = vector("list", Tn), V = vector("list", Tn))
  s$a[, 1] <- stats::rnorm(N)
  s$b[, 1] <- stats::rnorm(M, sd = 0.8)
  if (P > 0) s$beta[, 1] <- stats::rnorm(P, sd = 0.6)
  s$U[[1]] <- matrix(stats::rnorm(N * K, sd = 0.9), N, K)
  s$V[[1]] <- matrix(stats::rnorm(M * K, sd = 0.7), M, K)
  for (t in seq_len(Tn)[-1]) {
    s$a[, t] <- s$a[, t - 1] + stats::rnorm(N, sd = 0.35)
    s$b[, t] <- s$b[, t - 1] + stats::rnorm(M, sd = 0.30)
    if (P > 0) s$beta[, t] <- s$beta[, t - 1] + stats::rnorm(P, sd = 0.15)
    s$U[[t]] <- s$U[[t - 1]] + matrix(stats::rnorm(N * K, sd = 0.30), N, K)
    s$V[[t]] <- s$V[[t - 1]] + matrix(stats::rnorm(M * K, sd = 0.28), M, K)
  }

  X_row  <- if (Pr > 0) lapply(seq_len(Tn), function(t) matrix(stats::rnorm(N * Pr), N, Pr)) else NULL
  X_col  <- if (Pc > 0) lapply(seq_len(Tn), function(t) matrix(stats::rnorm(M * Pc), M, Pc)) else NULL
  X_dyad <- if (Pd > 0) lapply(seq_len(Tn), function(t) array(stats::rnorm(N * M * Pd), c(N, M, Pd))) else NULL

  Y <- lapply(seq_len(Tn), function(t) {
    y <- .fitted_period(
      s$a[, t], s$b[, t], s$U[[t]], s$V[[t]],
      if (P > 0) s$beta[, t] else numeric(0),
      if (is.null(X_row)) NULL else X_row[[t]],
      if (is.null(X_col)) NULL else X_col[[t]],
      if (is.null(X_dyad)) NULL else X_dyad[[t]]
    ) + matrix(stats::rnorm(N * M, sd = sd_eps), N, M)
    if (na_frac > 0) y[sample(N * M, round(na_frac * N * M))] <- NA
    dimnames(y) <- list(paste0("r", seq_len(N)), paste0("c", seq_len(M)))
    y
  })

  list(sim = s, Y = Y, X_row = X_row, X_col = X_col, X_dyad = X_dyad,
       N = N, M = M, K = K, Tn = Tn, P = P,
       dims = list(Pr = Pr, Pc = Pc, Pd = Pd))
}

#' A random parameter list matching a simulated panel's shape
random_params <- function(d, seed = 99) {
  set.seed(seed)
  list(
    a = matrix(stats::rnorm(d$N * d$Tn), d$N, d$Tn),
    b = matrix(stats::rnorm(d$M * d$Tn), d$M, d$Tn),
    beta = matrix(stats::rnorm(d$P * d$Tn), d$P, d$Tn),
    U = lapply(seq_len(d$Tn), function(t) matrix(stats::rnorm(d$N * d$K), d$N, d$K)),
    V = lapply(seq_len(d$Tn), function(t) matrix(stats::rnorm(d$M * d$K), d$M, d$K))
  )
}

#' Neutral penalties, the spec's starting values
unit_penalty <- function() {
  stats::setNames(rep(1, 5), c("beta", "a", "b", "U", "V"))
}

#' Long-format edge panel from a simulated one, for the public interface
as_edge_panel <- function(d, years = NULL) {
  years <- years %||% seq_len(d$Tn)
  do.call(rbind, lapply(seq_len(d$Tn), function(t) {
    g <- expand.grid(i = seq_len(d$N), j = seq_len(d$M))
    keep <- !is.na(d$Y[[t]][cbind(g$i, g$j)])
    data.frame(
      year = years[t],
      node_row = sprintf("%03d", g$i[keep]),
      node_col = sprintf("%03d", g$j[keep]),
      value = d$Y[[t]][cbind(g$i, g$j)][keep],
      stringsAsFactors = FALSE
    )
  }))
}
