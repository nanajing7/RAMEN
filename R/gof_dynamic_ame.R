# Goodness of fit for fit_dynamic_ame().
#
# Scope: one fit in, one description out. Nothing here refits a model.
#
# That line is deliberate and it decides what belongs. "Does the multiplicative
# term earn its place?" and "is the random walk the right transition law?" are
# model-comparison questions: answering either needs a second model estimated on
# the same data, so neither is here. What is here is whether the fit that was
# produced describes the panel it was produced from.
#
# Five groups:
#
#   A  conditional replication   does the fitted structure, with fresh
#                                observation errors, reproduce marginal
#                                heterogeneity, adjacent-period similarity, and
#                                the size of period-to-period change?
#   B  pair structure            is there residual cross-dyad structure that an
#                                additive fit plus independent noise does not
#                                produce?
#   C  fit accounting            what does each block explain, in sample?
#   D  estimation diagnostics    did the estimator converge without hitting a
#                                limit?
#   F  residual screening        is there row, column, temporal, or
#                                distributional structure left over?
#
# ── What A is, and what it is not ────────────────────────────────────────────
#
# A is NOT a posterior predictive check, and its intervals are not posterior
# predictive intervals. Each replicate holds `U` and `V` at their estimates,
# applies the limited dispersion adjustment `.dispersion_scales()` supplies for
# `beta`, `a` and `b`, and adds fresh Gaussian observation errors. The interval
# is therefore narrow in an amount that cannot be quantified, and its centre is
# built from the same data the observed statistic comes from. Both defects push
# in directions that are not known in advance, so the intervals are read as
# descriptive references, never as tests:
#
#   An observed value inside the interval is compatible with the fitted
#   conditional structure; a value outside indicates a discrepancy worth
#   inspecting.
#
# Consequences, enforced throughout: `inside` is never a red flag,
# `tail_fraction` is a fraction of replicates and not a p-value, `sd.cell` is
# flagged `by_construction` because the replicate variance is
# `var(fitted) + sigma_eps^2` while `sigma_eps^2` is itself the residual mean
# square (`.eb_update_variances()`), so the leading terms cancel, and
# `sd.delta.cell` is flagged because on a common cell set it is an exact
# function of `sd.cell` at both periods and `cor.lag1` between them. Neither
# flagged statistic is independent evidence; both are reported because scale and
# movement in outcome units are what a reader can interpret.
#
# `cor.lag1` is adjacent-panel similarity, not a test of the random-walk prior:
# stable covariates, persistent row or column heterogeneity, and a latent
# structure held fixed across replicates all raise it.


# ── Statistic helpers ────────────────────────────────────────────────────────

#' Reading order for the statistics
#'
#' Overall scale, then the two margins, then time. Alphabetical order splits
#' the margins apart and puts a temporal statistic first, which is not how any
#' of it is read.
#'
#' @keywords internal
#' @noRd
.GOF_STAT_ORDER <- c("sd.cell", "sd.rowmean", "sd.colmean", "cor.lag1",
                     "sd.delta.cell", "sd.rowcor", "sd.colcor")

#' Correlation matrix, skipping the pairwise machinery when nothing is missing
#'
#' `use = "pairwise.complete.obs"` visits every pair of variables separately,
#' because different pairs can have different sets of jointly observed rows. It
#' does so even when the matrix is complete and every pair would use every row,
#' which costs about two and a half times a plain `cor()` for no gain: with
#' nothing missing the two agree to floating-point noise.
#'
#' @param X Numeric matrix, variables in columns.
#' @return The variable-by-variable correlation matrix.
#' @keywords internal
#' @noRd
.cor_safe <- function(X) {
  if (anyNA(X)) stats::cor(X, use = "pairwise.complete.obs") else stats::cor(X)
}

#' Strictly upper triangle of a square matrix, as a vector
#'
#' @param m A square matrix.
#' @return The strictly upper triangular entries.
#' @keywords internal
#' @noRd
.gof_upper <- function(m) m[upper.tri(m)]

#' Standard deviation of the profile correlations along one margin
#'
#' Correlates the rows of `Y` with each other (`by = "row"`) or the columns with
#' each other, and returns how widely those correlations spread.
#'
#' A purely additive panel cannot spread them. If `Y_ij = a_i + b_j` then row
#' `i` and row `i'` differ by the constant `a_i - a_i'`, and a correlation is
#' unchanged by adding a constant to either variable, so every pair correlates
#' at one and the spread is zero. Noise pulls the correlations off one, but
#' symmetrically; structure that varies the *shape* of a row's profile, not just
#' its level, is what makes the spread large.
#'
#' @param Y An `N x M` outcome matrix, `NA` where unobserved.
#' @param by `"row"` correlates rows, `"col"` correlates columns.
#' @param min_overlap Minimum jointly observed cells a pair needs before its
#'   correlation is used. Pairs below it are dropped, not imputed: a correlation
#'   from three shared cells is noise wearing the units of a statistic.
#' @return A named numeric vector with `value` and `n_used` (pairs retained).
#' @keywords internal
#' @noRd
.stat_profile_cor <- function(Y, by = c("row", "col"), min_overlap = 5L) {
  by <- match.arg(by)
  X <- if (by == "row") t(Y) else Y      # the margin of interest in columns
  if (is.null(ncol(X)) || ncol(X) < 2L)
    return(c(value = NA_real_, n_used = 0))

  C <- suppressWarnings(.cor_safe(X))

  if (anyNA(X)) {
    # crossprod of the observation indicator counts, for every pair of
    # variables, the rows where both were seen.
    ov <- crossprod((!is.na(X)) + 0)
    C[ov < min_overlap] <- NA_real_
  }

  u <- .gof_upper(C)
  u <- u[is.finite(u)]
  if (length(u) < 2L) return(c(value = NA_real_, n_used = length(u)))
  c(value = stats::sd(u), n_used = length(u))
}

#' Margin means, dropping wholly unobserved rows or columns
#'
#' A row with no observed cell has no mean. `rowMeans(na.rm = TRUE)` returns
#' `NaN` there, which would propagate into the spread; it is removed and counted
#' instead.
#'
#' @param Y An outcome matrix.
#' @param by `"row"` or `"col"`.
#' @return A named numeric vector with `value` (the spread) and `n_used`.
#' @keywords internal
#' @noRd
.stat_margin_sd <- function(Y, by = c("row", "col")) {
  by <- match.arg(by)
  m <- if (by == "row") rowMeans(Y, na.rm = TRUE) else colMeans(Y, na.rm = TRUE)
  m <- m[is.finite(m)]
  if (length(m) < 2L) return(c(value = NA_real_, n_used = length(m)))
  c(value = stats::sd(m), n_used = length(m))
}


# ── Observed statistics ──────────────────────────────────────────────────────

#' Goodness-of-fit statistics for a panel of networks
#'
#' Computes the seven statistics [gof_dynamic_ame()] compares, directly from a
#' panel. No fitted model is involved, so the same numbers can be taken from a
#' simulated panel, a competing model's replicates, or the data itself.
#'
#' Five are per period:
#'
#' \describe{
#'   \item{`sd.cell`}{Spread of the cells. A scale reference.}
#'   \item{`sd.rowmean`}{Spread of the row means: how unequal the row nodes are.}
#'   \item{`sd.colmean`}{Spread of the column means: how unequal the column
#'     nodes are.}
#'   \item{`sd.rowcor`}{Spread of the correlations between row profiles. Large
#'     when rows differ in the *shape* of their profile across columns, which no
#'     additive structure produces; see `.stat_profile_cor()`.}
#'   \item{`sd.colcor`}{The same along the other margin.}
#' }
#'
#' Two describe a transition, and are labelled by the period they arrive at, so
#' they plot on the same axis as the rest:
#'
#' \describe{
#'   \item{`cor.lag1`}{Correlation between consecutive panels: how similar
#'     adjacent periods are.}
#'   \item{`sd.delta.cell`}{Spread of the cell-by-cell changes: how much the
#'     panel moves, in the units of the outcome.}
#' }
#'
#' Both transition statistics use only cells observed in both periods. On a
#' common cell set the two are related through `sd.cell`, by
#' `var(Y_t - Y_{t-1}) = var(Y_t) + var(Y_{t-1}) - 2 cov(Y_t, Y_{t-1})`, so
#' `sd.delta.cell` carries no information the other three lack. It is reported
#' because a movement in outcome units is legible where a correlation is not.
#'
#' @param Y_list List of `T` outcome matrices, `NA` where unobserved.
#' @param min_overlap Minimum jointly observed cells for a pair of rows or
#'   columns to contribute to `sd.rowcor` or `sd.colcor`.
#' @param which Statistics to compute. Defaults to all seven; the two
#'   correlation-spread statistics are much the most expensive, so the
#'   replication group asks only for the ones it uses.
#' @param periods Optional labels, one per element of `Y_list`.
#'
#' @return A data frame with `statistic`, `period`, `value` and `n_used`.
#'
#' @examples
#' \dontrun{
#' gof_stats_ame(fit$Y_list)
#' }
#' @export
gof_stats_ame <- function(Y_list,
                          min_overlap = 5L,
                          which = c("sd.cell", "sd.rowmean", "sd.colmean",
                                    "cor.lag1", "sd.delta.cell",
                                    "sd.rowcor", "sd.colcor"),
                          periods = NULL) {
  if (!is.list(Y_list) || !length(Y_list))
    stop("Y_list must be a non-empty list of outcome matrices.", call. = FALSE)
  which <- match.arg(which, several.ok = TRUE)

  Tn <- length(Y_list)
  if (is.null(periods)) periods <- names(Y_list) %||% as.character(seq_len(Tn))
  periods <- as.character(periods)
  if (length(periods) != Tn)
    stop("periods must have one label per period.", call. = FALSE)

  rows <- list()
  add <- function(stat, per, v) {
    rows[[length(rows) + 1L]] <<- data.frame(
      statistic = stat, period = per,
      value = unname(v[["value"]]), n_used = unname(v[["n_used"]]),
      stringsAsFactors = FALSE)
  }

  for (t in seq_len(Tn)) {
    Y <- Y_list[[t]]
    if ("sd.cell" %in% which) {
      ok <- sum(!is.na(Y))
      add("sd.cell", periods[t],
          c(value = if (ok < 2L) NA_real_ else stats::sd(Y, na.rm = TRUE),
            n_used = ok))
    }
    if ("sd.rowmean" %in% which) add("sd.rowmean", periods[t],
                                     .stat_margin_sd(Y, "row"))
    if ("sd.colmean" %in% which) add("sd.colmean", periods[t],
                                     .stat_margin_sd(Y, "col"))
    if ("sd.rowcor" %in% which)  add("sd.rowcor", periods[t],
                                     .stat_profile_cor(Y, "row", min_overlap))
    if ("sd.colcor" %in% which)  add("sd.colcor", periods[t],
                                     .stat_profile_cor(Y, "col", min_overlap))
  }

  if (Tn > 1L && any(c("cor.lag1", "sd.delta.cell") %in% which)) {
    for (t in 2:Tn) {
      prev <- as.vector(Y_list[[t - 1L]])
      cur <- as.vector(Y_list[[t]])
      keep <- !is.na(prev) & !is.na(cur)
      n <- sum(keep)
      if ("cor.lag1" %in% which) {
        v <- if (n < 3L) NA_real_ else
          suppressWarnings(stats::cor(prev[keep], cur[keep]))
        add("cor.lag1", periods[t],
            c(value = if (is.finite(v)) v else NA_real_, n_used = n))
      }
      if ("sd.delta.cell" %in% which) {
        add("sd.delta.cell", periods[t],
            c(value = if (n < 2L) NA_real_ else stats::sd(cur[keep] - prev[keep]),
              n_used = n))
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ── A: conditional replication ───────────────────────────────────────────────

#' Statistics of one conditional replicate
#'
#' @param fit A `dynamic_ame` object.
#' @param disp Output of `.dispersion_scales()`.
#' @param which Statistics to compute.
#' @param min_overlap Passed through.
#' @return The data frame `gof_stats_ame()` returns.
#' @keywords internal
#' @noRd
.gof_one_replicate <- function(fit, disp, which, min_overlap) {
  tr <- .draw_dispersed_trajectories(fit, disp)
  mu <- .fitted_from_params(fit, tr$beta, tr$a, tr$b)
  Ys <- .simulate_conditional_panel(mu, fit$Y_list, fit$sigma_eps2)
  gof_stats_ame(Ys, min_overlap = min_overlap, which = which,
                periods = as.character(fit$years))
}

#' Reference distribution for a set of statistics
#'
#' Runs `sim_one()` `nsim` times, stacks the results, and summarises each
#' statistic-by-period cell. The tail fraction is two-sided about the reference
#' median, matching the absence of any directional claim in group A; group B
#' overrides it, because its claim has a direction.
#'
#' @param sim_one A function of no arguments returning one replicate's statistic
#'   table.
#' @param observed The observed statistic table.
#' @param nsim,conf_level,seed,n_cores Simulation controls.
#' @param two_sided Tail fraction convention.
#' @param keep_draws Return every replicate's statistics.
#' @return A list with `summary` (one row per statistic and period) and,
#'   optionally, `draws`.
#' @keywords internal
#' @noRd
.gof_reference <- function(sim_one, observed, nsim, conf_level, seed,
                           n_cores = 1L, two_sided = TRUE, keep_draws = FALSE) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, nsim)

  one <- function(b) { set.seed(rep_seeds[b]); sim_one() }
  reps <- if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE) &&
              .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nsim), one, mc.cores = n_cores)
  } else {
    lapply(seq_len(nsim), one)
  }

  key <- paste(observed$statistic, observed$period, sep = "\r")
  # one column per replicate, rows aligned to `observed`
  draws <- vapply(reps, function(d)
    d$value[match(key, paste(d$statistic, d$period, sep = "\r"))],
    numeric(nrow(observed)))
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = nrow(observed))

  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  q <- function(p) apply(draws, 1, stats::quantile, probs = p, na.rm = TRUE)
  med <- apply(draws, 1, stats::median, na.rm = TRUE)

  tail_fraction <- vapply(seq_len(nrow(observed)), function(i) {
    d <- draws[i, ]
    d <- d[is.finite(d)]
    o <- observed$value[i]
    if (!length(d) || !is.finite(o)) return(NA_real_)
    if (two_sided) mean(abs(d - med[i]) >= abs(o - med[i])) else mean(d >= o)
  }, numeric(1))

  lo <- q(probs[1]); hi <- q(probs[2])
  out <- list(
    summary = data.frame(
      statistic = observed$statistic,
      period = observed$period,
      observed = observed$value,
      ref_median = as.numeric(med),
      ref_lo = as.numeric(lo),
      ref_hi = as.numeric(hi),
      ref_sd = apply(draws, 1, stats::sd, na.rm = TRUE),
      inside = observed$value >= lo & observed$value <= hi,
      tail_fraction = tail_fraction,
      n_used = observed$n_used,
      stringsAsFactors = FALSE),
    n_ok = sum(vapply(reps, function(d) !all(is.na(d$value)), logical(1))))
  if (keep_draws) out$draws <- draws
  out
}


# ── B: matched-moment additive reference ─────────────────────────────────────

#' Additive reference panels
#'
#' Rebuilds the panel from the fitted covariate and additive terms and replaces
#' everything else with noise carrying the same first two moments:
#'
#' \deqn{Y^{null}_t = cov_t + add_t + \bar r_t + \epsilon_t,
#'       \quad \epsilon_t \sim N(0, s_t^2)}
#'
#' with \eqn{\bar r_t} and \eqn{s_t^2} the mean and variance of
#' \eqn{Y_t - cov_t - add_t} over the observed cells. Matching both moments is
#' what makes the reference differ from the data in one respect only: the
#' pair-level variation is present in the same quantity, but unstructured.
#'
#' The mean matters less than it looks for the two statistics used here, since a
#' correlation is unchanged by a constant added to every cell; the centring of
#' \eqn{s_t^2} is the part that bites, an uncentred second moment being larger
#' and so inflating the reference.
#'
#' `beta`, `a` and `b` are not perturbed. The comparison is conditional on the
#' estimated covariate and additive structure, which is also why the reference
#' is narrower than one accounting for their estimation error, and why `excess`
#' is reported as a magnitude rather than as a test statistic.
#'
#' @param fit A `dynamic_ame` object.
#' @return A list with `base` (covariate plus additive plus the residual mean,
#'   per period) and `sd` (the matched noise scale, per period).
#' @keywords internal
#' @noRd
.gof_null_parts <- function(fit) {
  dec <- decompose_fit(fit)
  Tn <- length(dec)

  base <- vector("list", Tn)
  sdev <- numeric(Tn)
  for (t in seq_len(Tn)) {
    b <- dec[[t]]$covariate + dec[[t]]$additive
    r <- fit$Y_list[[t]] - b
    ok <- is.finite(r)
    if (sum(ok) < 2L) {
      base[[t]] <- b
      sdev[t] <- 0
      next
    }
    rbar <- mean(r[ok])
    base[[t]] <- b + rbar
    sdev[t] <- sqrt(mean((r[ok] - rbar)^2))
  }
  names(base) <- names(dec)
  list(base = base, sd = sdev)
}

#' One draw from the additive reference
#'
#' @param parts Output of `.gof_null_parts()`.
#' @param Y_list Original outcomes, for the missing pattern.
#' @param which,min_overlap,periods Passed to `gof_stats_ame()`.
#' @return The data frame `gof_stats_ame()` returns.
#' @keywords internal
#' @noRd
.gof_one_null <- function(parts, Y_list, which, min_overlap, periods) {
  Yn <- lapply(seq_along(parts$base), function(t) {
    B <- parts$base[[t]]
    Y <- B + matrix(stats::rnorm(length(B), sd = parts$sd[t]),
                    nrow(B), ncol(B))
    Y[is.na(Y_list[[t]])] <- NA_real_       # `B` carries its own NA already
    dimnames(Y) <- dimnames(Y_list[[t]])
    Y
  })
  names(Yn) <- names(parts$base)
  gof_stats_ame(Yn, min_overlap = min_overlap, which = which, periods = periods)
}


# ── C: in-sample accounting ──────────────────────────────────────────────────

#' Covariance accounting and nested reconstruction R-squared
#'
#' `Y` is the sum of the four blocks, so covariance shares
#' `cov(block, Y) / var(Y)` add to one exactly and split the between-block
#' covariance between the blocks that share it. Variance shares do neither.
#'
#' The nested R-squared values are reconstructions from a single fit, not three
#' fits: `beta`, `a` and `b` were estimated with the latent block in place, so
#' dropping it does not produce the additive model's own estimates. The numbers
#' describe how the fitted decomposition apportions in-sample fit. They are not
#' evidence that the latent block is worth estimating, which is a
#' model-comparison question and needs a second fit.
#'
#' @param fit A `dynamic_ame` object.
#' @return A list with `decomposition` and `increment` data frames.
#' @keywords internal
#' @noRd
.gof_accounting <- function(fit) {
  dec <- decompose_fit(fit)
  Tn <- length(dec)
  per <- names(dec)

  rows <- list()
  pool <- list()

  for (t in seq_len(Tn)) {
    Y <- fit$Y_list[[t]]
    d <- dec[[t]]
    resid <- Y - d$fitted
    ok <- is.finite(Y) & is.finite(d$covariate) & is.finite(d$additive) &
      is.finite(d$latent) & is.finite(resid)
    n <- sum(ok)
    if (n < 3L) next

    y <- Y[ok]
    blocks <- list(covariate = d$covariate[ok], additive = d$additive[ok],
                   latent = d$latent[ok], residual = resid[ok])
    vy <- stats::var(y)
    for (nm in names(blocks)) {
      rows[[length(rows) + 1L]] <- data.frame(
        period = per[t], block = nm,
        cov_share = if (vy > 0) stats::cov(blocks[[nm]], y) / vy else NA_real_,
        n_used = n, stringsAsFactors = FALSE)
    }
    pool[[length(pool) + 1L]] <- data.frame(
      y = y, cov = blocks$covariate, add = blocks$additive,
      lat = blocks$latent, stringsAsFactors = FALSE)
  }

  decomposition <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(decomposition)) rownames(decomposition) <- NULL

  increment <- NULL
  if (length(pool)) {
    p <- do.call(rbind, pool)
    sst <- sum((p$y - mean(p$y))^2)
    r2 <- function(f) if (sst > 0) 1 - sum((p$y - f)^2) / sst else NA_real_
    # No rung carries an intercept -- the level of the outcome lives in the
    # additive block -- so a rung below it can score worse than the mean and
    # come out negative. The differences are the quantity to read.
    v <- c(`covariate + additive` = r2(p$cov + p$add),
           `covariate + additive + latent` = r2(p$cov + p$add + p$lat))
    # With no covariates the first rung would be "predict zero", which is not a
    # model anyone fitted; the ladder then starts at the additive block.
    if (nrow(fit$beta) > 0L) {
      v <- c(covariate = r2(p$cov), v)
    } else {
      names(v)[1] <- "additive"
      names(v)[2] <- "additive + latent"
    }
    increment <- data.frame(model = names(v), r2 = as.numeric(v),
                            delta_r2 = c(NA_real_, diff(as.numeric(v))),
                            n_used = nrow(p), stringsAsFactors = FALSE)
    rownames(increment) <- NULL
  }

  list(decomposition = decomposition, increment = increment)
}


# ── D: estimation diagnostics ────────────────────────────────────────────────

#' Did the estimator finish cleanly
#'
#' `converged` is the answer to "did the outer loop meet its criteria". It is
#' the only iteration-related flag, deliberately.
#'
#' `outer_iterations` cannot be compared against `outer_max_iter` to decide
#' whether the budget ran out. Under `accelerate = "squarem"` the run has two
#' phases: the extrapolation gets `floor(outer_max_iter / 3)` steps at roughly
#' three passes each, and then the plain certification is given the whole of
#' `outer_max_iter` again, for the reason set out in `.ame_outer_eb()` — sharing
#' one budget made the accelerated mode converge worse than no acceleration.
#' `outer_iterations` counts every pass across both phases, so it runs to about
#' `4/3` of the cap on a healthy fit and exceeding the cap means nothing.
#' `accelerate` is reported alongside so the two numbers can be read together.
#'
#' `multistart_spread` is the range of the objective across the initial
#' candidates. Those are run through the inner block descent at *neutral*
#' penalties, before any empirical-Bayes update (see `fit_dynamic_ame()`), so it
#' measures how much the answer depends on where the search began. It does not
#' say the latent dimension is too large, and it does not count basins of the
#' converged solution. There is no calibrated threshold for it, so it is
#' reported and never flagged; a large value is a reason to refit under a
#' different seed and compare the fitted values, the objective, and the
#' coefficient trajectories.
#'
#' @param fit A `dynamic_ame` object.
#' @return A one-row data frame.
#' @keywords internal
#' @noRd
.gof_estimator <- function(fit) {
  cv <- fit$convergence
  st <- fit$settings
  mq <- cv$multistart_Q
  mq <- mq[is.finite(mq)]
  capped <- cv$penalties_capped

  data.frame(
    converged = isTRUE(cv$converged),
    outer_iterations = cv$outer_iterations %||% NA_integer_,
    outer_max_iter = st$outer_max_iter %||% NA_integer_,
    accelerate = st$accelerate %||% NA_character_,
    penalties_capped = if (length(capped))
      paste(sort(capped), collapse = ", ") else "",
    sigma_eps_floored = isTRUE(fit$sigma_eps2 <= (st$eps_var %||% 0)),
    multistart_spread = if (length(mq) > 1L) diff(range(mq)) else NA_real_,
    stringsAsFactors = FALSE
  )
}


# ── F: residual screening ────────────────────────────────────────────────────

#' Residual screening summaries
#'
#' Group A asks whether the model can generate as much row heterogeneity as the
#' data show. This asks the complementary question: whether any is left over
#' once the four blocks are removed. These are descriptions, not tests.
#'
#' @param fit A `dynamic_ame` object.
#' @param n_qq Points retained for the normal quantile plot.
#' @param n_bins Fitted-value bins for the variance-versus-level summary.
#' @return A list with `by_period`, `qq` and `scale` data frames.
#' @keywords internal
#' @noRd
.gof_residual <- function(fit, n_qq = 200L, n_bins = 10L) {
  res <- stats::residuals(fit)
  fv <- stats::fitted(fit)
  per <- names(res)
  Tn <- length(res)

  moments <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 4L) return(c(skew = NA_real_, kurt = NA_real_))
    z <- (x - mean(x)) / stats::sd(x)
    c(skew = mean(z^3), kurt = mean(z^4) - 3)
  }

  lag1 <- rep(NA_real_, Tn)
  if (Tn > 1L) {
    for (t in 2:Tn) {
      a <- as.vector(res[[t - 1L]]); b <- as.vector(res[[t]])
      k <- is.finite(a) & is.finite(b)
      if (sum(k) >= 3L) {
        v <- suppressWarnings(stats::cor(a[k], b[k]))
        lag1[t] <- if (is.finite(v)) v else NA_real_
      }
    }
  }

  by_period <- do.call(rbind, lapply(seq_len(Tn), function(t) {
    R <- res[[t]]
    m <- moments(as.vector(R))
    rm_ <- rowMeans(R, na.rm = TRUE); cm_ <- colMeans(R, na.rm = TRUE)
    data.frame(
      period = per[t],
      sd_rowmean_resid = stats::sd(rm_[is.finite(rm_)]),
      sd_colmean_resid = stats::sd(cm_[is.finite(cm_)]),
      cor_lag1_resid = lag1[t],
      sd_resid = stats::sd(R, na.rm = TRUE),
      skew = unname(m["skew"]),
      excess_kurtosis = unname(m["kurt"]),
      n_used = sum(is.finite(R)),
      stringsAsFactors = FALSE)
  }))
  rownames(by_period) <- NULL

  all_r <- unlist(lapply(res, as.vector), use.names = FALSE)
  all_f <- unlist(lapply(fv, as.vector), use.names = FALSE)
  keep <- is.finite(all_r) & is.finite(all_f)
  all_r <- all_r[keep]; all_f <- all_f[keep]

  qq <- NULL
  scale_df <- NULL
  if (length(all_r) >= 20L) {
    # store a thinned quantile pairing rather than every residual: the plot
    # needs the shape, and the object should not grow with the panel
    p <- stats::ppoints(min(n_qq, length(all_r)))
    z <- (all_r - mean(all_r)) / stats::sd(all_r)
    qq <- data.frame(theoretical = stats::qnorm(p),
                     sample = stats::quantile(z, probs = p, names = FALSE),
                     stringsAsFactors = FALSE)

    br <- stats::quantile(all_f, probs = seq(0, 1, length.out = n_bins + 1L),
                          names = FALSE)
    br <- unique(br)
    if (length(br) > 2L) {
      g <- cut(all_f, breaks = br, include.lowest = TRUE, labels = FALSE)
      scale_df <- do.call(rbind, lapply(sort(unique(g[!is.na(g)])), function(k) {
        i <- which(g == k)
        data.frame(bin = k, fitted_mid = stats::median(all_f[i]),
                   resid_sd = stats::sd(all_r[i]), n_used = length(i),
                   stringsAsFactors = FALSE)
      }))
      rownames(scale_df) <- NULL
    }
  }

  list(by_period = by_period, qq = qq, scale = scale_df)
}


# ── Entry point ──────────────────────────────────────────────────────────────

#' Goodness of fit for a dynamic bipartite AME model
#'
#' Describes one fit: whether it reproduces the panel it was estimated from,
#' what each block explains, whether the estimator finished cleanly, and whether
#' structure is left in the residuals. Nothing here refits a model, and nothing
#' here is a hypothesis test.
#'
#' @section What the groups answer:
#'
#' \describe{
#'   \item{`$replication`}{Simulates panels that hold the fitted structure and
#'     redraw the observation errors, and compares five panel statistics against
#'     the resulting spread. See the caveat below: these are reference
#'     intervals, not predictive intervals.}
#'   \item{`$latent_reference`}{Compares the spread of row- and column-profile
#'     correlations against panels rebuilt from the covariate and additive terms
#'     with the pair-level variation replaced by noise of the same mean and
#'     variance. A large positive `excess` says the data carry cross-dyad
#'     structure that an additive fit plus independent noise does not reproduce.}
#'   \item{`$decomposition`, `$increment`}{In-sample accounting.}
#'   \item{`$estimator`}{Convergence, and whether any limit was binding.}
#'   \item{`$residual`}{Row, column, temporal and distributional screening of
#'     the residuals.}
#' }
#'
#' @section Reading the replication intervals:
#'
#' They are descriptive conditional-reference intervals and **do not have
#' nominal coverage**. An observed value inside the interval is compatible with
#' the fitted conditional structure; a value outside indicates a discrepancy
#' worth inspecting. Every replicate holds `U` and `V` at their estimates and
#' gives `beta`, `a` and `b` only the limited dispersion adjustment the package
#' supplies, so the interval is narrow by an unquantified amount; its centre is
#' also built from the data the observed value comes from. `inside` is therefore
#' not a pass mark and never appears among the warnings, and `tail_fraction` is
#' the fraction of replicates at least as far from the reference median as the
#' observation, not a p-value.
#'
#' Two statistics are marked `by_construction` and are not independent evidence.
#' `sd.cell` is nearly matched by construction, the replicate variance being
#' `var(fitted) + sigma_eps^2` while `sigma_eps^2` is the residual mean square,
#' so the leading terms cancel. `sd.delta.cell` is an exact function of `sd.cell`
#' at both periods and the `cor.lag1` between them. Both are reported because
#' scale and movement in the outcome's own units are legible in a way the
#' remaining statistics are not.
#'
#' @section What this does not do:
#'
#' It is not a posterior predictive check. It does not test whether the
#' random-walk transition is correct: `cor.lag1` measures adjacent-panel
#' similarity, which stable covariates and a latent structure held fixed across
#' replicates both raise. It does not establish that the multiplicative term is
#' worth estimating, and `$latent_reference` is not a latent-factor significance
#' test, since unmodelled dyadic dependence, heteroskedasticity, or a
#' misspecified covariate can also produce excess. Those questions need a second
#' model fitted to the same data and belong to model comparison. Zero mass,
#' degree distributions and four-cycle counts are not checked: a Gaussian
#' outcome model does not claim to reproduce them.
#'
#' @param fit A `dynamic_ame` object from [fit_dynamic_ame()].
#' @param nsim Simulated panels per reference distribution. `0` skips both
#'   simulation groups, leaving the accounting and diagnostics, which need no
#'   simulation.
#' @param conf_level Width of the reference intervals.
#' @param min_overlap Minimum jointly observed cells for a pair of rows or
#'   columns to contribute to `sd.rowcor` or `sd.colcor`.
#' @param keep_draws Retain every replicate's statistics, which allows a
#'   different `conf_level` to be applied afterwards without simulating again.
#' @param seed Base seed. The two simulation groups are seeded independently of
#'   [bootstrap_ame()]; nothing here is shared with a bootstrap run.
#' @param n_cores Workers for the simulation groups. Both are cheap on panels of
#'   moderate size; the cost that grows is `sd.rowcor`, at order
#'   `T * N^2 * M` per replicate.
#' @param verbose Report progress.
#'
#' @return An object of class `gof_dynamic_ame`: a list of tidy data frames,
#'   `replication`, `latent_reference`, `increment`, `decomposition`,
#'   `estimator`, `residual`, `skipped` and `settings`.
#'
#' @seealso [gof_stats_ame()] for the statistics on their own,
#'   [gof_plot_ame()] to plot the result.
#'
#' @examples
#' \dontrun{
#' fit <- fit_dynamic_ame(edge_panel = my_panel, K = 2)
#' g <- gof_dynamic_ame(fit, nsim = 500)
#' g
#' gof_plot_ame(g, "replication")
#' }
#' @export
gof_dynamic_ame <- function(fit,
                            nsim = 500,
                            conf_level = 0.95,
                            min_overlap = 5L,
                            keep_draws = FALSE,
                            seed = 1,
                            n_cores = 1,
                            verbose = FALSE) {
  .validate_dynamic_ame(fit)
  if (length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between 0 and 1.", call. = FALSE)
  if (length(nsim) != 1L || nsim < 0 || nsim != as.integer(nsim))
    stop("nsim must be a non-negative integer.", call. = FALSE)
  nsim <- as.integer(nsim)
  if (nsim > 0L && nsim < 2L)
    stop("nsim must be 0 or at least 2.", call. = FALSE)

  Tn <- length(fit$years)
  periods <- as.character(fit$years)
  skipped <- list()
  note <- function(stat, per, why)
    skipped[[length(skipped) + 1L]] <<- data.frame(
      statistic = stat, period = per, reason = why, stringsAsFactors = FALSE)

  a_stats <- c("sd.cell", "sd.rowmean", "sd.colmean")
  if (Tn > 1L) {
    a_stats <- c(a_stats, "cor.lag1", "sd.delta.cell")
  } else {
    note("cor.lag1", NA_character_, "a single period has no transition")
    note("sd.delta.cell", NA_character_, "a single period has no transition")
  }
  b_stats <- c("sd.rowcor", "sd.colcor")

  replication <- NULL
  latent_reference <- NULL
  draws <- NULL

  if (nsim > 0L) {
    disp <- .dispersion_scales(fit)

    if (verbose) cat(sprintf("Conditional replication: %d panels\n", nsim))
    obs_a <- gof_stats_ame(fit$Y_list, min_overlap = min_overlap,
                           which = a_stats, periods = periods)
    ra <- .gof_reference(
      function() .gof_one_replicate(fit, disp, a_stats, min_overlap),
      obs_a, nsim, conf_level, seed, n_cores,
      two_sided = TRUE, keep_draws = keep_draws)
    replication <- ra$summary
    replication$by_construction <-
      replication$statistic %in% c("sd.cell", "sd.delta.cell")
    replication <- replication[, c("statistic", "period", "observed",
                                   "ref_median", "ref_lo", "ref_hi", "inside",
                                   "tail_fraction", "by_construction",
                                   "n_used")]

    if (verbose) cat(sprintf("Additive reference: %d panels\n", nsim))
    obs_b <- gof_stats_ame(fit$Y_list, min_overlap = min_overlap,
                           which = b_stats, periods = periods)
    parts <- .gof_null_parts(fit)
    rb <- .gof_reference(
      function() .gof_one_null(parts, fit$Y_list, b_stats, min_overlap,
                               periods),
      obs_b, nsim, conf_level, seed + 1L, n_cores,
      two_sided = FALSE, keep_draws = keep_draws)

    # the fitted model's own reconstruction of the same statistic: a point, with
    # no interval, because its dispersion would inherit A's defect
    dec <- decompose_fit(fit)
    recon <- gof_stats_ame(lapply(dec, `[[`, "fitted"),
                           min_overlap = min_overlap, which = b_stats,
                           periods = periods)
    key <- paste(rb$summary$statistic, rb$summary$period, sep = "\r")

    latent_reference <- data.frame(
      statistic = rb$summary$statistic,
      period = rb$summary$period,
      observed = rb$summary$observed,
      null_median = rb$summary$ref_median,
      null_sd = rb$summary$ref_sd,
      null_lo = rb$summary$ref_lo,
      null_hi = rb$summary$ref_hi,
      excess = (rb$summary$observed - rb$summary$ref_median) /
        rb$summary$ref_sd,
      tail_fraction = rb$summary$tail_fraction,
      in_sample_reconstruction =
        recon$value[match(key, paste(recon$statistic, recon$period,
                                     sep = "\r"))],
      n_used = rb$summary$n_used,
      stringsAsFactors = FALSE)

    for (i in which(!is.finite(replication$observed)))
      note(replication$statistic[i], replication$period[i],
           "not computable from the observed panel")
    for (i in which(!is.finite(latent_reference$observed)))
      note(latent_reference$statistic[i], latent_reference$period[i],
           sprintf("fewer than two row/column pairs reached min_overlap = %d",
                   min_overlap))

    if (keep_draws) draws <- list(replication = ra$draws,
                                  latent_reference = rb$draws)
  } else {
    note(paste(c(a_stats, b_stats), collapse = ", "), NA_character_,
         "nsim = 0: no reference distributions were simulated")
  }

  acc <- .gof_accounting(fit)
  resid <- .gof_residual(fit)

  structure(
    list(
      replication = replication,
      latent_reference = latent_reference,
      increment = acc$increment,
      decomposition = acc$decomposition,
      estimator = .gof_estimator(fit),
      residual = resid,
      draws = draws,
      skipped = if (length(skipped)) do.call(rbind, skipped) else NULL,
      settings = list(nsim = nsim, conf_level = conf_level,
                      min_overlap = min_overlap, seed = seed,
                      keep_draws = keep_draws,
                      dims = c(N = nrow(fit$a), M = nrow(fit$b), TT = Tn,
                               K = ncol(fit$U[[1]]), P = nrow(fit$beta)),
                      dispersion_blocks = c("beta", "a", "b"))
    ),
    class = "gof_dynamic_ame")
}


# ── Display ──────────────────────────────────────────────────────────────────

#' Print a goodness-of-fit result
#'
#' Warnings come first and are drawn only from the estimation diagnostics and
#' the skipped list. A statistic outside its reference interval is not among
#' them: those intervals have no nominal coverage, so falling outside one is
#' something to look at, not a failure. `multistart_spread` is likewise reported
#' but never flagged, having no calibrated threshold.
#'
#' @param x A `gof_dynamic_ame` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.gof_dynamic_ame <- function(x, ...) {
  d <- x$settings$dims
  cat("Goodness of fit for a dynamic bipartite AME model\n")
  cat(sprintf("Panel: %d x %d over %d period(s) | K = %d | P = %d\n",
              d[["N"]], d[["M"]], d[["TT"]], d[["K"]], d[["P"]]))
  cat(sprintf("Reference panels: %d | level: %.0f%%\n",
              x$settings$nsim, 100 * x$settings$conf_level))

  e <- x$estimator
  flags <- character(0)
  if (!e$converged)
    flags <- c(flags, sprintf(
      "the outer loop did not meet its convergence criterion in %d pass(es)",
      e$outer_iterations))
  # `outer_iterations` is not compared against `outer_max_iter`: under squarem
  # it totals two phases against a per-phase cap and exceeds it on a healthy
  # fit. `converged` is the flag; see `.gof_estimator()`.
  if (nzchar(e$penalties_capped))
    flags <- c(flags, sprintf("penalties held at the ceiling: %s",
                              e$penalties_capped))
  if (e$sigma_eps_floored)
    flags <- c(flags, "the residual variance sits on its floor")
  if (!is.null(x$skipped))
    flags <- c(flags, sprintf("%d statistic(s) could not be computed; see $skipped",
                              nrow(x$skipped)))

  if (length(flags)) {
    cat("\nWarnings\n")
    for (f in flags) cat("  - ", f, "\n", sep = "")
  } else {
    cat("\nNo estimation warnings.\n")
  }

  if (!is.null(x$replication)) {
    cat("\nConditional replication",
        sprintf("(%d of %d inside the reference interval)\n",
                sum(x$replication$inside, na.rm = TRUE),
                sum(!is.na(x$replication$inside))))
    cat("  These are descriptive reference intervals without nominal coverage.\n")
    cat("  Inside means compatible with the fitted conditional structure;\n")
    cat("  outside marks a discrepancy worth inspecting, not a failure.\n")
    cat("  sd.cell and sd.delta.cell are matched largely by construction.\n")
    s <- x$replication
    agg <- do.call(rbind, lapply(split(s, s$statistic), function(g)
      data.frame(statistic = g$statistic[1],
                 periods = nrow(g),
                 inside = sum(g$inside, na.rm = TRUE),
                 median_tail_fraction = stats::median(g$tail_fraction,
                                                      na.rm = TRUE),
                 by_construction = g$by_construction[1],
                 stringsAsFactors = FALSE)))
    print(agg, row.names = FALSE, digits = 3)
  }

  if (!is.null(x$latent_reference)) {
    cat("\nPair structure against the matched-moment additive reference\n")
    cat("  excess = (observed - null median) / null sd. Large and positive\n")
    cat("  means cross-dyad structure the additive reference does not produce.\n")
    s <- x$latent_reference
    agg <- do.call(rbind, lapply(split(s, s$statistic), function(g)
      data.frame(statistic = g$statistic[1],
                 periods = sum(is.finite(g$excess)),
                 min_excess = min(g$excess, na.rm = TRUE),
                 median_excess = stats::median(g$excess, na.rm = TRUE),
                 max_excess = max(g$excess, na.rm = TRUE),
                 stringsAsFactors = FALSE)))
    print(agg, row.names = FALSE, digits = 3)
  }

  if (!is.null(x$increment)) {
    cat("\nIn-sample reconstruction (one fit, nested blocks -- not three fits)\n")
    print(x$increment[, c("model", "r2", "delta_r2")], row.names = FALSE,
          digits = 3)
  }

  cat("\nEstimation\n")
  print(e, row.names = FALSE, digits = 4)
  cat("  multistart_spread measures sensitivity to the starting point only.\n")

  invisible(x)
}


#' Appearance of goodness-of-fit figures
#'
#' Collects the drawing parameters that [gof_plot_ame()] would otherwise fix
#' inside its layers -- line widths, point size and shape, colours, line types
#' and the reference band -- so they can be set from outside, the way type
#' sizes are set through `theme()`. Font sizes are deliberately not here:
#' `+ ggplot2::theme()` already does that, and a second route to the same
#' setting would only let the two disagree.
#'
#' Every argument defaults to `NULL`, meaning "use the figure's own default",
#' so `gof_style()` with no arguments reproduces the default look exactly and
#' only what is named changes. Defaults differ between figures on purpose --
#' the replication reference is neutral grey because it is this model, the
#' latent reference is blue because it is a different one -- which is why a
#' `NULL` here defers to the figure rather than fixing one value for all.
#'
#' Colours and line types are named by ROLE, not by legend text, so renaming
#' panels or legend entries never disconnects a colour from its series:
#'
#' \describe{
#'   \item{`observed`}{the observed statistic (replication, latent), and the
#'     points of the residual quantile plot}
#'   \item{`reference`}{the reference median line (replication, latent), and
#'     the 45-degree line of the residual quantile plot}
#'   \item{`covariate`, `additive`, `latent`, `residual`}{the four blocks of
#'     the decomposition figure}
#' }
#'
#' Roles a figure does not draw are ignored by it, so one style object can be
#' passed to all four figures.
#'
#' @param linewidth Width of the observed line, and of the block lines in the
#'   decomposition figure.
#' @param ref_linewidth Width of the reference line.
#' @param point_size Size of the plotted points.
#' @param point_shape Shape of the plotted points, as a ggplot2 shape code.
#' @param colours Named character vector of colours, keyed by role.
#' @param linetypes Named vector of line types, keyed by role. Numeric codes
#'   are accepted and translated to their names.
#' @param band_fill Fill of the reference band. In the latent figure it
#'   follows the `reference` colour unless set here.
#' @param band_alpha Opacity of the reference band, between 0 and 1.
#'
#' @return An object of class `gof_style`.
#'
#' @examples
#' \dontrun{
#' st <- gof_style(linewidth = 1.4, point_size = 3,
#'                 colours = c(observed = "black"),
#'                 linetypes = c(reference = "dotted"))
#' gof_plot_ame(gof, "replication", style = st) +
#'   ggplot2::theme(strip.text = ggplot2::element_text(size = 18))
#' }
#' @export
gof_style <- function(linewidth = NULL, ref_linewidth = NULL,
                      point_size = NULL, point_shape = NULL,
                      colours = NULL, linetypes = NULL,
                      band_fill = NULL, band_alpha = NULL) {

  roles <- c("observed", "reference", "covariate", "additive", "latent",
             "residual")

  positive <- function(x, nm) {
    if (!is.null(x) && (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
                        x <= 0))
      stop(nm, " must be a single positive number.", call. = FALSE)
  }
  positive(linewidth, "linewidth")
  positive(ref_linewidth, "ref_linewidth")
  positive(point_size, "point_size")

  if (!is.null(point_shape) &&
      (length(point_shape) != 1L || !(is.numeric(point_shape) ||
                                        is.character(point_shape))))
    stop("point_shape must be a single ggplot2 shape.", call. = FALSE)

  if (!is.null(band_alpha) &&
      (!is.numeric(band_alpha) || length(band_alpha) != 1L ||
       band_alpha < 0 || band_alpha > 1))
    stop("band_alpha must be a single number between 0 and 1.", call. = FALSE)

  if (!is.null(band_fill) && (!is.character(band_fill) ||
                              length(band_fill) != 1L))
    stop("band_fill must be a single colour.", call. = FALSE)

  by_role <- function(x, nm) {
    if (is.null(x)) return(NULL)
    if (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))))
      stop(nm, " must be named by role, e.g. c(observed = \"black\").",
           call. = FALSE)
    bad <- setdiff(names(x), roles)
    if (length(bad))
      stop(nm, ": unknown role(s) ", paste(bad, collapse = ", "),
           ". Roles are ", paste(roles, collapse = ", "), ".", call. = FALSE)
    x
  }
  colours <- by_role(colours, "colours")
  if (!is.null(colours) && !is.character(colours))
    stop("colours must be a character vector.", call. = FALSE)

  linetypes <- by_role(linetypes, "linetypes")
  if (is.numeric(linetypes)) {
    lt_names <- c("blank", "solid", "dashed", "dotted", "dotdash",
                  "longdash", "twodash")
    if (any(linetypes < 0 | linetypes > 6 | linetypes != round(linetypes)))
      stop("numeric linetypes must be integers from 0 to 6.", call. = FALSE)
    linetypes <- stats::setNames(lt_names[linetypes + 1L], names(linetypes))
  }

  structure(list(linewidth = linewidth, ref_linewidth = ref_linewidth,
                 point_size = point_size, point_shape = point_shape,
                 colours = colours, linetypes = linetypes,
                 band_fill = band_fill, band_alpha = band_alpha),
            class = "gof_style")
}


#' Resolve a role-keyed style vector against a figure's own defaults
#'
#' Keeps only the roles the figure draws, so a style written for another
#' figure passes through harmlessly.
#'
#' @param user Named vector from `gof_style()`, or `NULL`.
#' @param defaults Named vector of the figure's defaults.
#' @return `defaults` with the user's entries substituted.
#' @keywords internal
#' @noRd
.style_pick <- function(user, defaults) {
  if (is.null(user)) return(defaults)
  hit <- intersect(names(user), names(defaults))
  defaults[hit] <- user[hit]
  defaults
}


#' Plot a goodness-of-fit result
#'
#' Every panel shows the observed statistic as a solid line with filled points
#' against a shaded reference band and its dashed median. Colour and line type
#' both carry the distinction, so the figure survives greyscale printing and
#' colour vision deficiency.
#'
#' @param gof A `gof_dynamic_ame` object.
#' @param what Which figure: `"replication"`, `"latent"`, `"decomposition"` or
#'   `"residual"`.
#' @param statistics For `what = "replication"`, which statistics to draw.
#'   `NULL` draws every one except `sd.delta.cell`, which is an exact function
#'   of `sd.cell` at both periods and the `cor.lag1` between them and so adds a
#'   panel without adding information. Name it explicitly to see it — in the
#'   outcome's own units it is the more legible of the two, so it is worth a
#'   look when the temporal panels are the point.
#' @param labels Panel titles, as a named character vector keyed by statistic,
#'   e.g. `c(sd.rowmean = "Country heterogeneity")`. Partial: anything not
#'   named keeps its statistic's own name. Panels are titled with the raw
#'   statistic name by default, so that the figure, `gof$replication` and the
#'   documentation all use one vocabulary; rename them here when the figure has
#'   to stand on its own, where "row node" and "column node" have names of
#'   their own. Any note the panel carries is appended after the new title, the
#'   note being about the statistic rather than about what it is called.
#' @param style A [gof_style()] object setting line widths, point size and
#'   shape, colours, line types and the reference band. `NULL` draws the
#'   default look.
#' @param ... Ignored.
#' @return A `ggplot` object. Type sizes and everything else are set
#'   afterwards the usual way: `+ ggtitle()`, `+ labs()`, `+ theme()`. Panel
#'   titles are baked into the data, which is why `labels` exists; drawing
#'   parameters are baked into the layers, which is why `style` exists.
#' @export
gof_plot_ame <- function(gof, what = c("replication", "latent",
                                       "decomposition", "residual"),
                         statistics = NULL, labels = NULL, style = NULL, ...) {
  if (!inherits(gof, "gof_dynamic_ame"))
    stop("gof must be a gof_dynamic_ame object.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for gof_plot_ame().", call. = FALSE)
  what <- match.arg(what)
  # ggplot2 re-exports rlang's `.data`; taking it here keeps the pronoun
  # available without making ggplot2 a hard dependency of the package
  .data <- ggplot2::.data

  if (is.null(style)) style <- gof_style()
  if (!inherits(style, "gof_style"))
    stop("style must be made by gof_style().", call. = FALSE)
  # a style field, or the figure's own default when the field was left NULL
  sv <- function(field, default) style[[field]] %||% default

  num_period <- function(d) {
    p <- suppressWarnings(as.numeric(d$period))
    if (anyNA(p)) factor(d$period, levels = unique(d$period)) else p
  }

  if (!is.null(labels)) {
    if (!is.character(labels) || is.null(names(labels)) || anyNA(names(labels)))
      stop("labels must be a named character vector, e.g. ",
           "c(sd.rowmean = \"Country heterogeneity\").", call. = FALSE)
  }
  relabel <- function(stat) {
    out <- as.character(stat)
    if (!is.null(labels)) {
      hit <- out %in% names(labels)
      out[hit] <- unname(labels[out[hit]])
    }
    out
  }
  base <- ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0),
                   legend.position = "top")

  if (what == "replication") {
    d <- gof$replication
    if (is.null(d)) stop("no replication result; refit with nsim > 0.",
                         call. = FALSE)
    keep <- statistics %||% setdiff(unique(d$statistic), "sd.delta.cell")
    missing_stat <- setdiff(keep, unique(d$statistic))
    if (length(missing_stat))
      stop("not in this result: ", paste(missing_stat, collapse = ", "),
           ". Available: ", paste(unique(d$statistic), collapse = ", "),
           call. = FALSE)
    d <- d[d$statistic %in% keep, , drop = FALSE]
    if (!nrow(d)) stop("no statistics left to draw.", call. = FALSE)
    d$x <- num_period(d)

    # The two flagged statistics are flagged for different reasons and the
    # labels have to say which. `sd.cell` is all but guaranteed to sit inside
    # its band, the replicate variance being var(fitted) + sigma_eps^2 against
    # an observed var(Y) whose leading terms are the same two. `sd.delta.cell`
    # carries no information the other three lack, but it is free to land a
    # long way from the band -- as it does whenever cor.lag1 does. Calling
    # both "matched by construction" reads as a promise the second one does
    # not make.
    note <- c(sd.cell = "  (matched by construction)",
              sd.delta.cell = "  (determined by sd.cell and cor.lag1)")
    d$facet <- paste0(relabel(d$statistic),
                      ifelse(is.na(note[d$statistic]), "", note[d$statistic]))

    # Alphabetical order interleaves the margins with the temporal statistics.
    # Read it instead as the model is built: overall scale, then the two
    # margins, then time. Anything unrecognised keeps its own order at the end.
    d$facet <- factor(d$facet, levels = unique(
      d$facet[order(match(d$statistic, .GOF_STAT_ORDER), d$statistic)]))

    cols <- .style_pick(style$colours,
                        c(observed = "#D55E00", reference = "grey30"))
    ltys <- .style_pick(style$linetypes,
                        c(observed = "solid", reference = "dashed"))
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data$x)) +
        # the band goes through a fill scale rather than a bare colour so that
        # it appears in the legend: a grey area no key accounts for is the one
        # thing on the figure a reader cannot look up
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ref_lo,
                                          ymax = .data$ref_hi,
                                          fill = "95% reference interval"),
                             alpha = sv("band_alpha", 1)) +
        ggplot2::geom_line(ggplot2::aes(y = .data$ref_median,
                                        colour = "reference median",
                                        linetype = "reference median"),
                           linewidth = sv("ref_linewidth", 0.5)) +
        ggplot2::geom_line(ggplot2::aes(y = .data$observed, colour = "observed",
                                        linetype = "observed"),
                           linewidth = sv("linewidth", 0.5)) +
        ggplot2::geom_point(ggplot2::aes(y = .data$observed,
                                         colour = "observed"),
                            size = sv("point_size", 1.6),
                            shape = sv("point_shape", 19)) +
        ggplot2::scale_colour_manual(
          NULL, values = c(observed = unname(cols["observed"]),
                           `reference median` = unname(cols["reference"]))) +
        ggplot2::scale_linetype_manual(
          NULL, values = c(observed = unname(ltys["observed"]),
                           `reference median` = unname(ltys["reference"]))) +
        ggplot2::scale_fill_manual(
          NULL, values = c(`95% reference interval` =
                             sv("band_fill", "grey85"))) +
        ggplot2::facet_wrap(~ facet, scales = "free_y") +
        # no shared y label: each panel carries its own units, so one name over
        # all of them would have to be vague enough to be useless
        ggplot2::labs(x = "period", y = NULL) +
        base)
  }

  if (what == "latent") {
    d <- gof$latent_reference
    if (is.null(d)) stop("no latent reference; refit with nsim > 0.",
                         call. = FALSE)
    d$x <- num_period(d)
    lev <- intersect(.GOF_STAT_ORDER, d$statistic)
    d$statistic <- factor(relabel(d$statistic), levels = relabel(lev))

    cols <- .style_pick(style$colours,
                        c(observed = "#D55E00", reference = "#0072B2"))
    ltys <- .style_pick(style$linetypes,
                        c(observed = "solid", reference = "solid"))
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data$x)) +
        # The reference is a DIFFERENT MODEL, not this model's uncertainty, so
        # it gets a colour of its own rather than neutral grey, and the band
        # and its median share that colour so the two read as one object. A
        # dark dashed line inside a grey band reads instead as a second data
        # series competing with the observed one, which is what it is not.
        # The band follows the reference colour unless a fill is set, so
        # recolouring the reference keeps it one object.
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$null_lo,
                                          ymax = .data$null_hi,
                                          fill = "95% additive reference"),
                             alpha = sv("band_alpha", 0.30)) +
        ggplot2::geom_line(ggplot2::aes(y = .data$null_median,
                                        colour = "additive reference median",
                                        linetype = "additive reference median"),
                           linewidth = sv("ref_linewidth", 0.4)) +
        ggplot2::geom_line(ggplot2::aes(y = .data$observed,
                                        colour = "observed",
                                        linetype = "observed"),
                           linewidth = sv("linewidth", 0.7)) +
        ggplot2::geom_point(ggplot2::aes(y = .data$observed,
                                         colour = "observed"),
                            size = sv("point_size", 1.7),
                            shape = sv("point_shape", 19)) +
        ggplot2::scale_colour_manual(
          NULL, values = c(observed = unname(cols["observed"]),
                           `additive reference median` =
                             unname(cols["reference"]))) +
        ggplot2::scale_linetype_manual(
          NULL, values = c(observed = unname(ltys["observed"]),
                           `additive reference median` =
                             unname(ltys["reference"]))) +
        ggplot2::scale_fill_manual(
          NULL, values = c(`95% additive reference` =
                             sv("band_fill", unname(cols["reference"])))) +
        ggplot2::facet_wrap(~ statistic, scales = "free_y") +
        ggplot2::labs(x = "period", y = "spread of profile correlations") +
        base)
  }

  if (what == "decomposition") {
    d <- gof$decomposition
    if (is.null(d)) stop("no decomposition available.", call. = FALSE)
    d$x <- num_period(d)
    blocks <- sort(unique(as.character(d$block)))

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$cov_share,
                                         colour = .data$block,
                                         linetype = .data$block)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
      ggplot2::geom_line(linewidth = sv("linewidth", 0.5)) +
      ggplot2::geom_point(size = sv("point_size", 1.4),
                          shape = sv("point_shape", 19)) +
      ggplot2::labs(x = "period", y = "cov(block, Y) / var(Y)",
                    colour = NULL, linetype = NULL) +
      base

    # Only replace a scale the user actually asked to change, and fill any
    # block they left out from the palette ggplot would have used, so a
    # partial style changes the named blocks and nothing else.
    if (any(names(style$colours) %in% blocks)) {
      dflt <- stats::setNames(scales::hue_pal()(length(blocks)), blocks)
      p <- p + ggplot2::scale_colour_manual(
        NULL, values = .style_pick(style$colours, dflt))
    }
    if (any(names(style$linetypes) %in% blocks)) {
      dflt <- stats::setNames(scales::linetype_pal()(length(blocks)), blocks)
      p <- p + ggplot2::scale_linetype_manual(
        NULL, values = .style_pick(style$linetypes, dflt))
    }
    return(p)
  }

  d <- gof$residual$qq
  if (is.null(d)) stop("no residual quantiles available.", call. = FALSE)
  cols <- .style_pick(style$colours,
                      c(observed = "#D55E00", reference = "grey60"))
  ltys <- .style_pick(style$linetypes,
                      c(observed = "solid", reference = "dashed"))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = unname(cols["reference"]),
                         linetype = unname(ltys["reference"]),
                         linewidth = sv("ref_linewidth", 0.5)) +
    ggplot2::geom_point(size = sv("point_size", 1),
                        shape = sv("point_shape", 19),
                        colour = unname(cols["observed"])) +
    ggplot2::labs(x = "normal quantile", y = "standardised residual quantile") +
    base
}
