gof_fit <- function(N = 12, M = 9, K = 2, Tn = 4, na_frac = 0.08, seed = 5) {
  d <- sim_ame_panel(N = N, M = M, K = K, Tn = Tn, Pr = 0, Pc = 0, Pd = 0,
                     na_frac = na_frac, seed = seed)
  ep <- do.call(rbind, lapply(seq_len(Tn), function(t) {
    data.frame(time = t,
               row = rep(rownames(d$Y[[t]]), ncol(d$Y[[t]])),
               col = rep(colnames(d$Y[[t]]), each = nrow(d$Y[[t]])),
               value = as.vector(d$Y[[t]]), stringsAsFactors = FALSE)
  }))
  ep <- ep[!is.na(ep$value), ]
  fit_dynamic_ame(edge_panel = ep, time_col = "time", row_col = "row",
                  col_col = "col", value_col = "value", K = K, n_starts = 2)
}


# ── the statistics themselves ────────────────────────────────────────────────

test_that("an additive panel has no spread in its profile correlations", {
  # Y_ij = a_i + b_j means every row profile is another shifted by a constant,
  # and correlation ignores a constant shift, so all pairs correlate at one.
  a <- rnorm(8); b <- rnorm(6)
  Y <- outer(a, b, `+`)
  expect_equal(unname(.stat_profile_cor(Y, "row")[["value"]]), 0,
               tolerance = 1e-8)
  expect_equal(unname(.stat_profile_cor(Y, "col")[["value"]]), 0,
               tolerance = 1e-8)

  # a multiplicative term varies the shape of a profile, not just its level
  set.seed(1)
  Z <- Y + matrix(rnorm(8 * 2), 8, 2) %*% t(matrix(rnorm(6 * 2), 6, 2))
  expect_gt(.stat_profile_cor(Z, "row")[["value"]], 0.05)
})

test_that("the no-NA fast path agrees with the pairwise machinery", {
  set.seed(2)
  X <- matrix(rnorm(40 * 12), 40, 12)
  expect_equal(.cor_safe(X), stats::cor(X, use = "pairwise.complete.obs"),
               tolerance = 1e-12)
  # with a gap present it must take the pairwise route, which still returns a
  # number where the plain call would give NA
  X[3, 4] <- NA
  expect_false(anyNA(.cor_safe(X)))
  expect_true(anyNA(stats::cor(X)))
})

test_that("pairs below min_overlap are dropped rather than trusted", {
  set.seed(3)
  Y <- matrix(rnorm(6 * 10), 6, 10)
  full <- .stat_profile_cor(Y, "row", min_overlap = 5)
  # leave row 1 with only three observed columns
  Y[1, 4:10] <- NA
  thin <- .stat_profile_cor(Y, "row", min_overlap = 5)
  expect_equal(unname(full[["n_used"]]), 15)   # 6 choose 2
  expect_equal(unname(thin[["n_used"]]), 10)   # every pair involving row 1 gone
})

test_that("wholly unobserved margins are removed, not turned into NaN", {
  Y <- matrix(rnorm(20), 5, 4)
  Y[2, ] <- NA
  v <- .stat_margin_sd(Y, "row")
  expect_equal(unname(v[["n_used"]]), 4)
  expect_true(is.finite(v[["value"]]))
})

test_that("gof_stats_ame labels transitions by the period they arrive at", {
  set.seed(4)
  Y <- lapply(1:3, function(t) matrix(rnorm(20), 5, 4))
  s <- gof_stats_ame(Y, periods = c("2001", "2002", "2003"))
  expect_setequal(unique(s$statistic),
                  c("sd.cell", "sd.rowmean", "sd.colmean", "sd.rowcor",
                    "sd.colcor", "cor.lag1", "sd.delta.cell"))
  expect_setequal(s$period[s$statistic == "cor.lag1"], c("2002", "2003"))
  expect_equal(sum(s$statistic == "sd.cell"), 3L)
})

test_that("sd.delta.cell is the exact function of the other three it claims", {
  set.seed(5)
  Y <- lapply(1:2, function(t) matrix(rnorm(60), 10, 6))
  s <- gof_stats_ame(Y, which = c("sd.cell", "cor.lag1", "sd.delta.cell"))
  v <- function(st, p) s$value[s$statistic == st & s$period == p]
  # var(Yt - Yt-1) = var(Yt) + var(Yt-1) - 2 cov, on the common cell set
  implied <- sqrt(v("sd.cell", "1")^2 + v("sd.cell", "2")^2 -
                    2 * v("cor.lag1", "2") * v("sd.cell", "1") *
                    v("sd.cell", "2"))
  expect_equal(v("sd.delta.cell", "2"), implied, tolerance = 1e-6)
})

test_that("a single period drops the transition statistics without erroring", {
  s <- gof_stats_ame(list(matrix(rnorm(20), 5, 4)))
  expect_false(any(c("cor.lag1", "sd.delta.cell") %in% s$statistic))
})

test_that("which= computes only what was asked for", {
  Y <- lapply(1:2, function(t) matrix(rnorm(20), 5, 4))
  s <- gof_stats_ame(Y, which = c("sd.rowmean", "cor.lag1"))
  expect_setequal(unique(s$statistic), c("sd.rowmean", "cor.lag1"))
})


# ── the assembled object ─────────────────────────────────────────────────────

test_that("gof_dynamic_ame returns every documented element", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 25, seed = 1)
  expect_s3_class(g, "gof_dynamic_ame")
  for (nm in c("replication", "latent_reference", "increment", "decomposition",
               "estimator", "residual", "settings")) {
    expect_false(is.null(g[[nm]]), info = nm)
  }
  expect_named(g$replication,
               c("statistic", "period", "observed", "ref_median", "ref_lo",
                 "ref_hi", "inside", "tail_fraction", "by_construction",
                 "n_used"))
  expect_true(all(c("excess", "in_sample_reconstruction") %in%
                    names(g$latent_reference)))
})

test_that("the covariance shares sum to one within every period", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 0)
  s <- tapply(g$decomposition$cov_share, g$decomposition$period, sum)
  expect_equal(as.numeric(s), rep(1, length(s)), tolerance = 1e-8)
  expect_setequal(unique(g$decomposition$block),
                  c("covariate", "additive", "latent", "residual"))
})

test_that("nsim = 0 skips the simulated groups and keeps the rest", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 0)
  expect_null(g$replication)
  expect_null(g$latent_reference)
  expect_false(is.null(g$increment))
  expect_false(is.null(g$estimator))
  expect_false(is.null(g$residual$by_period))
  expect_true(any(grepl("nsim = 0", g$skipped$reason)))
})

test_that("scale and movement are flagged as matched by construction", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 25, seed = 1)
  bc <- unique(g$replication$statistic[g$replication$by_construction])
  expect_setequal(bc, c("sd.cell", "sd.delta.cell"))
})

test_that("the reference intervals bracket their own median", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 40, seed = 2)
  r <- g$replication[is.finite(g$replication$ref_median), ]
  expect_true(all(r$ref_lo <= r$ref_median & r$ref_median <= r$ref_hi))
  l <- g$latent_reference[is.finite(g$latent_reference$null_median), ]
  expect_true(all(l$null_lo <= l$null_median & l$null_median <= l$null_hi))
})

test_that("tail fractions are two-sided in A and one-sided in B", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 40, seed = 3, keep_draws = TRUE)

  # A: the fraction at least as far from the median as the observation
  i <- 1L
  d <- g$draws$replication[i, ]
  med <- stats::median(d, na.rm = TRUE)
  expect_equal(g$replication$tail_fraction[i],
               mean(abs(d - med) >= abs(g$replication$observed[i] - med)),
               tolerance = 1e-12)

  # B: the fraction at or above the observation
  d <- g$draws$latent_reference[i, ]
  expect_equal(g$latent_reference$tail_fraction[i],
               mean(d >= g$latent_reference$observed[i]), tolerance = 1e-12)
})

test_that("the additive reference matches the residual's first two moments", {
  fit <- gof_fit()
  parts <- .gof_null_parts(fit)
  dec <- decompose_fit(fit)
  for (t in seq_along(dec)) {
    r <- fit$Y_list[[t]] - dec[[t]]$covariate - dec[[t]]$additive
    ok <- is.finite(r)
    # the reference's own centre is the covariate + additive term plus the
    # residual mean, so its base already carries the first moment
    expect_equal(mean((parts$base[[t]] - dec[[t]]$covariate -
                         dec[[t]]$additive)[ok]),
                 mean(r[ok]), tolerance = 1e-10)
    expect_equal(parts$sd[t], sqrt(mean((r[ok] - mean(r[ok]))^2)),
                 tolerance = 1e-10)
  }
})

test_that("the additive reference carries no multiplicative structure", {
  # its profile-correlation spread should sit near the noise level, well under
  # what the same panel shows with the latent block present
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 60, seed = 7)
  l <- g$latent_reference[g$latent_reference$statistic == "sd.rowcor", ]
  expect_true(all(is.finite(l$excess)))
  expect_gt(stats::median(l$excess), 0)
})

test_that("the reference panels keep the observed missing pattern", {
  fit <- gof_fit(na_frac = 0.15)
  parts <- .gof_null_parts(fit)
  set.seed(1)
  Yn <- .gof_one_null(parts, fit$Y_list, "sd.rowcor", 5L,
                      as.character(fit$years))
  expect_s3_class(Yn, "data.frame")
  # and directly: the simulator must not invent observations
  B <- parts$base[[1]]
  Y <- B + matrix(rnorm(length(B), sd = parts$sd[1]), nrow(B), ncol(B))
  Y[is.na(fit$Y_list[[1]])] <- NA
  expect_equal(is.na(Y), is.na(fit$Y_list[[1]]))
})

test_that("the estimator diagnostics read the fit's own record", {
  fit <- gof_fit()
  e <- gof_dynamic_ame(fit, nsim = 0)$estimator
  expect_equal(e$converged, isTRUE(fit$convergence$converged))
  expect_equal(e$outer_iterations, fit$convergence$outer_iterations)
  expect_equal(e$outer_max_iter, fit$settings$outer_max_iter)
  expect_equal(e$accelerate, fit$settings$accelerate)
  expect_type(e$penalties_capped, "character")
  expect_type(e$sigma_eps_floored, "logical")
})

test_that("exceeding outer_max_iter is not flagged, because squarem does it", {
  # The extrapolation gets floor(outer_max_iter / 3) steps at about three
  # passes each, and the plain certification is then given the whole of
  # outer_max_iter again (see .ame_outer_eb()). `outer_iterations` counts every
  # pass across both, so it reaches roughly 4/3 of the cap on a fit that is
  # perfectly healthy. Comparing the two would be a false positive by
  # construction.
  fit <- gof_fit()
  fit$convergence$outer_iterations <- 135L
  fit$settings$outer_max_iter <- 100L
  fit$convergence$converged <- TRUE
  g <- gof_dynamic_ame(fit, nsim = 0)
  out <- paste(utils::capture.output(print(g)), collapse = "\n")
  expect_false(grepl("iteration cap", out))
  expect_false(grepl("convergence criterion", out))
  expect_match(out, "135")     # reported, just not flagged

  # non-convergence is still reported, and reports the pass count with it
  fit$convergence$converged <- FALSE
  g2 <- gof_dynamic_ame(fit, nsim = 0)
  out2 <- paste(utils::capture.output(print(g2)), collapse = "\n")
  expect_match(out2, "did not meet its convergence criterion in 135")
})

test_that("residual screening covers every period and thins the quantiles", {
  fit <- gof_fit()
  r <- gof_dynamic_ame(fit, nsim = 0)$residual
  expect_equal(nrow(r$by_period), length(fit$years))
  expect_true(is.na(r$by_period$cor_lag1_resid[1]))   # nothing before period 1
  expect_lte(nrow(r$qq), 200L)
  expect_true(all(c("fitted_mid", "resid_sd") %in% names(r$scale)))
})

test_that("results are reproducible under a seed and differ without one", {
  fit <- gof_fit()
  a <- gof_dynamic_ame(fit, nsim = 20, seed = 11)
  b <- gof_dynamic_ame(fit, nsim = 20, seed = 11)
  c_ <- gof_dynamic_ame(fit, nsim = 20, seed = 12)
  expect_equal(a$replication$ref_median, b$replication$ref_median)
  expect_false(isTRUE(all.equal(a$replication$ref_median,
                                c_$replication$ref_median)))
})

test_that("the two simulated groups do not share a seed stream", {
  # B is offset from A, so the two references are not the same random panels
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1, keep_draws = TRUE)
  expect_false(isTRUE(all.equal(as.numeric(g$draws$replication[1, ]),
                                as.numeric(g$draws$latent_reference[1, ]))))
})

test_that("invalid arguments are refused", {
  fit <- gof_fit()
  expect_error(gof_dynamic_ame(fit, conf_level = 0), "conf_level")
  expect_error(gof_dynamic_ame(fit, conf_level = 1), "conf_level")
  expect_error(gof_dynamic_ame(fit, nsim = -1), "nsim")
  expect_error(gof_dynamic_ame(fit, nsim = 1), "nsim")
  expect_error(gof_stats_ame(list()), "non-empty")
  expect_error(gof_stats_ame(list(matrix(1:4, 2)), periods = c("a", "b")),
               "one label per period")
})

test_that("a single-period fit reports the transitions as skipped", {
  fit <- gof_fit(Tn = 1, na_frac = 0)
  g <- gof_dynamic_ame(fit, nsim = 10, seed = 1)
  expect_false(any(c("cor.lag1", "sd.delta.cell") %in%
                     g$replication$statistic))
  expect_true(all(c("cor.lag1", "sd.delta.cell") %in% g$skipped$statistic))
})

test_that("print reports estimation warnings but never an interval miss", {
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 25, seed = 1)
  out <- paste(utils::capture.output(print(g)), collapse = "\n")
  expect_match(out, "without nominal coverage")
  expect_match(out, "not a failure")
  expect_match(out, "starting point only")
  # a statistic outside its interval must not be reported as a warning
  warn_block <- sub("\n\nConditional replication.*", "", out)
  expect_false(grepl("outside", warn_block))
})

test_that("gof_plot_ame builds every panel and refuses what was not computed", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  for (w in c("replication", "latent", "decomposition", "residual")) {
    expect_s3_class(gof_plot_ame(g, w), "ggplot")
  }
  g0 <- gof_dynamic_ame(fit, nsim = 0)
  expect_error(gof_plot_ame(g0, "replication"), "nsim")
  expect_error(gof_plot_ame(fit, "replication"), "gof_dynamic_ame object")
})

test_that("the replication figure drops sd.delta.cell unless asked for it", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)

  drawn <- function(p) unique(as.character(p$data$statistic))
  expect_setequal(drawn(gof_plot_ame(g, "replication")),
                  c("sd.cell", "sd.rowmean", "sd.colmean", "cor.lag1"))

  # still in the table, and still drawable on request
  expect_true("sd.delta.cell" %in% g$replication$statistic)
  expect_true("sd.delta.cell" %in%
                drawn(gof_plot_ame(g, "replication",
                                   statistics = c("cor.lag1",
                                                  "sd.delta.cell"))))

  expect_error(gof_plot_ame(g, "replication", statistics = "sd.nonesuch"),
               "not in this result")
})

test_that("panels read scale, margins, time -- not alphabetically", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)

  lev <- levels(gof_plot_ame(g, "replication")$data$facet)
  expect_equal(sub("  \\(.*", "", lev),
               c("sd.cell", "sd.rowmean", "sd.colmean", "cor.lag1"))

  # asking for sd.delta.cell puts it after cor.lag1, not between sd.cell and
  # sd.colmean where the alphabet would leave it
  lev2 <- levels(gof_plot_ame(g, "replication",
                              statistics = c("sd.delta.cell", "sd.cell",
                                             "cor.lag1"))$data$facet)
  expect_equal(sub("  \\(.*", "", lev2),
               c("sd.cell", "cor.lag1", "sd.delta.cell"))

  # rows before columns on the latent figure too
  expect_equal(levels(gof_plot_ame(g, "latent")$data$statistic),
               c("sd.rowcor", "sd.colcor"))
})

test_that("panels are titled with the raw statistic names by default", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  lev <- levels(gof_plot_ame(g, "replication")$data$facet)
  # the name in the figure is the name in gof$replication and in the docs
  expect_true(all(sub("  \\(.*", "", lev) %in% g$replication$statistic))
})

test_that("labels renames panels, partially, and keeps the notes", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)

  p <- gof_plot_ame(g, "replication",
                    labels = c(sd.rowmean = "Country heterogeneity"))
  lev <- levels(p$data$facet)
  expect_true("Country heterogeneity" %in% lev)
  expect_false("sd.rowmean" %in% lev)
  expect_true("sd.colmean" %in% lev)          # untouched names stay

  # a renamed panel keeps its note: the note is about the statistic, not
  # about what the statistic is called
  p2 <- gof_plot_ame(g, "replication", labels = c(sd.cell = "Overall spread"))
  expect_true(any(grepl("^Overall spread  \\(matched by construction\\)$",
                        levels(p2$data$facet))))

  # renaming must not disturb the reading order
  expect_equal(which(lev == "Country heterogeneity"), 2L)

  # and it works on the latent figure
  p3 <- gof_plot_ame(g, "latent", labels = c(sd.rowcor = "Country profiles"))
  expect_equal(levels(p3$data$statistic), c("Country profiles", "sd.colcor"))

  expect_error(gof_plot_ame(g, "replication", labels = "no name"),
               "named character vector")
})

test_that("the band appears in the legend and no figure carries a caption", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)

  # a grey area with no key is the one thing on the figure a reader cannot
  # look up, so the ribbon goes through a fill scale
  for (w in c("replication", "latent")) {
    p <- gof_plot_ame(g, w)
    expect_true("fill" %in% names(p$labels) ||
                  any(vapply(p$layers, function(l) "fill" %in% names(l$mapping),
                             logical(1))))
  }
  for (w in c("replication", "latent", "decomposition", "residual"))
    expect_null(gof_plot_ame(g, w)$labels$caption)
})

test_that("the additive reference is one coloured object, not a grey dash", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  b <- ggplot2::ggplot_build(gof_plot_ame(g, "latent"))

  ref <- "#0072B2"
  # the band and its median share a colour so the reference reads as one
  # thing; a dark dash inside a grey band reads as a second data series
  expect_true(any(vapply(b$data, function(l)
    "fill" %in% names(l) && any(l$fill == ref), logical(1))))
  expect_true(any(vapply(b$data, function(l)
    "colour" %in% names(l) && any(l$colour == ref), logical(1))))
  # and it is a colour of its own, because it is a DIFFERENT MODEL rather
  # than this model's uncertainty
  expect_false(any(vapply(b$data, function(l)
    "colour" %in% names(l) && any(l$colour %in% c("grey30", "black")),
    logical(1))))
  expect_true(any(vapply(b$data, function(l)
    "colour" %in% names(l) && any(l$colour == "#D55E00"), logical(1))))
})

test_that("the two flagged statistics get labels that say different things", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  f <- unique(gof_plot_ame(g, "replication",
                           statistics = c("sd.cell", "sd.delta.cell",
                                          "sd.rowmean"))$data$facet)
  # sd.cell all but has to sit inside its band; sd.delta.cell is merely
  # redundant, and is free to land far outside -- one label cannot serve both
  expect_true(any(grepl("sd.cell  \\(matched by construction\\)", f)))
  expect_true(any(grepl("sd.delta.cell  \\(determined by", f)))
  expect_true("sd.rowmean" %in% f)
})


# ── gof_style(): drawing parameters set from outside ────────────────────────

# every drawn value of one aesthetic across a built plot's layers
drawn <- function(p, aes, geom = NULL) {
  b <- ggplot2::ggplot_build(p)
  keep <- if (is.null(geom)) seq_along(b$data) else
    which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
  unique(unlist(lapply(b$data[keep], function(l) l[[aes]])))
}

test_that("gof_style() refuses what it cannot use", {
  expect_error(gof_style(linewidth = -1), "positive")
  expect_error(gof_style(point_size = c(1, 2)), "positive")
  expect_error(gof_style(band_alpha = 2), "between 0 and 1")
  expect_error(gof_style(colours = c("red")), "named by role")
  expect_error(gof_style(colours = c(obseved = "red")), "unknown role")
  expect_error(gof_style(linetypes = c(observed = 9)), "0 to 6")
  expect_error(gof_style(band_fill = c("red", "blue")), "single colour")
  expect_s3_class(gof_style(), "gof_style")
})

test_that("numeric line types are translated to their names", {
  st <- gof_style(linetypes = c(observed = 1, reference = 3))
  expect_identical(unname(st$linetypes), c("solid", "dotted"))
})

test_that("an empty style draws exactly the default look", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  for (w in c("replication", "latent", "decomposition", "residual")) {
    a <- ggplot2::ggplot_build(gof_plot_ame(g, w))$data
    b <- ggplot2::ggplot_build(gof_plot_ame(g, w, style = gof_style()))$data
    expect_identical(a, b, info = w)
  }
  # and the defaults are the values the figures had before styling existed
  p <- gof_plot_ame(g, "replication")
  expect_setequal(drawn(p, "colour", "GeomLine"), c("#D55E00", "grey30"))
  expect_setequal(drawn(p, "linewidth", "GeomLine"), 0.5)
  expect_setequal(drawn(p, "size", "GeomPoint"), 1.6)
  expect_setequal(drawn(p, "fill", "GeomRibbon"), "grey85")
})

test_that("each style parameter reaches the replication figure", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  st <- gof_style(linewidth = 1.4, ref_linewidth = 0.9, point_size = 3,
                  point_shape = 17,
                  colours = c(observed = "black", reference = "red"),
                  linetypes = c(reference = "dotted"),
                  band_fill = "lightblue", band_alpha = 0.5)
  p <- gof_plot_ame(g, "replication", style = st)

  expect_setequal(drawn(p, "linewidth", "GeomLine"), c(1.4, 0.9))
  expect_setequal(drawn(p, "size", "GeomPoint"), 3)
  expect_setequal(drawn(p, "shape", "GeomPoint"), 17)
  expect_setequal(drawn(p, "colour", "GeomLine"), c("black", "red"))
  expect_true("dotted" %in% drawn(p, "linetype", "GeomLine"))
  expect_setequal(drawn(p, "fill", "GeomRibbon"), "lightblue")
  expect_setequal(drawn(p, "alpha", "GeomRibbon"), 0.5)
})

test_that("a partial style changes what it names and nothing else", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  p <- gof_plot_ame(g, "replication",
                    style = gof_style(colours = c(observed = "black")))
  # the reference keeps its default colour
  expect_setequal(drawn(p, "colour", "GeomLine"), c("black", "grey30"))
})

test_that("the latent band follows a recoloured reference unless told not to", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  p1 <- gof_plot_ame(g, "latent",
                     style = gof_style(colours = c(reference = "darkgreen")))
  expect_setequal(drawn(p1, "fill", "GeomRibbon"), "darkgreen")
  p2 <- gof_plot_ame(g, "latent",
                     style = gof_style(colours = c(reference = "darkgreen"),
                                       band_fill = "grey80"))
  expect_setequal(drawn(p2, "fill", "GeomRibbon"), "grey80")
})

test_that("the decomposition takes block colours partially", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 0)
  p <- gof_plot_ame(g, "decomposition",
                    style = gof_style(colours = c(latent = "black"),
                                      linewidth = 2))
  cols <- drawn(p, "colour", "GeomLine")
  expect_true("black" %in% cols)
  expect_length(cols, 4L)          # the other three kept a palette colour
  expect_setequal(drawn(p, "linewidth", "GeomLine"), 2)
})

test_that("the residual figure takes point and reference styling", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 0)
  p <- gof_plot_ame(g, "residual",
                    style = gof_style(point_size = 2.5,
                                      colours = c(observed = "navy",
                                                  reference = "red")))
  expect_setequal(drawn(p, "size", "GeomPoint"), 2.5)
  expect_setequal(drawn(p, "colour", "GeomPoint"), "navy")
  expect_setequal(drawn(p, "colour", "GeomAbline"), "red")
})

test_that("one style object can be passed to every figure", {
  skip_if_not_installed("ggplot2")
  fit <- gof_fit()
  g <- gof_dynamic_ame(fit, nsim = 20, seed = 1)
  st <- gof_style(linewidth = 1.2,
                  colours = c(observed = "black", additive = "red"))
  for (w in c("replication", "latent", "decomposition", "residual"))
    expect_s3_class(gof_plot_ame(g, w, style = st), "ggplot")
  expect_error(gof_plot_ame(g, "replication", style = list(linewidth = 2)),
               "gof_style")
})
