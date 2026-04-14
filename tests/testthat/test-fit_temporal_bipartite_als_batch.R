test_that("fit_temporal_bipartite_als_batch returns list of fits", {

  df1 <- expand.grid(
    year = 2020:2021,
    node_row = c("OH01", "OH02", "TX01", "TX02"),
    node_col = c("311111", "333120"),
    stringsAsFactors = FALSE
  )

  df2 <- df1

  set.seed(123)
  df1$value <- sample(1:50, nrow(df1), replace = TRUE)
  df2$value <- sample(1:50, nrow(df2), replace = TRUE)

  panel_list <- list(panel1 = df1, panel2 = df2)

  res <- fit_temporal_bipartite_als_batch(
    panel_list = panel_list,
    K = 1,
    lambda = 0.1,
    gamma = 1,
    max_iter = 50,
    tol = 1e-4,
    first_period_n_starts = 2,
    verbose = FALSE
  )

  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_true(all(names(res) == c("panel1", "panel2")))
})
