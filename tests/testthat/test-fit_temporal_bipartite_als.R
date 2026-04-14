test_that("fit_temporal_bipartite_als returns expected outputs", {
  df <- expand.grid(
    year = 2020:2021,
    node_row = c("OH01", "OH02", "TX01", "TX02"),
    node_col = c("311111", "333120"),
    stringsAsFactors = FALSE
  )

  set.seed(123)
  df$value <- sample(1:50, nrow(df), replace = TRUE)

  res <- fit_temporal_bipartite_als(
    edge_panel = df,
    K = 1,
    lambda = 0.1,
    gamma = 1,
    max_iter = 50,
    tol = 1e-4,
    first_period_n_starts = 2,
    verbose = FALSE
  )

  expect_type(res, "list")
  expect_true(all(c("results", "node_df", "coef_df", "params") %in% names(res)))
  expect_true(is.data.frame(res$node_df))
})
