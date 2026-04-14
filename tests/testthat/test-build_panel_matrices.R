test_that("build_panel_matrices returns expected structure", {
  df <- data.frame(
    year = c(2020, 2020, 2021, 2021),
    node_row = c("OH01", "OH01", "OH01", "OH01"),
    node_col = c("311111", "333120", "311111", "333120"),
    value = c(10, 20, 15, 25)
  )

  res <- build_panel_matrices(df)

  expect_type(res, "list")
  expect_true(all(c("years", "Y_list", "row_ids", "col_ids") %in% names(res)))
  expect_equal(length(res$years), 2)
  expect_equal(length(res$Y_list), 2)

  expect_true(all(res$row_ids == "row_OH01"))
  expect_true(all(res$col_ids == c("col_311111", "col_333120")))
})
