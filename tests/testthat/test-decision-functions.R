context("Decision functions")

test_that("prob_best works", {

  mat <- matrix(rnorm(13e2), 100, 13)
  pb_mat <- prob_best(mat)

  expect_equal(ncol(mat), length(pb_mat))
  expect_equal(1, sum(pb_mat))

})
