context("Decision functions")

test_that("prob_best works", {

  mat <- matrix(rnorm(13e2), 100, 13)
  pb_mat <- prob_best(mat)

  expect_equal(ncol(mat), length(pb_mat))
  expect_equal(1, sum(pb_mat))

  mat <- MASS::mvrnorm(1000, mu = c(1, rep(0, 9)), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  expect_equal(all(rep(pb_mat[1], 9) > pb_mat[2:10]), T)

  mat <- MASS::mvrnorm(1000, mu = c(rep(0, 4), 1, rep(0, 5)), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  expect_equal(all(rep(pb_mat[5], 9) > c(pb_mat[1:4], pb_mat[6:10])), T)

  mat <- MASS::mvrnorm(1000, mu = c(rep(0, 9), 1), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  expect_equal(all(rep(pb_mat[10], 9) > pb_mat[1:9]), T)

})



test_that("brar", {

  # assess brar updates when there should be a clear winner and loser
  mat <- MASS::mvrnorm(1000, mu = c(1, 3, rep(1, 7), 0), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  pb_var <- diag(var(mat))
  ss <- rep(100, 10)

  new_probs <- brar(pb_mat, ss, pb_var, 0.01, 1/10)

  expect_equal(sum(new_probs), 1, tolerance = 0.001)
  expect_equal(new_probs[10], 0, tolerance = 0.001)
  expect_equal(new_probs[1], 0.1, tolerance = 0.001)
  expect_equal(all(rep(new_probs[2], 9) > c(new_probs[1], new_probs[3:10])))


  # assess brar updates all should be approx the same
  mat <- MASS::mvrnorm(1000, mu = rep(1, 10), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  pb_var <- diag(var(mat))
  ss <- rep(100, 10)

  new_probs <- brar(pb_mat, ss, pb_var, 0.01, NULL)
  expect_equal(abs(new_probs - 0.1), rep(0, 10), tolerance = 0.01)


})