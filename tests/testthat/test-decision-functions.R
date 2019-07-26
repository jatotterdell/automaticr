context("Decision functions")

test_that("prob_best works", {

  mat <- matrix(rnorm(13e2), 100, 13)
  pb_mat <- prob_best(mat)

  expect_equal(ncol(mat), length(pb_mat))
  expect_equal(1, sum(pb_mat))

  mat <- MASS::mvrnorm(1000, mu = c(1, rep(0, 9)), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  expect_equal(all(pb_mat[1] > pb_mat[2:10]), T)

  mat <- MASS::mvrnorm(1000, mu = c(rep(0, 4), 1, rep(0, 5)), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  expect_equal(all(pb_mat[5] > c(pb_mat[1:4], pb_mat[6:10])), T)

  mat <- MASS::mvrnorm(1000, mu = c(rep(0, 9), 1), Sigma = diag(1, 10, 10))
  pb_mat <- prob_best(mat)
  expect_equal(all(pb_mat[10] > pb_mat[1:9]), T)

})

test_that("prob_rank works", {

  mat <- MASS::mvrnorm(1000, mu = c(1, 2, 3, 4), Sigma = diag(1, 4, 4))
  pr <- prob_rank(mat)

  expect_equal(all(pr[1,4] > pr[1, 1:3]), T)
  expect_equal(all(pr[2,3] > pr[2, c(1, 2, 4)]), T)
  expect_equal(all(pr[3,2] > pr[3, c(1, 3, 4)]), T)
  expect_equal(all(pr[4,1] > pr[4, c(2, 3, 4)]), T)

})

test_that("prob_all_noninferior works", {

  mat <- MASS::mvrnorm(1000, mu = c(1.2, 1, 1, 1), Sigma = diag(1, 4, 4))
  pr1 <- prob_all_noninferior(mat, mat[,1], 0.5)
  pr2 <- prob_all_noninferior(mat, mat[,1], 0.1)
  # larger noninferiority margin for 1 than 2
  expect_gt(pr1, pr2)

  mat1 <- MASS::mvrnorm(1000, mu = c(1.2, 1, 1, 1), Sigma = diag(1, 4, 4))
  mat2 <- MASS::mvrnorm(1000, mu = c(1.4, 1, 1, 1), Sigma = diag(1, 4, 4))
  pr1 <- prob_all_noninferior(mat1, mat1[,1], 0.3)
  pr2 <- prob_all_noninferior(mat2, mat2[,1], 0.3)
  # smaller mean for 1 than 2
  expect_gt(pr1, pr2)


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