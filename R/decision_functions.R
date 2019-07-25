#' Probability a column is maximum
#'
#' @param mat Matrix of MC draws
#'
#' @return Numeric vector giving probability each column is the maximum
#' @export
prob_best <- function(mat) {
  as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
}


#' Probability best of linear transformation
#'
#' For a supplied linear transform, A, of MC draws, determine
#' the probability best.
#'
#' @param mat Matrix of MC draws
#' @param A Transformation matrix
#' @return Numeric vector giving probability best each transformed column is maximum
#' @export
prob_best_lin <- function(mat, A) {
  prob_best(mat %*% t(A))
}


#' Probability of each rank for each arm.
#'
#' Rank each row and return the probability that each arm is ranked at each value
#'
#' @param mat Matrix of MC draws
#' @return Numeric square matrix giving probability each arm is ranked in each place.
#'   Rows are ranks and columns are arms. For example, the first row in the matrix gives the probability
#'   each arm is the best. The last row gives the probability each arm is the worst, etc.
#' @export
prob_rank <- function(mat) {
  k <- ncol(mat)
  apply(apply(apply(mat, 1, function(x) rank(-x)), 1, function(z) table(factor(z, levels = 1:k))), 2, prop.table)
}


#' Probability all columns are non-inferior to best column
#'
#' @param mat Matrix of MC draws
#' @param colvec_best A vector of draws from the best
#' @param eps The reference value for non-inferiority
#' @return A real number giving the probability that all arms
#'   are non-inferior to colvec_best.
#' @export
prob_all_noninferior <- function(mat, colvec_best, eps) {
  dmat <- sweep(mat, 1, colvec_best, "-")
  mean(apply(dmat, 1, function(x) all(x > -eps)))
}

#' Probability all columns are equivalent - spherical
#'
#' Spherical test that all columns are within bound delta of global mean
#'
#' @param mat Matrix of MC draws
#' @param delta Equivalence bound
#' @return The probability of spherical equivalence
#' @export
prob_equi <- function(mat, delta, ...) {
  dmat <- pairwise_diff(mat, ...)
  mean(apply(dmat, 1, crossprod) < delta^2)
}


#' Calculate matrix of pairwise differences of MC draws
#'
#' @param mat The matrix of draws
#' @param trans A transformation function, default is identity \code{I()}.
#' @return A matrix of pairwise differences
#' @export
#' @importFrom utils combn
pairwise_diff <- function(mat, trans = I) {
  trans <- match.fun(trans)
  pair_comp <- combn(ncol(mat), 2)
  trans_mat <- trans(mat)
  pair_mat <- apply(pair_comp, 2, function(x) trans_mat[,x[1]] - trans_mat[,x[2]])
  colnames(pair_mat) <- apply(pair_comp, 2, paste, collapse = "-")
  return(pair_mat)
}


#' Make a pairwise comparison between parameter draws
#'
#' @param mat Matrix of parameter draws
#' @param eps Reference value, default is \code{0}.
#' @return Vector of probabilities that pairwise difference is greater than \code{eps}.
#' @export
pairwise_comp <- function(mat, eps = 0) {
  pmat <- pairwise_diff(mat)
  apply(pmat, 2, function(x) mean(x > eps))
}


#' Calcualte probability that each column is superior to a reference column.
#'
#' @param mat Matrix of parameter draws
#' @param col The reference column integer
#' @param eps The reference value for superiority, default is \code{0}
#' @return A vector of probabilities
#' @export
col_comp <- function(mat, col, eps = 0) {
  out <- colMeans(mat[, -col] > mat[, col] + eps)
  return(out)
}
