#' Probability a column is maximum
#'
#' @param mat Matrix of MC draws
#'
#' @return Numeric vector giving probability each column is the maximum
#' @export
prob_best <- function(mat) {
  prop.table(table(factor(max.col(mat), levels = 1:ncol(mat))))
}

prob_all_equivalent <- function(mat, eps) {
  a <- rep(1/ncol(mat), ncol(mat))
  mean(apply(sweep(mat, 1, drop(mat %*% a), `-`), 1,
             function(x) all(abs(x) <= eps)))
}

prob_all_noninferior <- function(mat, colvec_best, eps) {
  dmat <- sweep(mat, 1, colvec_best, "-")
  mean(apply(dmat, 1, function(x) all(x > -eps)))
}

#' @importFrom utils combn
pairwise_diff <- function(mat, trans = I) {
  trans <- match.fun(trans)
  pair_comp <- combn(ncol(mat), 2)
  trans_mat <- trans(mat)
  pair_mat <- apply(pair_comp, 2, function(x) trans_mat[,x[1]] - trans_mat[,x[2]])
  colnames(pair_mat) <- apply(pair_comp, 2, paste, collapse = "-")
  return(pair_mat)
}

pairwise_comp <- function(pmat, eps) {
  apply(pmat, 2, function(x) mean(x > eps))
}
