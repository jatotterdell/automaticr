#' Probability a column is maximum
#'
#' @param mat Matrix of MC draws
#'
#' @return Numeric vector giving probability each column is the maximum
#' @export
prob_best <- function(mat) {
  as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
}


#' Bayesian response adaptive randomisation
#'
#' @param pbest Probability arm is best
#' @param sampsize Current sample size allocated to arm
#' @param variance Current posterior variance
#' @param no_alloc_thres Threshold for setting allocation to zero
#' @param fix_ctrl Fix allocation to control by this amount
#'
#' @return A numeric vector giving the BRAR allocation probabilities
#' @export
brar <- function(pbest, sampsize, variance, no_alloc_thres, fix_ctrl = NULL) {
  stopifnot(all(pbest >= 0))
  stopifnot(all(sampsize > 0))
  m <- length(pbest)
  r <- sqrt(pbest * variance / sampsize)
  r[which(pbest < no_alloc_thres)] <- 0
  w <- r / sum(r)

  if(is.null(fix_ctrl)) {
    return(w)
  } else {
    if(!(fix_ctrl > 0 & fix_ctrl < 1)) stop("fix_ctrl must be between 0 and 1.")
    # Re-distribute
    w[1] <- fix_ctrl
    w[-1] <- w[-1] / sum(w[-1]) * (1 - w[1])
    return(w)
  }
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
  out <- colMeans(mat[, -col] > mat[, col])
  return(out)
}