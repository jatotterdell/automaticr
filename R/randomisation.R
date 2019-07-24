#' Mass-weighted urn randomisation
#'
#' @param target_alloc The target allocation ratios
#' @param sample_size The number of allocations to generate
#' @param alpha Parameter to control imbalance between arms
#' @return A list detailing the mass-weighted-urn process.
#' @export
mass_weighted_urn_design <- function(
  target_alloc,
  sample_size,
  alpha = 4
) {
  arms <- length(target_alloc)
  prob_alloc <- target_alloc / sum(target_alloc)
  # Masses
  x <- matrix(0, sample_size + 1, arms)
  x[1, ] <- alpha * prob_alloc
  # Sample size
  n <- matrix(0, sample_size + 1, arms)
  # Random number
  y <- runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size + 1, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation Predictability
  g <- rep(0, sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)

  imbalance_cap <- sqrt(sum(((alpha - 1)*(1 - prob_alloc) + (arms - 1))^2))

  for(i in 2:(sample_size + 1)) {
    # Update allocation probabilities
    p[i - 1, ] <- pmax(alpha * prob_alloc - n[i - 1, ] + (i - 1)*prob_alloc, 0)
    p[i - 1, ] <- p[i - 1, ] / sum(p[i - 1, ])
    trt[i-1] <- findInterval(y[i - 1], c(0, cumsum(p[i - 1, ])))
    # Update sample sizes
    n[i, ] <- n[i - 1, ]
    n[i, trt[i-1]] <- n[i, trt[i-1]] + 1
    # Update urn masses
    x[i, trt[i-1]] <- x[i - 1, trt[i-1]] - 1 + prob_alloc[trt[i-1]]
    x[i, -trt[i-1]] <- x[i - 1, -trt[i-1]] + prob_alloc[-trt[i-1]]
    # Calculate imbalance
    d[i - 1] <- sqrt(sum((n[i, ] - (i - 1)*prob_alloc)^2))
    # Calculate allocation predictability
    g[i] <- d[i - 1] / alpha
  }
  return(list(
    max_imbalance_bound = imbalance_cap,
    imbalance = d,
    alloc_predict = g,
    rand_num = y,
    trt = trt,
    mass = x,
    sample_size = n,
    selection_prob = p))
}

#' Bayesian adaptive randomisation method 1
#'
#' @param p The probability each arm is best
#' @param no_alloc_thres The threshold for setting allocations to zero
#' @param fix_ctrl The fixed allocation amount to control
#' @return A vector of allocation probabilities for each arm
#' @export
brar1 <- function(p, no_alloc_thres = 0.01, fix_ctrl = NULL) {
  stopifnot(all(p >= 0))
  p[which(p < no_alloc_thres)] <- 0
  r <- sqrt(p)
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

#' Bayesian adaptive randomisation method 2
#'
#' @param p The probability each arm is best
#' @param n The current total sample size
#' @param N The total maximum sample size
#' @param no_alloc_thres The threshold for setting allocations to zero
#' @param fix_ctrl The fixed allocation amount to control
#' @return A vector of allocation probabilities for each arm
#' @export
brar2 <- function(p, n, N, no_alloc_thres = 0.01, fix_ctrl = NULL) {
  stopifnot(all(p >= 0))
  p[which(p < no_alloc_thres)] <- 0
  r <- sqrt(p ^ (n / N))
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

#' Bayesian adaptive randomisation method 3
#'
#' @param p The probability each arm is best
#' @param n The current sample size in each arm
#' @param V The current posterior variance in each arm
#' @param no_alloc_thres The threshold for setting allocations to zero
#' @param fix_ctrl The fixed allocation amount to control
#' @param min_alloc The minimum allocation to an arm
#' @return A vector of allocation probabilities for each arm
#' @export
brar3 <- function(p, n, V, no_alloc_thres = 0.01, fix_ctrl = NULL, min_alloc = 0) {
  stopifnot(all(p >= 0))
  stopifnot(all(n > 0))
  m <- length(p)

  # Set min allocation and dropping
  w <- rep(min_alloc, m)
  w[which(p < no_alloc_thres)] <- 0
  w_rem <- 1 - sum(w)
  p[which(p < no_alloc_thres)] <- 0

  r <- sqrt(p * V / n)
  w <- w + w_rem * r / sum(r)

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

#' Update allocation probabilities if some must be set to zero
#'
#' @param p The allocation probabilities
#' @param zero_ind The indices of values to be set to zero
#' @return Updated allocation probabilities
#' @export
update_alloc <- function(p, zero_ind) {
  p[zero_ind] <- 0
  p[-zero_ind] <- p[-zero_ind] / sum(p[-zero_ind])
  return(p)
}