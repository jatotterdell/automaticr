#' Return list of trial parameters
#'
#' @return A \code{list} giving trial parameters.
#' @export
get_trial_params <- function() {
  trial_params <- list(
    n_arms = 13,
    n_messages = 4,
    n_timings = 3,
    sup_thres = 0.95,
    zero_alloc_thres = 0.01,
    fix_ctrl_alloc = 1/13,
    init_allocations = rep(1/13, 13)
  )
  return(trial_params)
}

#' Return intervention mapping data.frame
#'
#' @return A \code{data.frame} giving the intervention to arm mapping.
#' @export
get_intervention_map <- function() {
  automaticr:::intervention_map
}

#' Return a \code{data.frame} giving the interim schedule and seeds.
#'
#' @return A \code{data.frame} giving the interim schedule and seeds.
#' @export
get_interim_schedule <- function() {
  automaticr:::interim_schedule
}

#' Return a numeric vector giving the seeds for an interim.
#'
#' @param interim The interim number for which seeds are required.
#' @return A numeric vector giving the allocation seed then the Stan seed.
#' @export
get_interim_seeds <- function(interim) {
  automaticr:::interim_schedule[interim, c("alloc_seed", "stan_seed")]
}
