#' Intialise the interim log
#'
#' Intialises the interim log based on \code{trial_params} setting
#' dates, number of treatment arms and initial randomisation
#' probabilities.
#' Should only be called on first build.
#'
#' @param file The path for the interim log file
#'
#' @return NULL
#' @export
#'
#' @importFrom dplyr %>%
init_interim_log <- function(file) {
  if(file.exists(file)) {
    stop("File already exists.")
  }
  warning(paste("Initialising", file))
  pars <- automaticr::get_trial_params()
  tb <- tibble::tibble(
    interim_date = Sys.Date(),
    interim_num  = 0,
    n_analysed    = 0,
    arm          = sprintf("%02d", 1:pars$n_arms),
    alloc_prob   = pars$init_allocations,
    is_alloc     = TRUE
  ) %>%
    tidyr::gather(variable, value, -(interim_date:arm)) %>%
    tidyr::unite(temp, variable, arm) %>%
    tidyr::spread(temp, value)

  readr::write_csv(tb, file)

  return(invisible(NULL))
}


#' Read raw data
#'
#' @param file The path to the database data file.
#'
#' @return A \code{tibble}.
#' @export
read_raw_data <- function(file) {
  if(!file.exists(file)) {
    warning(paste(file, "file not found."))
    quit(status = 0)
  } else {
    dat <- readr::read_csv(
      file,
      na = c("null", "NA", ""))
  }
  if(nrow(dat) == 0) {
    warning(paste(file, "has no records."))
    quit(status = 0)
  }
  return(dat)
}


#' Process raw data
#'
#' Needs to check raw data for any coding issues etc.
#'
#' @param raw_dat The raw data.
#' @param ref_date The reference date from which to calculate how much time has passed since the vaccine due date.
#'   This should probably be the date on which the function was called, i.e. \code{Sys.Date()}.
#'
#' @return A \code{tibble}.
#' @export
process_raw_data <- function(raw_dat, ref_date) {
  # raw_dat = dat

  if(!lubridate::is.Date(ref_date))
    stop("ref_date is not a date.")
  if(nrow(raw_dat) == 0)
    stop("raw_dat has zero rows.")

  dat <- dplyr::mutate(raw_dat,
    current_eligible_vaccination =
    factor(current_eligible_vaccination,
          levels = paste0(c(2,4,6,12,18,48), "m")),
    randomisation_outcome =
    factor(randomisation_outcome,
          levels = 1:13, labels = sprintf("%02d", 1:13)),
    date_vaccine_due = lubridate::as_date(date_vaccine_due),
    date_of_vaccination_administration = lubridate::as_date(date_of_vaccination_administration),
    time_to_vax = difftime(date_of_vaccination_administration, date_vaccine_due, units = "days"),
    time_since_due = difftime(ref_date, date_vaccine_due, units = "days"),
    vax_past_due = time_since_due > 28,
    on_time = dplyr::case_when(
    time_to_vax <= 28 ~ 1,
    time_to_vax > 28 ~ 0,
    is.na(time_to_vax) & vax_past_due ~ 0,
    TRUE ~ NA_real_)
  )

  return(dat)
}


#' Return index records (i.e. first record) for \code{parent_id}
#'
#' Primary analysis will be based on only index records for \code{parent_id}.
#'
#' @param dat The processed raw data; a \code{tibble}.
#'
#' @return A \code{tibble} subset of \code{dat}
#'   with one row per \code{parent_id} corresponding to the earliest
#'   due vaccination.
#' @export
get_index_data <- function(dat) {
  dat %>%
    dplyr::group_by(parent_id) %>%
    dplyr::top_n(1, date_vaccine_due)
}


#' Aggregate the index data set
#'
#' For modelling purposes, the data should be aggregated according to
#' \code{randomisation_outcome} and \code{current_eligible_vaccination}
#' which will be covariates in the model.
#'
#' @param dat The index data set as output from \code{get_index_data}.
#'
#' @return A \code{tibble} of aggregated data from \code{dat}.
#' @export
aggregate_data <- function(dat) {
  # agg_dat <- dplyr::ungroup(dplyr::arrange(
  #   dplyr::summarise(
  #     dplyr::group_by(
  #       dplyr::filter(dat, !is.na(on_time)),
  #       randomisation_outcome, current_eligible_vaccination),
  #     y = sum(on_time), trials = dplyr::n()),
  # randomisation_outcome, current_eligible_vaccination))

  agg_dat <- dat %>% dplyr::filter(!is.na(on_time)) %>%
    dplyr::group_by(randomisation_outcome, current_eligible_vaccination) %>%
    dplyr::summarise(y = sum(on_time), trials = dplyr::n()) %>%
    dplyr::arrange(randomisation_outcome, current_eligible_vaccination) %>%
    dplyr::ungroup()
  dplyr::left_join(agg_dat, automaticr::get_intervention_map(), by = "randomisation_outcome")
}

#' Eigen-decompose a constraining prior to create identifiable parameters
#'
#' @param X The design matrix
#' @param S The constraint matrix
#' @return The constrained design matrix
#' @export
constrain_design <- function(X, S) {
  d <- ncol(X)
  r <- Matrix::rankMatrix(S)
  return(X %*% eigen(S)$vector[, 1:r])
}


#' Make data list for input into Stan model
#'
#' @param agg_dat Aggregated data as outcome form aggregate_data.
#' @return A list of data for input int \code{rstan::sampling}.
#' @export
make_model_data <- function(agg_dat) {
  W <- model.matrix( ~ 0 + droplevels(current_eligible_vaccination), data = agg_dat)
  Sw <- diag(1, ncol(W)) - 1 / ncol(W)
  W_con <- constrain_design(W, Sw)

  des <- automaticr::get_design()
  # Rows of design matrix ordered by arm number,
  # so just repeat the appropriate number of times
  X <- des$X[dplyr::pull(agg_data, arm) + 1, ]

  return(list(
    N = length(agg_dat$y),
    K = ncol(X),
    L = ncol(W_con),
    M = nrow(des$Q),
    X = X,
    W = W_con,
    Xpred = des$X,
    Q = des$Q,
    y = agg_dat$y,
    n = agg_dat$trials,
    S = des$S,
    prior_only = 0)
  )
}


#' Calculate probability best
#'
#' @param draws MC model draws for arm means
#' @param sampsize Sample size for each arm
#' @return A \code{tibble} giving posterior quantities of interest for each arm.
#' @export
get_posterior_quantities <- function(draws, sampsize) {
  pars <- automaticr::get_trial_params()
  arm_post_summ <- tibble(
    arm = 0:(pars$n_arms - 1),
    sampsize = sampsize,
    mean = apply(draws, 2, mean),
    variance = diag(var(draws)),
    pbest = prob_best(draws),
    pbesttrt = c(NA, prob_best(draws[, -1])),
    pbeatctr = c(NA, col_comp(draws, 1, pars$delta)),
    inactive = c(FALSE, pbeatctr[-1] < pars$inactive_thres_ctr | pbest[-1] < pars$inactive_thres_sup),
    alloc_prob = automaticr::brar(pbest, sampsize, variance, inactive, pars$fix_ctrl_alloc)
  ) %>%
  bind_cols(as_tibble(t(HDInterval::hdi(draws))))
  # arm_post_summ <- dplyr::mutate(
  #   arm_post_summ,
  #   alloc_prob = brar(
  #     pbest,
  #     sampsize,
  #     variance),
  #   is_alloc = as.numeric(alloc_prob > 0))
  return(arm_post_summ)
}


#' Read interim log file.
#'
#' @param file Path to interim log.
#'
#' @return A \code{data.table} of interim log.
#' @export
read_interim_log <- function(file) {
  if(!file.exists(file)) {
    warning(paste(file, "not found. Exiting session."))
    quit(status = 0)
  } else {
    return(readr::read_csv(file))
  }
}


#' Update the interim log file
#'
#' The purpose of the interim log file is to track the current
#' interim, allocation probabilities, and allocation indicators,
#' as well as record when interims were undertaken and with
#' what sample size.
#'
#' @param file The path to the interim log file
#' @param date The date of the interim analysis
#' @param interim The interim sequence number
#' @param n_analysed The sample size analysed
#' @param prob_alloc The allocation probability vector
#' @param is_alloc Flag indicating arms which will receive allocations up to the next interim.
#'
#' @return NULL
#' @export
update_interim_log <- function(file, date, interim, n_analysed, prob_alloc, inactive) {
  if(!file.exists(file)) {
    warning(paste(file, "not found. Exiting session."))
    quit(status = 0)
  }

  pars <- automaticr::get_trial_params()
  tb <- tibble::tibble(
    interim_date = date,
    interim_num  = interim,
    n_analysed    = n_analysed,
    arm          = sprintf("%02d", 1:pars$n_arms),
    alloc_prob   = prob_alloc,
    inactive     = inactive
  )
  tb <- tidyr::spread(
    tidyr::unite(
      tidyr::gather(tb, variable, value, -(interim_date:arm)),
      temp, variable, arm),
    temp, value)
  readr::write_csv(tb, file, append = TRUE)
  return(invisible(NULL))
}


#' Generate an allocation sequence
#'
#' @param num_alloc Number of allocations to generate, i.e. length of sequence
#' @param alloc_prob Allocation probability to each arm
#' @param seed Seed for PRNG
#'
#' @return A numeric vector sequence of arm allocations
#' @export
generate_allocation_sequence <- function(num_alloc, alloc_prob, seed) {
  set.seed(seed)
  trial_par <- automaticr::get_trial_params()
  alloc_seq <- sample.int(trial_par$n_arms, num_alloc, prob = alloc_prob, replace = TRUE)
  rm(.Random.seed, envir = .GlobalEnv) # Decouple future RN.
  return(alloc_seq)
}


#' Write allocation sequence to file
#'
#' @param file Filepath for allocation sequence file
#' @param alloc_seq The allocation sequence
#' @param interim Integer giving interim from which allocation was generated
#'
#' @return NULL
#' @export
write_allocation_sequence <- function(file, alloc_seq, interim) {
  seq_dat <- tibble::tibble(
    TRTNO = sprintf("%02d_%04d", interim, 1:length(alloc_seq)),
    RANDOMISATION_OUTCOME = alloc_seq
  )
  readr::write_csv(seq_dat, file)
  message(paste("New allocation sequence written to", file))
  return(invisible(NULL))
}
