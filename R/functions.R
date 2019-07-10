#' Intialise the interim log
#'
#' Intialises the interim log based on \code{trial_params}.
#' Should only be called on first build.
#'
#' @param file The path for the interim log file
#'
#' @return NULL
#' @export
init_interim_log <- function(file) {
  if(file.exists(file)) {
    stop("File already exists.")
  }
  message(paste("Initialising", file))

  tb <- tibble::tibble(
    interim_date = Sys.Date(),
    interim_num  = 0,
    nrow_data    = 0,
    arm          = sprintf("%02d", 1:automaticr:::trial_params$n_arms),
    alloc_prob   = automaticr:::trial_params$init_allocations,
    is_alloc     = TRUE
  )
  tb <- tidyr::spread(
    tidyr::unite(
      tidyr::gather(tb, variable, value, -(interim_date:arm)),
      temp, variable, arm),
  temp, value)
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
    stop(paste(file, "file not found."))
  } else {
    dat <- readr::read_csv(
      file,
      na = c("null", "NA", ""))
  }
  if(nrow(dat) == 0)
    stop(paste(file, "has no records."))

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
#' @return A \code{data.table}.
#' @export
process_raw_data <- function(raw_dat, ref_date) {
  if(!lubridate::is.Date(ref_date))
    stop("ref_date is not a date.")
  if(nrow(raw_dat) == 0)
    stop("raw_dat has zero rows.")

  dat <- dplyr::mutate(raw_dat,
    date_vaccine_due = lubridate::as_date(date_vaccine_due),
    date_of_vaccination_administration = lubridate::as_date(date_of_vaccination_administration),
    time_to_vax = date_of_vaccination_administration - date_vaccine_due,
    time_since_due = ref_date - date_vaccine_due,
    vax_past_due = time_since_due > 28,
    on_time = time_to_vax <= 28,
    on_time = dplyr::if_else(vax_past_due, FALSE, NA)
  )
  return(dat)
}

#' Return index records for \code{parent_id}
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
  dplyr::top_n(dplyr::group_by(dat, parent_id), 1, date_vaccine_due)
}

#' Read interim log file.
#'
#' @param file Path to interim log.
#'
#' @return A \code{data.table} of interim log.
#' @export
read_interim_log <- function(file) {
  if(!file.exists(file)) {
    stop(paste(file, "not found."))
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
update_interim_log <- function(file, date, interim, n_analysed, prob_alloc, is_alloc) {
  if(!file.exists(file))
    stop(paste(file, "not found."))

  tb <- tibble::tibble(
    interim_date = date,
    interim_num  = interim,
    nrow_data    = n_analysed,
    arm          = sprintf("%02d", 1:automaticr:::trial_params$n_arms),
    alloc_prob   = prob_alloc,
    is_alloc     = is_alloc
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
  alloc_seq <- sample.int(automaticr:::trial_params$n_arms, num_alloc, prob = alloc_prob, replace = TRUE)
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
    seq_id = sprintf("%02d_%04d", interim, 1:length(alloc_seq)),
    randomisation_outcome = alloc_seq
  )
  readr::write_csv(seq_dat, file)
  message(paste("New allocation sequence written to", file))
  return(invisible(NULL))
}