#' Create simulated data for testing purposes
#'
#' Ideally this will cover all error/coding issues we
#' might expect to see so we can check our handling
#'
#' @param n Number of rows to simulate.
#' @param seed Seed for PRNG.
#' @param effect Rate used to generate exponential RV as offset on vax_due, see details.
#'
#' @details
#' \itemize{
#' \item Currently we compute the date on which the vaccination took place as
#' the due-date plus some offset governed by a draw from an exponential
#' distribution parameterised with rates from passed in by the calling
#' function. Defaults to linear declining mean (1/rate) from 35 to 28 days.
#' If the computed vaccination date occurs after the current date, it is
#' set to 'null'.
#' \item All other stochastic elements (vaccination due date, randomisation arm,
#' current vaccination type etc) are obtained via simple random sample with
#' replacement. Vax due-date is sampled from 90 days prior to the current date
#' up to the current date.
#' }
#'
#' @return A \code{data.frame} of simulated data.
#' @importFrom stats rexp
#' @export
generate_trial_data <- function(n, seed, effect = seq(1/35, 1/28, length.out = 13)) {
  set.seed(seed)
  l <- expand.grid(0:9, 0:9, 0:9, LETTERS, stringsAsFactors = FALSE)
  L <- apply(l, 1, function(x) paste0(rev(x),  collapse = ''))
  randomisation_outcome <- sample(1:13, n, replace = TRUE)
  vax_due <- sample(seq(Sys.Date() - 90, Sys.Date(), by = "1 day"), n, replace = T)
  sms_date <-  vax_due + ifelse(randomisation_outcome %in% 1:5, 0,
                                ifelse(randomisation_outcome %in% 6:9, -7, 14))
  vax_date <- as.character(vax_due +
                             stats::rexp(n, effect[randomisation_outcome]))
  # is 'null' what we anticipate from redcap/middleware?
  vax_date[vax_date >= Sys.Date()] <- "null"
  tibble::tibble(
    parent_id = sample(L, n, replace = TRUE),
    # so child and clinic id can be same as parent id? or these are just
    # dummy values that we do not rely on?
    child_id = sample(L, n, replace = TRUE),
    clinic_id = sample(L, n, replace = TRUE),
    parent_postcode = sample(6000:6999, n, replace = TRUE),
    current_eligible_vaccination = sample(c("2m", "4m", "6m", "12m"), n, replace = TRUE),
    child_date_of_birth = ifelse(runif(n) < 0.1, "null", as.character(lubridate::as_date(Sys.Date() - rnorm(n, 90, 10)))),
    date_vaccine_due = as.character(vax_due),
    randomisation_outcome = randomisation_outcome,
    expected_date_sms_sent = as.character(sms_date),
    actual_date_sms_sent = as.character(sms_date),
    sms_delivery_failure_indicator = "success",
    date_of_vaccination_administration = vax_date,
    product_name_of_vaccine_administered="product x",
    opted_out = FALSE
  )
}

#' Create messy (with errors) simulated data for testing purposes
#'
#' Ideally this will cover all error/coding issues we
#' might expect to see so we can check our handling
#'
#' @param n Number of rows to simulate.
#' @param seed Seed for PRNG.
#' @param effect Rate used to generate exponential RV as offset on vax_due, see details.
#'
#' @details
#' \itemize{
#' \item Currently we compute the date on which the vaccination took place as
#' the due-date plus some offset governed by a draw from an exponential
#' distribution parameterised with rates from passed in by the calling
#' function. Defaults to linear declining mean (1/rate) from 35 to 28 days.
#' If the computed vaccination date occurs after the current date, it is
#' set to 'null'.
#' \item All other stochastic elements (vaccination due date, randomisation arm,
#' current vaccination type etc) are obtained via simple random sample with
#' replacement. Vax due-date is sampled from 90 days prior to the current date
#' up to the current date.
#' }
#'
#' @return A \code{data.frame} of simulated data.
#' @importFrom stats rexp
#' @export
generate_messy_trial_data <- function(n, seed, effect = seq(1/35, 1/28, length.out = 13)) {
  set.seed(seed)
  l <- expand.grid(0:9, 0:9, 0:9, LETTERS, stringsAsFactors = FALSE)
  L <- apply(l, 1, function(x) paste0(rev(x),  collapse = ''))
  randomisation_outcome <- sample(1:13, n, replace = TRUE)
  vax_due <- sample(seq(Sys.Date() - 90, Sys.Date() + 90, by = "1 day"), n, replace = T)
  sms_date <-  vax_due + ifelse(randomisation_outcome %in% 1:5, 0,
                                ifelse(randomisation_outcome %in% 6:9, -7, 14))
  vax_date <- as.character(vax_due +
                             sample(c(-1,1), n, replace = TRUE)*stats::rexp(n, effect[randomisation_outcome]))
  vax_date <- ifelse(runif(n) < 0.05, stringr::str_remove_all(vax_date, "-"), vax_date)
  # is 'null' what we anticipate from redcap/middleware?
  vax_date[vax_date >= Sys.Date()] <- "null"
  tibble::tibble(
    parent_id = sample(L, n, replace = TRUE),
    # so child and clinic id can be same as parent id? or these are just
    # dummy values that we do not rely on?
    child_id = sample(L, n, replace = TRUE),
    clinic_id = sample(L, n, replace = TRUE),
    parent_postcode = sample(6000:6999, n, replace = TRUE),
    current_eligible_vaccination = sample(c("2m", "4m", "6m", "12m"), n, replace = TRUE),
    child_date_of_birth = ifelse(runif(n) < 0.1, "null", as.character(lubridate::as_date(Sys.Date() - rnorm(n, 90, 10)))),
    date_vaccine_due = as.character(vax_due),
    randomisation_outcome = randomisation_outcome,
    expected_date_sms_sent = as.character(sms_date),
    actual_date_sms_sent = as.character(sms_date),
    sms_delivery_failure_indicator = "success",
    date_of_vaccination_administration = vax_date,
    product_name_of_vaccine_administered="product x",
    opted_out = FALSE
  )
}

