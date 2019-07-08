#' Create simulated data for testing purposes
#'
#' Ideally this will cover all error/coding issues we
#' might expect to see so we can check our handling
#'
#' @param n Number of rows to simulate.
#' @param seed Seed for PRNG.
#'
#' @return A \code{data.frame} of simulated data.
#' @importFrom stats rexp
#' @export
generate_trial_data <- function(n, seed) {
  set.seed(seed)
  effect <- seq(1/35, 1/28, length.out = 13)
  l <- expand.grid(0:9, 0:9, LETTERS, stringsAsFactors = FALSE)
  L <- apply(l, 1, function(x) paste0(rev(x),  collapse = ''))
  randomisation_outcome <- sample(1:13, n, replace = TRUE)
  vax_due <- sample(seq(Sys.Date() - 90, Sys.Date(), by = "1 day"), n, replace = T)
  sms_date <-  vax_due + ifelse(randomisation_outcome %in% 1:5, 0,
                                ifelse(randomisation_outcome %in% 6:9, -7, 14))
  vax_date <- as.character(vax_due +
                             stats::rexp(n, effect[randomisation_outcome]))
  vax_date[vax_date >= Sys.Date()] <- "null"
  data.table::data.table(
    parent_id = sample(L, n, replace = TRUE),
    child_id = sample(L, n, replace = TRUE),
    clinic_id = sample(L, n, replace = TRUE),
    parent_postcode = sample(6000:6999, n, replace = TRUE),
    current_eligible_vaccination = sample(c("2m", "4m", "6m", "12m"), n, replace = TRUE),
    child_date_of_birth = "null",
    date_vaccine_due = as.character(vax_due),
    randomisation_outcome = randomisation_outcome,
    expected_date_sms_sent = as.character(sms_date),
    actual_date_sms_sent = as.character(sms_date),
    sms_delivery_failure_indicator = "success",
    date_of_vaccination_administration = vax_date,
    product_name_of_vaccine_administered="product x",
    opted_out = FALSE,
    key = c("parent_id", "date_vaccine_due")
  )
}
