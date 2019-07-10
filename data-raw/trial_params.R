trial_params <- list(
  n_arms = 13,
  sup_thres = 0.95,
  zero_alloc_thres = 0.01,
  fix_ctrl_alloc = 1/13,
  init_allocations = rep(1/13, 13)
)

interim_schedule <- data.frame(
  interim_num = 0:20,
  interim_n = seq(0, 1e4, 500),
  interim_new = c(NA, rep(500, 20)),
  alloc_seed = c(11871,
                 62871, 63579, 65231, 23841, 72859, 62859, 63511, 65213, 23897, 72801,
                 62873, 63577, 65233, 23843, 72857, 62857, 63513, 65217, 23899, 72803),
  stan_seed = c(11117,
                61871, 64579, 66231, 21841, 71859, 64511, 66211, 67213, 21897, 73801,
                61873, 64577, 66233, 21843, 71857, 64513, 66213, 67217, 21899, 73803)
)

design_data <- list(
  K = 2,
  N = 13,
  L1 = 3,
  L2 = 4,
  L3 = 12,
  X = cbind(c(1, rep(0, 12)), c(0, rep(1, 12))),
  Z1 = rbind(0, kronecker(rep(1,4), diag(1,3))),
  Z2 = rbind(0, kronecker(diag(1, 4), rep(1, 3))),
  Z3 = rbind(0, diag(1, 12))
)

usethis::use_data(
  trial_params,
  interim_schedule,
  design_data,
  internal = TRUE, overwrite = TRUE)