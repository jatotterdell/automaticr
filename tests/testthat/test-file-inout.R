context("File input and output")

tmp <- tempfile()
tmp_interim <- tempfile()
sim_dat <- generate_trial_data(1000, round(runif(0:9999)))

setup({
  write.csv(sim_dat, tmp, row.names = F)
})
teardown({
  unlink(tmp)
  unlink(tmp_interim)
})

test_that("read accumulating data", {
  dat <- read_raw_data(tmp)
  expect_equal(nrow(dat), nrow(sim_dat))

  proc_dat <- process_raw_data(dat)
  expect_equal(nrow(dat), nrow(sim_dat))
})

test_that("read/write interim log", {
  expect_error(read_interim_log(tmp_interim), "not found")

  init_interim_log(tmp_interim)
  expect_error(init_interim_log(tmp_interim), "already exists")

  interim_log <- read_interim_log(tmp_interim)
  expect_equal(nrow(interim_log), 1)
})