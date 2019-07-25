context("File input and output")

tmp <- tempfile()
tmp_interim <- tempfile()
tmp_alloc <- tempfile()
sim_dat <- generate_trial_data(1000, sample(0:999, size = 1))

setup({
  write.csv(sim_dat, tmp, row.names = F)
})
teardown({
  unlink(tmp)
  unlink(tmp_interim)
  unlink(tmp_alloc)
})

test_that("read/process/aggregate accumulating data", {
  dat <- read_raw_data(tmp)
  expect_equal(nrow(dat), nrow(sim_dat))

  proc_dat <- process_raw_data(dat, Sys.Date())
  expect_equal(nrow(dat), nrow(sim_dat))

  agg_dat <- aggregate_data(dat = proc_dat)
})

test_that("read/write/update interim log", {
  expect_error(read_interim_log(file = tmp_interim), "not found")

  init_interim_log(file = tmp_interim)
  expect_error(init_interim_log(tmp_interim), "already exists")

  interim_log <- read_interim_log(tmp_interim)
  expect_equal(nrow(interim_log), 1)

  ap <- runif(13)
  ap <- ap / sum(ap)
  ap[ap < 0.1] <- 0
  ap <- ap / sum(ap)
  update_interim_log(tmp_interim, Sys.Date(), 1, 500, ap, as.numeric(ap > 0))
  interim_log2 <- read_interim_log(tmp_interim)
  expect_equal(nrow(interim_log2), 2)
})

test_that("read/write allocation sequence", {
  expect_error(generate_allocation_sequence(10, c(0.5,0.5), 1), "incorrect number")

  aseq <- generate_allocation_sequence(100, rep(1/13,13), 123)
  expect_equal(length(aseq), 100)

  write_allocation_sequence(tmp_alloc, aseq, 1)
  aseq_file <- readr::read_csv(tmp_alloc)
  expect_equal(ncol(aseq_file), 2)
  expect_equal(nrow(aseq_file), length(aseq))
})