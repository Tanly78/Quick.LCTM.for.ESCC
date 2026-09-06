# tests/testthat/test-zscore.R
test_that("zscore_within_cohort produces zero-mean", {
  set.seed(1)
  data <- data.frame(
    id = rep(1:10, each = 2),
    timepoint = rep(c("T1", "T2"), 10),
    ALB = rnorm(20, 40, 5),
    TC = rnorm(20, 4, 1),
    LYM = rnorm(20, 1.5, 0.5)
  )
  out <- zscore_within_cohort(data, markers = c("ALB", "TC", "LYM"))
  expect_true("ALB_z" %in% names(out))
  expect_true("TC_z" %in% names(out))
  m <- tapply(out$ALB_z, out$timepoint, mean)
  expect_lt(abs(m["T1"]), 1e-8)
  expect_lt(abs(m["T2"]), 1e-8)
})

test_that("flip_crp negates CRP_z", {
  data <- data.frame(CRP_z = c(1, -2, 0.5, 3))
  out <- flip_crp(data)
  expect_equal(out$CRP_z_flip, c(-1, 2, -0.5, -3))
})