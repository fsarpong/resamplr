test_that("bootstrap SE matches known SE for normal mean", {
  set.seed(1)
  n <- 50
  x <- rnorm(n, mean = 0, sd = 1)
  true_se <- 1 / sqrt(n)

  b <- boot_stat(x, mean, B = 500)
  se_hat <- boot_se(b)

  expect_true(abs(se_hat - true_se) < 0.1)
})
