test_that("permutation test type I error is close to alpha", {
  set.seed(1)
  nrep <- 200
  alpha <- 0.05
  rejections <- logical(nrep)
  for (i in seq_len(nrep)) {
    x <- rnorm(10)
    y <- rnorm(10)
    out <- perm_test(x, y, B = 999)
    rejections[i] <- (out$p_value < alpha)
  }
  type1 <- mean(rejections)
  expect_true(abs(type1 - alpha) < 0.03)
})
