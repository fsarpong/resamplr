test_that("kfold_cv returns finite estimate", {
  set.seed(1)
  n <- 40
  x <- rnorm(n)
  y <- 2 * x + rnorm(n)
  dat <- data.frame(x = x)

  model_fit <- function(data, y, ...) lm(y ~ x, data = data)
  predict_fun <- function(fit, data) predict(fit, newdata = data)
  loss <- function(y_true, y_pred) mean((y_true - y_pred)^2)

  out <- kfold_cv(dat, y, model_fit, predict_fun, loss, K = 5)
  expect_true(is.finite(out$cv_estimate))
})
