#' K-fold cross-validation
#'
#' Estimate predictive risk using K-fold cross-validation for a generic model.
#'
#' @param data Data frame or list of inputs.
#' @param y Response vector (numeric or factor).
#' @param model_fit Function \code{model_fit(data_train, y_train, ...)}
#'   returning a fitted model object.
#' @param predict_fun Function \code{predict_fun(fit, data_test)} returning
#'   predicted values.
#' @param loss Function \code{loss(y_true, y_pred)} returning mean loss.
#' @param K Number of folds.
#' @param stratify Logical; for classification, preserve class proportions.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments passed to \code{model_fit}.
#'
#' @return List with components \code{cv_estimate} and \code{fold_losses}.
#' @export
kfold_cv <- function(data, y,
                     model_fit,
                     predict_fun,
                     loss,
                     K = 5L,
                     stratify = FALSE,
                     seed = NULL,
                     ...) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  if (K < 2L || K > n) stop("K must be between 2 and n.")

  if (stratify && is.factor(y)) {
    # stratified folds
    folds <- rep(NA_integer_, n)
    for (lvl in levels(y)) {
      idx <- which(y == lvl)
      idx <- sample(idx)
      k_seq <- rep(seq_len(K), length.out = length(idx))
      folds[idx] <- k_seq
    }
  } else {
    folds <- sample(rep(seq_len(K), length.out = n))
  }

  fold_losses <- numeric(K)

  for (k in seq_len(K)) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)

    data_train <- data[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    data_test <- data[test_idx, , drop = FALSE]
    y_test <- y[test_idx]

    fit <- model_fit(data_train, y_train, ...)
    y_pred <- predict_fun(fit, data_test)
    fold_losses[k] <- loss(y_test, y_pred)
  }

  list(cv_estimate = mean(fold_losses),
       fold_losses = fold_losses,
       folds = folds)
}

#' Leave-one-out cross-validation
#'
#' @inheritParams kfold_cv
#' @export
loo_cv <- function(data, y,
                   model_fit,
                   predict_fun,
                   loss,
                   ...) {
  kfold_cv(data, y,
           model_fit = model_fit,
           predict_fun = predict_fun,
           loss = loss,
           K = length(y),
           stratify = FALSE,
           seed = NULL,
           ...)
}
