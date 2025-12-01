#' Two-sample permutation test
#'
#' Perform a permutation test for difference between two samples using
#' an arbitrary test statistic.
#'
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param stat Function \code{stat(x, y)} returning numeric scalar.
#'   Default is difference in means, \code{mean(x) - mean(y)}.
#' @param B Number of random permutations.
#' @param alternative One of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param exact Logical; if \code{TRUE} and \code{length(x)+length(y)} is small,
#'   enumerate all permutations to obtain an exact p-value.
#'
#' @return A list with components
#' \describe{
#'   \item{stat_obs}{Observed test statistic.}
#'   \item{perm_stats}{Vector of permutation statistics.}
#'   \item{p_value}{Estimated p-value (ASL).}
#' }
#' @export
perm_test <- function(x, y,
                      stat = function(x, y) mean(x) - mean(y),
                      B = 9999L,
                      alternative = c("two.sided", "less", "greater"),
                      exact = FALSE) {
  alternative <- match.arg(alternative)
  x <- as.numeric(x); y <- as.numeric(y)
  z <- c(x, y)
  n1 <- length(x); n <- length(z)

  stat_obs <- stat(x, y)

  if (exact && n <= 10) {
    # Enumerate all combinations of positions for group 1
    combs <- utils::combn(n, n1)
    B <- ncol(combs)
    perm_stats <- numeric(B)
    for (b in seq_len(B)) {
      idx <- combs[, b]
      xb <- z[idx]
      yb <- z[-idx]
      perm_stats[b] <- stat(xb, yb)
    }
  } else {
    perm_stats <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n)
      xb <- z[idx[1:n1]]
      yb <- z[idx[(n1 + 1):n]]
      perm_stats[b] <- stat(xb, yb)
    }
  }

  # Achieved significance level (ASL)
  if (alternative == "two.sided") {
    p <- mean(abs(perm_stats) >= abs(stat_obs))
  } else if (alternative == "greater") {
    p <- mean(perm_stats >= stat_obs)
  } else {
    p <- mean(perm_stats <= stat_obs)
  }

  list(stat_obs = stat_obs,
       perm_stats = perm_stats,
       p_value = p)
}
