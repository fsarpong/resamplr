#' Bootstrap statistic
#'
#' Generate bootstrap replicates of a statistic \eqn{T = t(X_1,\dots,X_n)}.
#' Ordinary bootstrap draws i.i.d. samples from the empirical distribution
#' \eqn{\hat F_n}; the balanced bootstrap ensures each index appears exactly
#' \code{B} times across all resamples.
#'
#' @param x A numeric vector (or object) containing the data.
#' @param stat A function \code{stat(x, ...)} returning a numeric scalar.
#' @param B Integer, number of bootstrap resamples.
#' @param type Character, one of \code{"ordinary"} or \code{"balanced"}.
#' @param ... Additional arguments passed to \code{stat}.
#'
#' @return A list with elements
#' \describe{
#'   \item{t0}{Observed statistic on the original data.}
#'   \item{tstar}{Numeric vector of length \code{B} of bootstrap replicates.}
#' }
#' @examples
#' x <- rnorm(50)
#' b <- boot_stat(x, stat = mean, B = 200)
#' boot_se(b)
#' boot_ci(b, method = "percentile")
#' @export
boot_stat <- function(x, stat, B = 999L,
                      type = c("ordinary", "balanced"), ...) {
  type <- match.arg(type)
  n <- length(x)
  if (n <= 1L) stop("Need at least 2 observations.")
  if (B <= 1L) stop("B must be > 1.")
  t0 <- stat(x, ...)
  tstar <- numeric(B)

  if (type == "ordinary") {
    for (b in seq_len(B)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      tstar[b] <- stat(x[idx], ...)
    }
  } else { # balanced
    # Create an n x B index matrix with each index used B times
    idx <- matrix(sample(rep(seq_len(n), times = B)),
                  nrow = n, ncol = B)
    for (b in seq_len(B)) {
      tstar[b] <- stat(x[idx[, b]], ...)
    }
  }

  structure(list(t0 = t0, tstar = tstar), class = "boot_resamplr")
}

#' Bootstrap standard error
#'
#' @param boot An object returned by \code{\link{boot_stat}}.
#' @return Numeric scalar, bootstrap estimate of the standard error of \eqn{T}.
#' @export
boot_se <- function(boot) {
  tstar <- boot$tstar
  stats::sd(tstar)
}

#' Bootstrap bias estimate
#'
#' @param boot An object returned by \code{\link{boot_stat}}.
#' @return Numeric scalar, bootstrap estimate of the bias \eqn{E^*[T^*]-T}.
#' @export
boot_bias <- function(boot) {
  mean(boot$tstar) - boot$t0
}

#' Bootstrap confidence interval
#'
#' Compute bootstrap confidence intervals: normal, basic, percentile,
#' or t-based (studentized) intervals.
#'
#' @param boot An object from \code{\link{boot_stat}}.
#' @param method Interval type: \code{"normal"}, \code{"basic"},
#'   \code{"percentile"}, or \code{"t"}.
#' @param level Confidence level between 0 and 1.
#' @param tstar_se Optional vector of bootstrap standard errors for
#'   studentized intervals (if \code{method = "t"}).
#'
#' @return A numeric vector of length 2: lower and upper endpoints.
#' @export
boot_ci <- function(boot,
                    method = c("normal", "basic", "percentile", "t"),
                    level = 0.95,
                    tstar_se = NULL) {
  method <- match.arg(method)
  alpha <- 1 - level
  t0 <- boot$t0
  tstar <- sort(boot$tstar)
  B <- length(tstar)

  if (method == "percentile") {
    k_lo <- max(1, floor((alpha / 2) * B))
    k_hi <- min(B, ceiling((1 - alpha / 2) * B))
    return(c(tstar[k_lo], tstar[k_hi]))
  }

  mean_t <- mean(tstar)
  se <- stats::sd(tstar)
  z <- stats::qnorm(1 - alpha / 2)

  if (method == "normal") {
    return(c(t0 - z * se, t0 + z * se))
  }

  if (method == "basic") {
    # basic: (2 T0 - T_(1-alpha/2), 2 T0 - T_(alpha/2))
    k_lo <- max(1, floor((alpha / 2) * B))
    k_hi <- min(B, ceiling((1 - alpha / 2) * B))
    return(c(2 * t0 - tstar[k_hi], 2 * t0 - tstar[k_lo]))
  }

  # t-based (studentized) interval: needs tstar_se
  if (is.null(tstar_se)) {
    stop("t-based interval requires tstar_se values.")
  }
  if (length(tstar_se) != B) stop("tstar_se must have length B.")
  t_stats <- (tstar - t0) / tstar_se
  t_stats <- sort(t_stats)
  k_lo <- max(1, floor((1 - alpha / 2) * B))
  k_hi <- min(B, ceiling((alpha / 2) * B))
  t_lo <- t_stats[k_lo]
  t_hi <- t_stats[k_hi]
  c(t0 - t_lo * se, t0 - t_hi * se)
}
