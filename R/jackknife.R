#' Jackknife statistic
#'
#' Compute leave-one-out replicates \eqn{T_{(i)}} and jackknife pseudo-values
#' for a statistic \eqn{T = t(X_1,\dots,X_n)}.
#'
#' @param x Numeric vector.
#' @param stat Function \code{stat(x, ...)} returning numeric scalar.
#' @param ... Additional arguments to \code{stat}.
#'
#' @return A list with components:
#' \describe{
#'   \item{t0}{Statistic on full data.}
#'   \item{t_i}{Vector of leave-one-out statistics \eqn{T_{(i)}}.}
#'   \item{pseudo}{Pseudo-values \eqn{p_i = n T - (n-1) T_{(i)}}.}
#' }
#' @export
jack_stat <- function(x, stat, ...) {
  n <- length(x)
  if (n <= 1L) stop("Need at least 2 observations.")
  t0 <- stat(x, ...)
  t_i <- numeric(n)
  for (i in seq_len(n)) {
    t_i[i] <- stat(x[-i], ...)
  }
  pseudo <- n * t0 - (n - 1) * t_i
  structure(list(t0 = t0, t_i = t_i, pseudo = pseudo),
            class = "jack_resamplr")
}

#' Jackknife standard error
#'
#' @param jack Object from \code{\link{jack_stat}}.
#' @return Jackknife estimate of standard error.
#' @export
jack_se <- function(jack) {
  t_i <- jack$t_i
  n <- length(t_i)
  mean_t <- mean(t_i)
  sqrt((n - 1) / n * sum((t_i - mean_t)^2))
}

#' Jackknife bias estimate
#'
#' @param jack Object from \code{\link{jack_stat}}.
#' @return Jackknife estimate of bias.
#' @export
jack_bias <- function(jack) {
  n <- length(jack$t_i)
  mean_p <- mean(jack$pseudo)
  mean_p - jack$t0
}
