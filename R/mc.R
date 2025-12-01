#' Monte Carlo integration with variance reduction
#'
#' Estimate \eqn{\theta = E[g(X)]} using Monte Carlo with optional variance
#' reduction: antithetic variables, control variates, importance sampling,
#' or stratified sampling.
#'
#' @param g Function \code{g(x)} to be averaged.
#' @param sampler Function to draw samples. For methods \code{"plain"},
#'   \code{"antithetic"} and \code{"control"}, this should generate
#'   \code{n} samples from the target distribution.
#' @param n Integer, number of samples.
#' @param method One of \code{"plain"}, \code{"antithetic"},
#'   \code{"control"}, \code{"importance"}, or \code{"stratified"}.
#' @param control List of method-specific options; see Details.
#'
#' @details
#' For \code{method = "control"}, supply \code{control$h} (control variate
#' function) and \code{control$mu_h} (known mean of \code{h(X)}).
#'
#' For \code{method = "importance"}, supply \code{control$sampler_q} to
#' draw samples from \eqn{q}, and \code{control$w} a function returning
#' the importance weight \eqn{f(x) / q(x)}.
#'
#' For \code{method = "stratified"}, \code{sampler} should accept an
#' argument \code{stratum} indicating which stratum to sample from, and
#' \code{control$strata} should give a vector of stratum labels and
#' \code{control$n_strata} a vector of sample sizes per stratum.
#'
#' @return A list with components \code{theta_hat} (estimate) and
#'   \code{se_hat} (Monte Carlo standard error).
#' @export
mc_integrate <- function(g, sampler, n,
                         method = c("plain", "antithetic",
                                    "control", "importance",
                                    "stratified"),
                         control = list()) {
  method <- match.arg(method)

  if (method == "plain") {
    x <- sampler(n)
    y <- g(x)
    return(.mc_summary(y))
  }

  if (method == "antithetic") {
    # Assume sampler generates U ~ Uniform(0,1)
    U <- sampler(n)
    U2 <- 1 - U
    y <- (g(U) + g(U2)) / 2
    return(.mc_summary(y))
  }

  if (method == "control") {
    h <- control$h
    mu_h <- control$mu_h
    if (is.null(h) || is.null(mu_h)) {
      stop("control$h and control$mu_h must be provided for control variates.")
    }
    x <- sampler(n)
    Y <- g(x)
    H <- h(x)
    beta <- stats::cov(Y, H) / stats::var(H)
    Y_cv <- Y - beta * (H - mu_h)
    return(.mc_summary(Y_cv))
  }

  if (method == "importance") {
    sampler_q <- control$sampler_q
    w_fun <- control$w
    if (is.null(sampler_q) || is.null(w_fun)) {
      stop("control$sampler_q and control$w must be provided for importance sampling.")
    }
    x <- sampler_q(n)
    w <- w_fun(x)
    y <- g(x) * w
    return(.mc_summary(y))
  }

  if (method == "stratified") {
    strata <- control$strata
    n_strata <- control$n_strata
    if (is.null(strata) || is.null(n_strata)) {
      stop("control$strata and control$n_strata must be provided for stratified sampling.")
    }
    if (length(strata) != length(n_strata)) stop("strata and n_strata lengths differ.")
    Y <- numeric(sum(n_strata))
    idx <- 1
    for (j in seq_along(strata)) {
      nj <- n_strata[j]
      xj <- sampler(nj, stratum = strata[j])
      yj <- g(xj)
      Y[idx:(idx + nj - 1L)] <- yj
      idx <- idx + nj
    }
    return(.mc_summary(Y))
  }

  stop("Unsupported method.")
}

.mc_summary <- function(y) {
  n <- length(y)
  theta_hat <- mean(y)
  se_hat <- stats::sd(y) / sqrt(n)
  list(theta_hat = theta_hat, se_hat = se_hat)
}
