#' Erlang Density Function
#'
#' Provides the probability density function at point x for an erlang distribution of parameters k and lambda. X, k, and l can be vectors.
#'
#' @param x vector of quantiles.
#' @param k shape parameter, positive integer.
#' @param l short for lambda, rate parameter, must be greater than zero.
#' @param log logical: if TRUE, log of probability is returned.
#' @returns The probability of the Erlang probability distribution at x with parameters K and lambda.
#'
#' @examples
#' derlang(1, 1, 1)
#'
#' @export
derlang <- function(x, k, l, log = FALSE) {
  return(derlang_c(x, k, l, log))
}

#' Erlang Cumulative Density Function
#'
#' Provides the cumulative density function at point q for an erlang distribution of parameters k and lambda. X, k, and l can be vectors.
#'
#' @param q vector of quantiles.
#' @param k shape parameter, positive integer.
#' @param l short for lambda, rate parameter, must be greater than zero.
#' @param lower.tail logical: if TRUE, returns densities from 0 to q, otherwise q to 1.
#' @param log.p logical: if TRUE, log of probability is returned.
#' @returns The cumulative probability of the Erlang probability distribution at based on quantile q with parameters K and lambda.
#'
#' @examples
#' perlang(1, 1, 1)
#'
#' @export
perlang <- function(q, k, l, lower.tail = TRUE, log.p = FALSE) {
  return(perlang_c(q, k, l, lower.tail, log.p))
}



#' Gamma-Gompertz Distribution Function
#'
#' Provides probability density function for Gamma-Gompertz distribution.
#'
#' @param x vector of quantiles.
#' @param b scale paramater, must be greater than 0.
#' @param sigma,beta shape parameters, must be greater than 0.
#' @param log logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles x.
#'
#' @examples
#' dgamgomp(1,1,1,1)
#'
#'
#' @export
dgamgomp <- function(x, b, sigma, beta, log = FALSE) {
  return(dgamgomp_c(x, b, sigma, beta, log))
}


#' Gamma-Gompertz Cumulative Distribution Function
#'
#' Provides cumulative density function for Gamma-Gompertz distribution.
#'
#' @param q vector of quantiles.
#' @param b scale paramater, must be greater than 0.
#' @param sigma,beta shape parameters, must be greater than 0.
#' @param lower.tail logical: if TRUE, returns densities from 0 to q, otherwise q to 1.
#' @param log.p logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles q.
#'
#' @examples
#' pgamgomp(1,1,1,1)
#'
#'
#' @export
pgamgomp <- function(q, b, sigma, beta, lower.tail = TRUE, log.p = FALSE) {
  return(pgamgomp_c(q, b, sigma, beta, lower.tail, log.p))
}



#' Log Cauchy Distribution Functions
#'
#' Provides probability distribution function for Log Cauchy distribution.
#'
#' @param x vector of quantiles.
#' @param mu location parameter, must be real.
#' @param sigma scale parameter, must be greater than 0.
#' @param log logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles x.
#'
#' @examples
#' dlogcauchy(1,1,1)
#'
#' @export
dlogcauchy <- function(x, mu, sigma, log = FALSE) {
  return(dlogcauchy_c(x, mu, sigma, log))
}

#' Log Cauchy Cumulative Distribution Functions
#'
#' Provides cumulative distribution function for Log Cauchy distribution.
#'
#' @param q vector of quantiles.
#' @param mu location parameter, must be real.
#' @param sigma scale parameter, must be greater than 0.
#' @param lower.tail logical: if TRUE, returns densities from 0 to q, otherwise q to 1.
#' @param log.p logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles q.
#'
#' @examples
#' plogcauchy(1,1,1)
#'
#' @export
plogcauchy <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
  return(plogcauchy_c(q, mu, sigma, lower.tail, log.p))
}


#' Hypertabastic Distribution Function
#'
#' Provides probability distribution function for Hypertabastic distribution.
#'
#' @param x vector of quantiles.
#' @param a alpha parameter. Must be greater than 0.
#' @param b beta parameter. Must be greater than 0.
#' @param log logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles x.
#'
#' @examples
#' dhypertab(1,1,1)
#'
#' @export
dhypertab <- function(x, a, b, log = FALSE) {
  return(dhypertab_c(x, a, b, log))
}


#' Hypertabastic Cumulative Distribution Function
#'
#' Provides cumulative distribution function for Hypertabastic distribution.
#'
#' @param q vector of quantiles.
#' @param a alpha parameter. Must be greater than 0.
#' @param b beta parameter. Must be greater than 0.
#' @param lower.tail logical: if TRUE, returns densities from 0 to q, otherwise q to 1.
#' @param log.p logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles q.
#'
#' @examples
#' phypertab(1,1,1)
#'
#' @export
phypertab <- function(q, a, b, lower.tail = TRUE, log.p = FALSE) {
  return(phypertab_c(q, a, b, lower.tail, log.p))
}


#' Inverse Lindley Distribution Function
#'
#' Providers probability distribution function for Inverse Lindley distribution.
#'
#' @param x vector of quantiles.
#' @param theta paramater, must be greater than 0.
#' @param log logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles x.
#'
#' @examples
#' dinvlind(1,1)
#'
#' @export
dinvlind <- function(x, theta, log = FALSE) {
  return(dinvlind_c(x, theta, log))
}

#' Inverse Lindley Distribution Function
#'
#' Providers probability distribution function for Inverse Lindley distribution.
#'
#' @param q vector of quantiles.
#' @param theta paramater, must be greater than 0.
#' @param lower.tail logical: if TRUE, returns densities from 0 to q, otherwise q to 1.
#' @param log.p logical: if TRUE, log of probability is returned.
#' @returns Probabilities for quantiles q.
#'
#' @examples
#' pinvlind(1,1)
#'
#' @export
pinvlind <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
  return(pinvlind_c(q, theta, lower.tail, log.p))
}

