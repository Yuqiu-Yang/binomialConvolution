#' Moment generating function of a binomial distribution
#'
#' This function returns the MGF of a binomial distribution
#' given the number of trials and probability of success
#'
#' @param n_trials The number of trials
#' @param success_prob The probability of success
#' @return The MGF of the specified binomial distribution
#' @export
binomial_mgf <- function(n_trials,
                         success_prob)
{
  mgf <- function(s)
  {
    return((1-success_prob + success_prob * exp(s))^n_trials)
  }
  return(Vectorize(mgf, vectorize.args = "s"))
}


#' Cumulant generating function of a binomial distribution
#'
#' This function returns the CGF of a binomial distribution
#' given the number of trials and probability of success
#'
#' @param n_trials The number of trials
#' @param success_prob The probability of success
#' @return The CGF of the specified binomial distribution
#' @export
binomial_cgf <- function(n_trials,
                         success_prob)
{
  cgf <- function(s)
  {
    return(n_trials * log(1-success_prob + success_prob * exp(s)))
  }
  return(Vectorize(cgf, vectorize.args = "s"))
}


#' Characteristic function of a binomial distribution
#'
#' This function returns the CF of a binomial distribution
#' given the number of trials and probability of success
#'
#' @param n_trials The number of trials
#' @param success_prob The probability of success
#' @return The CF of the specified binomial distribution
#' @export
binomial_cf <- function(n_trials,
                        success_prob)
{
  cf <- function(s)
  {
    return((1-success_prob + success_prob * exp(1i*s))^n_trials)
  }
  return(Vectorize(cf, vectorize.args = "s"))
}


#' First order derivative of the cumulant generating function of a binomial distribution
#'
#' This function returns the derivative of the CGF of a binomial distribution
#' given the number of trials and probability of success
#'
#' @param n_trials The number of trials
#' @param success_prob The probability of success
#' @return The derivative of the CGF of the specified binomial distribution
#' @export
d_binomial_cgf <- function(n_trials,
                           success_prob)
{
  dcgf <- function(s)
  {
    if(s >= 0)
    {
      return(n_trials * success_prob / ((1-success_prob)*exp(-s) + success_prob))
    }else{
      return(n_trials * success_prob * exp(s) / ((1-success_prob) + success_prob * exp(s)))
    }
  }
  return(Vectorize(dcgf, vectorize.args = "s"))
}


#' Second order derivative of the cumulant generating function of a binomial distribution
#'
#' This function returns the second order derivative of the CGF of a binomial distribution
#' given the number of trials and probability of success
#'
#' @param n_trials The number of trials
#' @param success_prob The probability of success
#' @return The second order derivative of the CGF of the specified binomial distribution
#' @export
dd_binomial_cgf <- function(n_trials,
                            success_prob)
{
  ddcgf <- function(s)
  {
    if(s >=0)
    {
      temp = success_prob / ((1-success_prob)*exp(-s) + success_prob)
    }else{
      temp = success_prob * exp(s) / ((1-success_prob) + success_prob * exp(s))
    }
    return(n_trials * temp * (1-temp))
  }
  return(Vectorize(ddcgf, vectorize.args = "s"))
}


#' Empirical estimates of the MGF of a Binomial Convolution distribution
#'
#' This function returns an empirical estimator of the MGF of
#' a binomial convolution distribution given a set of samples
#' @param samples A vector of IID samples
#' @return The EMGF of a binomial convolution distribution
#' @export
empirical_mgf <- function(samples)
{
  mgf <- function(s)
  {
    return(mean(exp(s * samples)))
  }
  return(Vectorize(mgf, vectorize.args = "s"))
}


#' Empirical estimates of the CF of a Binomial Convolution distribution
#'
#' This function returns an empirical estimator of the CF of
#' a binomial convolution distribution given a set of samples
#' @param samples A vector of IID samples
#' @return The ECF of a binomial convolution distribution
#' @export
empirical_cf <- function(samples)
{
  cf <- function(s)
  {
    return(mean(exp(1i*s * samples)))
  }
  return(Vectorize(cf, vectorize.args = "s"))
}


#' Empirical estimates of the CGF of a Binomial Convolution distribution
#'
#' This function returns an empirical estimator of the CGF of
#' a binomial convolution distribution given a set of samples
#' @param samples A vector of IID samples
#' @return The ECGF of a binomial convolution distribution
#' @export
empirical_cgf <- function(samples)
{
  emgf = empirical_mgf(samples=samples)
  cgf <- function(s)
  {
    return(log(emgf(s)))
  }
  return(Vectorize(cgf, vectorize.args = "s"))
}


#' The first six cumulants of a binomial distribution
#'
#' This function returns the first six cumulants of a
#' binomial distribution given the number of trials
#' and the success probability
#' @param n_trials A positive integer standing for the number of trials
#' @param success_prob A floats between 0 and 1 standing for the probability of success
#' @return The first six cumulants of the specified binomial distribution
#' @export
binomial_cumulants <- function(n_trials,
                               success_prob)
{
  k1 = n_trials * success_prob
  k2 = n_trials * success_prob * (1 - success_prob)
  k3 = n_trials * success_prob * (1 - success_prob) * (1 - 2 * success_prob)
  k4 = n_trials * success_prob * (1 - success_prob) * (1 - 6 * success_prob * (1 - success_prob))
  k5 = n_trials * success_prob * (1 - success_prob) * (1 - 2 * success_prob) * (1 - 12 * success_prob + 12 * success_prob^2)
  k6 = n_trials * success_prob * (1 - success_prob) * (1 - 30 * success_prob + 150 * success_prob^2 - 240 * success_prob^3 + 120 * success_prob^4)
  return(c(k1, k2, k3, k4, k5, k6))
}

