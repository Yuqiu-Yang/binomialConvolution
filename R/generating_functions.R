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


#' Characteristic function of a Binomial Convolution distribution
#'
#' This function returns the CF of a Binomial Convolution distribution
#' given the number of trials and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The CF of the specified binomial convolution distribution
#' @export
theoretical_cf <- function(n_trials,
                           success_probs)
{
  n_components = length(n_trials)
  cf_list = vector(mode="list",
                   length=n_components)
  for(component in 1 : n_components)
  {
    cf_list[[component]] = local({n_trials = n_trials[component];
                                  success_prob=success_probs[component];
                                  binomial_cf(n_trials=n_trials,
                                              success_prob=success_prob)})
  }

  cf <- function(s)
  {
    result = 1
    for(component in 1 : n_components)
    {
      result = result * cf_list[[component]](s)
    }
    return(result)
  }
  return(Vectorize(cf, vectorize.args = "s"))
}


#' Moment generating function of a Binomial Convolution distribution
#'
#' This function returns the MGF of a Binomial Convolution distribution
#' given the number of trials and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The MGF of the specified binomial convolution distribution
#' @export
theoretical_mgf <- function(n_trials,
                            success_probs)
{
  n_components = length(n_trials)
  mgf_list = vector(mode="list",
                    length=n_components)
  for(component in 1 : n_components)
  {
    mgf_list[[component]] = local({n_trials = n_trials[component];
                                  success_prob=success_probs[component];
                                  binomial_mgf(n_trials=n_trials,
                                         success_prob=success_prob)})
  }

  mgf <- function(s)
  {
    result = 1
    for(component in 1 : n_components)
    {
      result = result * mgf_list[[component]](s)
    }
    return(result)
  }
  return(Vectorize(mgf, vectorize.args = "s"))
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


#' Cumulant generating function of a Binomial Convolution distribution
#'
#' This function returns the CGF of a Binomial Convolution distribution
#' given the number of trials and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The CGF of the specified binomial convolution distribution
#' @export
theoretical_cgf <- function(n_trials,
                            success_probs)
{
  n_components = length(n_trials)
  cgf_list = vector(mode="list",
                    length=n_components)
  for(component in 1 : n_components)
  {
    cgf_list[[component]] = local({n_trials = n_trials[component];
                                  success_prob=success_probs[component];
                                  binomial_cgf(n_trials=n_trials,
                                               success_prob=success_prob)})
  }

  cgf <- function(s)
  {
    result = 0
    for(component in 1 : n_components)
    {
      result = result + cgf_list[[component]](s)
    }
    return(result)
  }
  return(Vectorize(cgf, vectorize.args = "s"))
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


#' First order derivative of the cumulant generating function of a Binomial Convolution distribution
#'
#' This function returns the derivative of the CGF of a Binomial Convolution distribution
#' given the number of trials and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The derivative of the CGF of the specified binomial convolution distribution
#' @export
theoretical_dcgf <- function(n_trials,
                             success_probs)
{
  n_components = length(n_trials)
  dcgf_list = vector(mode="list",
                     length=n_components)
  for(component in 1 : n_components)
  {
    dcgf_list[[component]] = local({n_trials = n_trials[component];
                                    success_prob=success_probs[component];
                                    d_binomial_cgf(n_trials=n_trials,
                                                   success_prob=success_prob)})
  }

  dcgf <- function(s)
  {
    result = 0
    for(component in 1 : n_components)
    {
      result = result + dcgf_list[[component]](s)
    }
    return(result)
  }
  return(Vectorize(dcgf, vectorize.args = "s"))
}


#' Second order derivative of the cumulant generating function of a Binomial Convolution distribution
#'
#' This function returns the second order derivative of the CGF of a Binomial Convolution distribution
#' given the number of trials and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The second order derivative of the CGF of the specified binomial convolution distribution
#' @export
theoretical_ddcgf <- function(n_trials,
                              success_probs)
{
  n_components = length(n_trials)
  ddcgf_list = vector(mode="list",
                      length=n_components)
  for(component in 1 : n_components)
  {
    ddcgf_list[[component]] = local({n_trials = n_trials[component];
                                     success_prob=success_probs[component];
                                     dd_binomial_cgf(n_trials=n_trials,
                                                     success_prob=success_prob)})
  }

  ddcgf <- function(s)
  {
    result = 0
    for(component in 1 : n_components)
    {
      result = result + ddcgf_list[[component]](s)
    }
    return(result)
  }
  return(Vectorize(ddcgf, vectorize.args = "s"))
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


#' The first six cumulants of a binomial convolution distribution
#'
#' This function returns the first six cumulants of a
#' binomial covolution distribution given the number of trials
#' and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The first six cumulants of the specified binomial convolution distribution
#' @export
cumulants <- function(n_trials,
                      success_probs)
{
  n_components = length(n_trials)

  result = numeric(6)
  for(component in 1 : n_components)
  {
    result = result + binomial_cumulants(n_trials=n_trials[component],
                                         success_prob=success_probs[component])
  }
  return(result)
}


#' The first six moments about 0 of a binomial convolution distribution
#'
#' This function returns the first six moments about 0 of a
#' binomial covolution distribution given the number of trials
#' and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The first six moments about 0 of the specified binomial convolution distribution
#' @export
non_central_moments <- function(n_trials,
                                success_probs)
{
  k = cumulants(n_trials=n_trials,
                success_probs=success_probs)
  v1 = k[1]
  v2 = k[2] + (k[1])^2
  v3 = k[3] + 3 * k[2] * k[1] + (k[1])^3
  v4 = k[4] + 4 * k[3] * k[1] + 3 * (k[2])^2 + 6 * (k[2]) * (k[1])^2 + (k[1])^4
  v5 = k[5] + 5 * k[4] * k[1] + 10 * (k[3]) * (k[2]) + 10 * (k[3]) * (k[1])^2 + 15 * (k[2])^2 * (k[1]) + 10 * (k[2]) * (k[1])^3 + (k[1])^5
  v6 = k[6] + 6 * k[5] * k[1] + 15 * (k[4]) * (k[2]) + 15 * (k[4]) * (k[1])^2 + 10 * (k[3])^2 + 60 * (k[3]) * (k[2]) * (k[1]) + 20 * (k[3]) * (k[1])^3 + 15 * (k[2])^3 + 45 * (k[2])^2 * (k[1])^2 + 15 * (k[2]) * (k[1])^4 + (k[1])^6

  return(c(v1, v2, v3, v4, v5, v6))
}


#' The first six central moments of a binomial convolution distribution
#'
#' This function returns the first six central moments of a
#' binomial covolution distribution given the number of trials
#' and the success probabilities
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The first six central moments of the specified binomial convolution distribution
#' @export
central_moments <- function(n_trials,
                            success_probs)
{
  k = cumulants(n_trials=n_trials,
                success_probs=success_probs)
  v1 = 0
  v2 = k[2]
  v3 = k[3]
  v4 = k[4] + 3 * k[2]^2
  v5 = k[5] + 10 * (k[3]) * (k[2])
  v6 = k[6] + 15 * (k[4]) * (k[2]) + 10 * (k[3])^2 + 15 * (k[2])^3

  return(c(v1, v2, v3, v4, v5, v6))
}
