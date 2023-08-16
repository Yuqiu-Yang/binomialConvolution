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
