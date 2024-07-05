#' Binomial Convolution random variables
#'
#' @param n_samples Number of samples to be generated
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return A vector of length n_samples
#' @export
rbinomconv <- function(n,
                       n_trials,
                       success_probs)
{
  result = simulate_binomial_convolution(n_samples=n,
                                         n_trials=n_trials,
                                         success_probs=success_probs)
  return(result)
}


#' Density of a Binomial Convolution
#'
#' @param x Vector of quantiles
#' @param n_samples Number of samples to be generated
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return A vector of probabilities
#' @export
dbinomconv <- function(x,
                       n_trials,
                       success_probs,
                       log=FALSE)
{
  pmf = binomial_convolution_pmf(n_trials=n_trials,
                                 success_probs=success_probs)
  result = pmf(x)
  if(log)
  {
    result = log(result)
  }
  return(result)
}


#' CDF of a Binomial Convolution
#'
#' @param q Vector of quantiles
#' @param n_samples Number of samples to be generated
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return A vector of cumulative probabilities
#' @export
pbinomconv <- function(q,
                       n_trials,
                       success_probs,
                       lower.tail=TRUE,
                       log.p=FALSE)
{
  pmf = binomial_convolution_pmf(n_trials=n_trials,
                                 success_probs=success_probs)
  s = cumsum(pmf(0:sum(n_trials)))
  ind_neg = which(q < 0)
  int_q = pmax(0, pmin(floor(q), sum(n_trials))) + 1
  result = s[int_q]
  result[ind_neg] = 0
  if(!lower.tail)
  {
    result = 1 - result
  }
  if(log.p)
  {
    result = log(result)
  }
  return(result)
}


#' Inverse CDF of a Binomial Convolution
#'
#' @param p Vector of probabilities
#' @param n_samples Number of samples to be generated
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return A vector of quantiles
#' @export
qbinomconv <- function(p,
                       n_trials,
                       success_probs,
                       lower.tail=TRUE,
                       log.p=FALSE)
{
  pmf = binomial_convolution_pmf(n_trials=n_trials,
                                 success_probs=success_probs)
  if(log.p)
  {
    p = exp(p)
  }
  p = prob_fence(p-.Machine$double.eps)
  s = cumsum(pmf(0:sum(n_trials)))
  if(lower.tail)
  {
    result = v_first_which(x=s, p=p)-1
  }else{
    s=1-s
    s[length(s)] = 0
    result = v_last_which(x=s, p=p)
  }

  return(result)
}


#' The pmf of a binomial convolution distribution
#'
#' This function computes the pmf of a specified
#' binomial convolution distribution
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @param computation_method The method to use to compute the pmf
#' @return The pmf of the specified binomial convolution distribution
#' @export
binomial_convolution_pmf <- function(n_trials,
                                     success_probs,
                                     computation_method="convolution")
{
  if(computation_method == "convolution")
  {
    pmf = exact_pmf_convolution(n_trials=n_trials,
                                success_probs=success_probs)
  }else if(computation_method == "recursive"){
    pmf = exact_pmf_recursive(n_trials=n_trials,
                              success_probs=success_probs)
  }else if(computation_method == "saddlepoint"){
    pmf = approximate_pmf_saddlepoint(n_trials=n_trials,
                                      success_probs=success_probs)
  }else if(computation_method == "kolmogorov"){
    pmf = approximate_pmf_kolmogorov(n_trials=n_trials,
                                     success_probs=success_probs)
  }else if(computation_method == "dft"){
    pmf = exact_pmf_dft(n_trials=n_trials,
                        success_probs=success_probs)
  }else{
    stop("Unknown method. Has to be one of 'convolution', 'recursive', 'saddlepoint', 'kolmogorov'")
  }
  return(pmf)
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

