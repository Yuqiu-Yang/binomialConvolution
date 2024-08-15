#' Simulate Binomial Convolution random variables
#'
#' This function will generate a vector of Binomial Convolution
#' random variables given the parameters of each binomial component:
#' Mj and Pj
#'
#' @param n_samples Number of samples to be generated
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return A vector of length n_samples
#' @export
simulate_binomial_convolution <- function(n_samples,
                                          n_trials,
                                          success_probs)
{
  # Check input
  if((n_samples%%1!=0) | (n_samples<=0)) stop("n_samples should be a positive integer")
  if((any(n_trials%%1!=0)) | (any(n_trials<=0))) stop("n_trials should be a vector of positive integers")
  if((any(success_probs>1)) | (any(success_probs<0))) stop("success_probs should be a vector of probabilities")
  if(length(n_trials) != length(success_probs)) stop("n_trials and success_probs should be of the same length")

  # Final sample is just an accumulative sum of binomials
  n_components = length(n_trials)
  result = numeric(length=n_samples)
  for(component in 1 : n_components)
  {
    result = result + rbinom(n=n_samples,
                             size=n_trials[component],
                             prob=success_probs[component])
  }
  return(result)
}

