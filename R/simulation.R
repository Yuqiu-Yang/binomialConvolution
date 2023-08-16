#' Simulate Binomial Convolution random variables
#'
#' This function will generate a vector of Binomial Convolutio
#' random variables given the parameters of each binomail component:
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

  # Create a n_sample by n_component matrix
  # Each column contains IID samples from
  # one component binomial distribution
  n_components = length(n_trials)
  latent_matrix = matrix(NA,
                         nrow=n_samples,
                         ncol=n_components)
  for(component in 1 : n_components)
  {
    latent_matrix[, component] = rbinom(n=n_samples,
                                        size=n_trials[component],
                                        prob=success_probs[component])
  }
  # The final sample is simple the row sum of the matrix
  result = rowSums(latent_matrix)
  return(result)
}
