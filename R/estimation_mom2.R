#' @export
is_legit_solution <- function(est)
{
  return(all((est >= 0) & (est <= 1)))
}


#' Helper function for method of moment for a 2-component Binomial Convolution model
#'
#' This function performs method of moment for a 2-component model
#' @param samples A vector of samples generated from a binomial convolution model
#' @param n_trials A vector specifying the number of trials
#' @returns A list of two solutions to the mom.
#' @export
mom_2_helper <- function(samples,
                         n_trials)
{
  mu = mean(samples)
  sig2 = var(samples)

  M1 = n_trials[1]
  M2 = n_trials[2]
  M = sum(n_trials)
  d = M1 * M2 * ((M - mu) * mu - M * sig2)
  if(d < 0)
  {
    return(list(est1=NA,
                est2=NA))
  }
  p11 = (M1 * mu - sqrt(d)) / (M1 * M)
  p12 = (M2 * mu + sqrt(d)) / (M2 * M)
  p21 = (M1 * mu + sqrt(d)) / (M1 * M)
  p22 = (M2 * mu - sqrt(d)) / (M2 * M)

  est1 = c(p11, p12)
  est2 = c(p21, p22)

  if(!is_legit_solution(est1))
  {
    est1 = NA
  }
  if(!is_legit_solution(est2))
  {
    est2 = NA
  }

  return(list(est1=est1,
              est2=est2))
}


#' Method of moment for a 2-component Binomial Convolution model
#'
#' This function performs method of moment for a 2-component model
#' @param samples A vector of samples generated from a binomial convolution model
#' @param n_trials A vector specifying the number of trials
#' @returns A vector of estimated probabilities
#' @export
mom_2 <- function(samples,
                  n_trials)
{
  est_candidates = mom_2_helper(samples=samples,
                                n_trials=n_trials)
  mu3 = mean((samples-mu)^3)

  ind = names(which(!is.na(est_candidates)))

  if(length(ind) == 0)
  {
    est = NA
  }else if(length(ind) == 1){
    est = est_candidates[[ind]]
  }else{
    est1_central_moments = central_moments(n_trials=n_trials,
                                           success_probs=est1)
    mu3_abs_diff_1 = abs(est1_central_moments[3] - mu3)
    est2_central_moments = central_moments(n_trials=n_trials,
                                           success_probs=est2)
    mu3_abs_diff_2 = abs(est2_central_moments[3] - mu3)

    if(mu3_abs_diff_1 < mu3_abs_diff_2)
    {
      est = est_candidates$est1
    }else{
      est = est_candidates$est2
    }
  }
  return(est)
}
