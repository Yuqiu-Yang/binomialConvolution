library(rootSolve)

#' Converts a vector of probabilities into parameters
#'
#' As a constraint optimization problem is in general more difficult
#' than an unconstraint one, we convert probabilities to numbers on the real line.
#' This function further considers the order. If there is no order imposed,
#' the probabilities are converted via the logit function.
#' On the other hand, the largest probability will be transformed via the logit.
#' The rest will be converted via the log of difference in logits.
#' @param probs A vector of probabilities
#' @param order_success_probs A vector specifying the order in descending order of probabilities
#' @return The transformed parameters of the vector of probabilities
#' @export
probs_to_par <- function(probs,
                         order_success_probs = NULL)
{
  if(is.null(order_success_probs))
  {
    par = log(probs/(1-probs))
  }else{
    n_trials = length(probs)
    probs = probs[order_success_probs]
    logits = log(probs/(1-probs))
    par = c(logits[1], log(-diff(logits)))
  }
  return(par)
}


#' Converts a vector of parameters into probabilities
#'
#' As a constraint optimization problem is in general more difficult
#' than an unconstraint one, we convert probabilities to numbers on the real line.
#' This function further considers the order. If there is no order imposed,
#' the probabilities are computed via the inverse logit function.
#' On the other hand, the first parameter will be transformed via the inverse logit.
#' The rest will be converted via the cumulative sum of exp of difference in logits and reordered.
#' @param par A vector of parameters
#' @param order_success_probs A vector specifying the order in descending order of probabilities
#' @return The transformed probabilities of the vector of parameters
#' @export
par_to_probs <- function(par,
                         order_success_probs = NULL)
{
  if(is.null(order_success_probs))
  {
    success_probs = 1/(1+exp(-par))
  }else{
    n_trials = length(par)
    if(n_trials > 1)
    {
      par[2 : n_trials] = -exp(par[2 : n_trials])
    }
    success_probs_ordered = cumsum(par)
    success_probs_ordered = 1/(1+exp(-success_probs_ordered))
    success_probs = success_probs_ordered[order(order_success_probs)]
  }
  return(success_probs)
}


#' MLE of a given sample from a Binomial Convolution distribution
#'
#' This function computes the MLE of a BC model given a sample and
#' the number of trials for each component. It also allows the users
#' to specify the order of probabilities as well as initial guesses.
#' @param samples A sample from the BC distribution
#' @param n_trials The number of trials for each component
#' @param computation_method The method to compute the PMF
#' @param return_mean If the results should be averaged over the sample
#' @param order_success_probs A vector specifying the order in descending order of probabilities
#' @param initial_probs The initial guess of the probabilies
#' @return The MLE and the corresponding log likelihood
#' @export
mle <- function(samples,
                n_trials,
                computation_method='convolution',
                return_mean=FALSE,
                order_success_probs = NULL,
                initial_probs = NULL)
{
  if(! is.null(initial_probs))
  {
    # If an inital guess is supplied
    # we convert it to parameters
    par = probs_to_par(probs=initial_probs,
                       order_success_probs=order_success_probs)
  }else{
    # If an initial guess is not supplied,
    # we use random guesses
    par = rnorm(length(n_trials))
  }

  # We run the optimization till it has converged
  has_converged = FALSE
  while(! has_converged)
  {
    results = optim(par=par,
                    fn=log_likelihood_order_constraint,
                    n_trials=n_trials,
                    samples=samples,
                    computation_method=computation_method,
                    return_mean=return_mean,
                    order_success_probs=order_success_probs,
                    control=list(fnscale=-1),
                    method = "BFGS")
    if(results$convergence != 0)
    {
      par = results$par
    }else{
      has_converged = TRUE
    }
  }

  # make sure the convergence message is 0
  # try number of iterations to 9999
  # number of function evaluations to 9999
  # Better initial values
  # Use nelder mead than other methods
  # Until convergence
  success_probs = par_to_probs(par=results$par,
                               order_success_probs=order_success_probs)
  return(list(probs=success_probs,
              ll=results$value))
}


#' Computes the log likelihood with possible order constraint
#'
#' This function computes the log likelihood with order constraint
#' @param par The parameters (not the probabilities)
#' @param n_trials The number of trials for each component
#' @param samples A sample from the BC distribution
#' @param computation_method The method to compute the PMF
#' @param return_mean If the results should be averaged over the sample
#' @param order_success_probs A vector specifying the order in descending order of probabilities
#' @return The log likelihood
#' @export
log_likelihood_order_constraint <- function(par,
                                            n_trials,
                                            samples,
                                            computation_method='convolution',
                                            return_mean=FALSE,
                                            order_success_probs = NULL)
{
  # The first entry of par is the largest probit of all the probs
  success_probs = par_to_probs(par=par,
                               order_success_probs=order_success_probs)
  return(log_likelihood(success_probs=success_probs,
                        n_trials=n_trials,
                        samples=samples,
                        computation_method=computation_method,
                        return_mean=return_mean))
}



#' Compute log likelihood
#'
#' This function computes the log likelihood given
#' a sample as well as the parameters
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param samples A vector of samples from a binomial convolution distribution
#' @param computation_method The method to use to compute the pmf
#' @param return_mean Whether or not to return average instead of the sum
#' @return the log likelihood
#' @export
log_likelihood <- function(success_probs,
                           n_trials,
                           samples,
                           computation_method='convolution',
                           return_mean=FALSE)
{
  pmf = binomial_convolution_pmf(n_trials=n_trials,
                                 success_probs=success_probs,
                                 computation_method=computation_method)
  if(return_mean)
  {
    l = mean(log(pmf(samples)))
  }else{
    l = sum(log(pmf(samples)))
  }
  return(l)
}


#' Compute the score function
#'
#' This function computes the score function given
#' a sample as well as the parameters
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param samples A vector of samples from a binomial convolution distribution
#' @param return_mean Whether or not to return average instead of the sum
#' @return the score function
#' @export
score <- function(success_probs,
                  n_trials,
                  samples,
                  return_mean=FALSE)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)
  N = length(samples)
  pmf = exact_pmf_dft(n_trials=n_trials,
                      success_probs=success_probs)
  omega = 2*pi / (sum_trials+1)

  x_list = dft_helper(n_trials=n_trials,
                      success_probs=success_probs)
  l = 0 : sum_trials
  gradient_matrix = matrix(0, nrow=n_components,
                           ncol=N)
  for(component in 1 : n_components)
  {
    temp_matrix = x_list$x
    temp_matrix[, component] = x_list$dx[, component]
    x = apply(temp_matrix, MARGIN=1, FUN=prod)
    for(n in 1 : N)
    {
      deno = Re(mean(exp(-1i*omega*(samples[n])*l) * x))
      gradient_matrix[component, n] = deno/pmf(samples[n])
    }
  }
  if(return_mean)
  {
    gradient_result = apply(gradient_matrix, MARGIN=1, FUN=mean)
  }else{
    gradient_result = apply(gradient_matrix, MARGIN=1, FUN=sum)
  }

  return(gradient_result)
}


#' Compute the hessian function
#'
#' This function computes the hessian function given
#' a sample as well as the parameters
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param samples A vector of samples from a binomial convolution distribution
#' @param return_mean Whether or not to return average instead of the sum
#' @return the hessian function
#' @export
hessian <- function(success_probs,
                    n_trials,
                    samples,
                    return_mean=FALSE)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)
  N = length(samples)
  pmf = exact_pmf_dft(n_trials=n_trials,
                      success_probs=success_probs)
  omega = 2*pi / (sum_trials+1)
  x_list = dft_helper(n_trials=n_trials,
                      success_probs=success_probs)
  l = 0 : sum_trials

  hessian_array = array(0, dim=c(n_components,
                                 n_components,
                                 N))
  for(component1 in 1 : n_components)
  {
    for(component2 in component1 : n_components)
    {
      if(component1 == component2)
      {
        temp_matrix1 = x_list$x
        temp_matrix1[, component1] = x_list$dx[, component1]
        x1 = apply(temp_matrix1, MARGIN=1, FUN=prod)
        temp_matrix2 = x_list$x
        temp_matrix2[, component1] = x_list$ddx[, component1]
        x2 = apply(temp_matrix2, MARGIN=1, FUN=prod)
        for(n in 1 : N)
        {
          deno1 = Re(mean(exp(-1i*omega*(samples[n])*l) * x1))^2

          deno2 = Re(mean(exp(-1i*omega*(samples[n])*l) * x2))

          hessian_array[component1, component2, n] = deno2/pmf(samples[n]) - deno1/(pmf(samples[n])^2)
        }
      }else{
        temp_matrix1 = x_list$x
        temp_matrix1[, component1] = x_list$dx[, component1]
        x1 = apply(temp_matrix1, MARGIN=1, FUN=prod)
        temp_matrix2 = x_list$x
        temp_matrix2[, component2] = x_list$dx[, component2]
        x2 = apply(temp_matrix2, MARGIN=1, FUN=prod)
        temp_matrix3 = x_list$x
        temp_matrix3[, component1] = x_list$dx[, component1]
        temp_matrix3[, component2] = x_list$dx[, component2]
        x3 = apply(temp_matrix3, MARGIN=1, FUN=prod)
        for(n in 1 : N)
        {
          deno1 = Re(mean(exp(-1i*omega*(samples[n])*l) * x1))

          deno2 = Re(mean(exp(-1i*omega*(samples[n])*l) * x2))

          deno3 = Re(mean(exp(-1i*omega*(samples[n])*l) * x3))

          hessian_array[component1, component2, n] = deno3/pmf(samples[n]) - deno1*deno2/(pmf(samples[n])^2)
          hessian_array[component2, component1, n] = hessian_array[component1, component2, n]
        }
      }
    }
  }
  hessian_result = rowSums(hessian_array, dims=2)
  if(return_mean)
  {
    hessian_result = hessian_result/N
  }
  return(hessian_result)
}





# library(plotly)
# #' @export
# objective_function_1 <- function(success_probs,
#                                  n_trials,
#                                  samples,
#                                  mode = "mgf",
#                                  weight_function=function(s){return(exp(-s^2))},
#                                  lower_limit=-5,
#                                  upper_limit=5)
# {
#   if(mode == "mgf")
#   {
#     egf = empirical_mgf(samples = samples)
#     gf = theoretical_mgf(n_trials = n_trials,
#                          success_probs = success_probs)
#   }else{
#     egf = empirical_cgf(samples = samples)
#     gf = theoretical_cgf(n_trials = n_trials,
#                           success_probs = success_probs)
#   }
#
#   integrand = discrepancy_integrand(theoretical_gf=gf,
#                                     empirical_gf=egf,
#                                     weight_function=weight_function)
#   discrepancy = compute_discrepancy(integrand=integrand,
#                                     lower_limit=lower_limit,
#                                     upper_limit=upper_limit)
#   return(discrepancy)
# }


# #' @export
# visualize_objective_function_1 <- function(success_probs,
#                                            n_trials,
#                                            samples,
#                                            weight_function=function(s){return(exp(-s^2))},
#                                            lower_limit=-5,
#                                            upper_limit=5,
#                                            grid_p1=seq(0, 1, length.out=50),
#                                            grid_p2=seq(0, 1, length.out=50))
# {
#   n_p1 = length(grid_p1)
#   n_p2 = length(grid_p2)
#   dis_m = matrix(0, nrow=n_p1, ncol=n_p2)
#   for(i in 1 : n_p1)
#   {
#     for(j in 1 : n_p2)
#     {
#       dis_m[i,j] = objective_function(ps=c(grid_p1[i], grid_p2[j]),
#                                       n_trials=n_trials,
#                                       samples=samples,
#                                       weight_function=weight_function,
#                                       lower_limit=lower_limit,
#                                       upper_limit=upper_limit)
#     }
#   }
#   plot_ly(x=grid_p1,
#           y=grid_p2,
#           z=dis_m,
#           contours=list(z=list(show=TRUE,
#                                start=0,
#                                end=1,
#                                size=0.02,
#                                color="red"))) %>% add_surface()
# }
#
