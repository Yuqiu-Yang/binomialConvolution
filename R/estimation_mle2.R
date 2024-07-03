#' Computes MLE for a 2-component BC model
#'
#' Due to the pecularity of a BC distribution, we single out the 2-component
#' model as the simplies case.
#' @param samples A sample from the BC distribution
#' @param n_trials The number of trials for each component
#' @param computation_method The method to compute the PMF
#' @param return_mean If the results should be averaged over the sample
#' @param initial_probs The initial guess of the probabilities
#' @param report_both If report both solutions
#' @return The MLE estimates as well as their corresponding log likelihood
#' @export
mle2 <- function(samples,
                 n_trials,
                 computation_method='convolution',
                 return_mean=FALSE,
                 initial_probs = NULL,
                 report_both=FALSE)
{
  mle_est1 = mle(samples = samples,
                 n_trials = n_trials,
                 computation_method = computation_method,
                 return_mean=return_mean,
                 order_success_probs = c(1,2),
                 initial_probs = NULL)

  mle_est2 = mle(samples = samples,
                 n_trials = n_trials,
                 computation_method = computation_method,
                 return_mean=return_mean,
                 order_success_probs = c(2,1),
                 initial_probs = NULL)
  if(report_both)
  {
    return(list(est1=mle_est1,
                est2=mle_est2))
  }else{
    l1 = mle_est1$ll
    l2 = mle_est2$ll
    if(l1 > l2)
    {
      return(mle_est1)
    }else{
      return(mle_est2)
    }
  }
}


#' Given a radius and an angle, compute the difference between the chisq stats
#'
#' This function is used for computing the joint confidence region based on
#' the likelihood ratio test.
#' @param r The radius
#' @param theta The angle
#' @param solution_to_use The MLE solution
#' @param solution_nllh The negative log likelihood at the MLE
#' @param n_trials The number of trials for each component
#' @param samples A sample from the BC distribution
#' @param sig_level The significance level
#' @return The difference between chisq stats
#' @export
diff_chisq <- Vectorize(function(r, theta,
                         solution_to_use,
                         solution_nllh,
                         n_trials,
                         samples,
                         sig_level=0.95)
  {
    current_point = prob_fence(solution_to_use + r * c(cos(theta), sin(theta)))

    current_nllh = -log_likelihood(current_point,n_trials,samples)

    v = 2 * (current_nllh - solution_nllh)

    return(v - qchisq(sig_level, 2))
  },
  vectorize.args = "r")


#' Given two points, sample points in the middle to trace out the contour of the confidence region
#'
#' @param theta1 The angle of the first point
#' @param theta2 The angle of the second point
#' @param r1 The radius of the first point
#' @param r2 The radius of the second point
#' @param solution_to_use The MLE solution
#' @param solution_nllh The negative log likelihood at the MLE
#' @param n_trials The number of trials for each component
#' @param samples A sample from the BC distribution
#' @return The contour traced out by points between the given two points
#' @export
trace_contour <- function(theta1, theta2,
                          r1, r2,
                          solution_to_use,
                          solution_nllh,
                          n_trials,
                          samples,
                          pots=c())
{
  if(abs(theta1-theta2)<(pi/64))
  {
    return(pots)
  }
  theta = (theta1 + theta2)/2
  tryCatch(
    {
      radius_info = rootSolve::uniroot.all(diff_chisq,
                                theta=theta,
                                solution_to_use=solution_to_use,
                                solution_nllh=solution_nllh,
                                n_trials=n_trials,
                                samples=samples,
                                sig_level=0.95,
                                lower=0,
                                upper=1)
    },
    error = function(cond){
      radius_info = list(root = 100)
    }
  )
  if(length(radius_info) < 1)
  {
    radius_info = list(root = 100)
  }else{
    if(class(radius_info) != "list")
    {
      radius_info = list(root = min(radius_info))
    }
  }
  current_r = radius_info$root
  pot = prob_fence((radius_info$root) * c(cos(theta), sin(theta)) + solution_to_use)
  pots = rbind(pots, pot)
  pots=trace_contour(theta1=theta, theta2=theta2,
                     r1=current_r,r2=r2,
                     solution_to_use=solution_to_use,
                     solution_nllh=solution_nllh,
                     n_trials=n_trials,
                     samples=samples,
                     pots=pots)
  pots=trace_contour(theta1=theta1, theta2=theta,
                     r1=r1, r2=current_r,
                     solution_to_use=solution_to_use,
                     solution_nllh=solution_nllh,
                     n_trials=n_trials,
                     samples=samples,
                     pots=pots)
  return(pots)
}


#' Use points to trace out the contour of the confidence region
#'
#' @param solution_to_use The MLE solution
#' @param solution_nllh The negative log likelihood at the MLE
#' @param n_trials The number of trials for each component
#' @param samples A sample from the BC distribution
#' @return The contour traced out by points
#' @export
trace_confidence_region <- function(solution_to_use,
                                    solution_nllh,
                                    n_trials,
                                    samples)
{
  # We first approximate the region using an ellipse
  hess = binomialConvolution::hessian(success_probs = solution_to_use,
                                      n_trials=n_trials,
                                      samples=samples)
  eigv = eigen(hess)$vectors
  # Then we find the four points along the major and minor axies
  theta_r_pair = matrix(0,
                        ncol=4,
                        nrow=4)
  theta_r_pair[,1] = c(vector_to_theta_r(eigv[,1])$theta,
                       vector_to_theta_r(eigv[,1])$theta+pi,
                       vector_to_theta_r(eigv[,2])$theta,
                       vector_to_theta_r(eigv[,2])$theta+pi)
  for(i in 1 : nrow(theta_r_pair))
  {
    theta = theta_r_pair[i, 1]
    tryCatch(
      {
        radius_info = rootSolve::uniroot.all(diff_chisq,
                                  theta=theta,
                                  solution_to_use=solution_to_use,
                                  solution_nllh=solution_nllh,
                                  n_trials=n_trials,
                                  samples=samples,
                                  sig_level=0.95,
                                  lower=0,
                                  upper=1)
      },
      error = function(cond){
        radius_info = list(root = 100)
      }
    )
    if(length(radius_info) < 1)
    {
      radius_info = list(root = 100)
    }else{
      radius_info = list(root = min(radius_info))
    }
    theta_r_pair[i,2] = radius_info$root
    theta_r_pair[i, 3:4] = prob_fence((theta_r_pair[i,2]) * c(cos(theta), sin(theta)) + solution_to_use)
  }

  pots = trace_contour(theta_r_pair[1,1], theta_r_pair[2,1],
                       theta_r_pair[1,2],theta_r_pair[2,2],
                       solution_to_use, solution_nllh, n_trials, samples, c())
  pots = trace_contour(theta_r_pair[1,1]+pi, theta_r_pair[2,1]+pi,
                       theta_r_pair[1,2],theta_r_pair[2,2],
                       solution_to_use, solution_nllh, n_trials, samples, pots)
  return(pots)
}


#' Compute the chisq stat based on profile likelihood
#'
#' @param n_trials The number of trials for each component
#' @param samples A sample from the BC distribution
#' @param fixed_success_probs Probabilities of the fixed components
#' @param fixed_success_probs_index The corresponding indecies of fixed_success_probs
#' @returns The Chisq statistics
#' @export
profile_likelihood_chisq <- function(n_trials,
                                     samples,
                                     fixed_success_probs,
                                     fixed_success_probs_index)
{
  n_grid <- 200
  fixed_success_prob_grid <- seq(0,1,length.out=n_grid)
  profile_ll <- numeric(n_grid)
  for (j in 1:n_grid) {
    profile_ll[j] <- -profile_log_likelihood(samples=samples,
                                           n_trials=n_trials,
                                           fixed_success_probs=fixed_success_prob_grid[j],
                                           fixed_success_probs_index=fixed_success_probs_index)$ll
  }
  return(list("prob_grid"=fixed_success_prob_grid,
              "chisq"=2*(profile_ll-min(profile_ll))))
}


