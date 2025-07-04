#' Computes some key quantities in the DFT
#'
#' This function computes some key quantities used repeatedly
#' in the DFT computation.
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return A list of quantities used in DFT computation
#' @export
dft_helper <- function(n_trials,
                       success_probs)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)
  omega = 2*pi/(1+sum_trials)
  x_matrix = dx_matrix = ddx_matrix = matrix(0+0*1i,
                                            nrow=1+sum_trials,
                                            ncol=n_components)
  for(l in 0 : sum_trials)
  {
    for(j in 1 : n_components)
    {
      pj = success_probs[j]
      Mj = n_trials[j]
      x_matrix[l+1, j] = (1-pj+pj*exp(1i*omega*l))^(Mj)
      dx_matrix[l+1, j] = exp(1i*omega*l)-1
      if(Mj > 1)
      {
        dx_matrix[l+1, j] = dx_matrix[l+1, j] * Mj * ((1-pj+pj*exp(1i*omega*l))^(Mj-1))
      }
      ddx_matrix[l+1, j] = (exp(1i*omega*l)-1)^2
      if(Mj < 2)
      {
        ddx_matrix[l+1, j] = 0 + 0*1i
      }else if(Mj == 2){
        ddx_matrix[l+1, j] = ddx_matrix[l+1, j] * (Mj * (Mj-1)) * ((1-pj+pj*exp(1i*omega*l))^(Mj-2))
      }else{
        ddx_matrix[l+1, j] = ddx_matrix[l+1, j] * 2
      }
    }
  }
  return(list(x=x_matrix,
              dx=dx_matrix,
              ddx=ddx_matrix))
}


#' Use convolution to compute the exact pmf of a binomial convolution distribution
#'
#' This function uses convolution to compute the pmf of a specified
#' binomial convolution distribution
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The pmf of the specified binomial convolution distribution
#' @export
exact_pmf_convolution <- function(n_trials,
                                  success_probs)
{
  # Instead of enumerating all possible combinations
  # such that the sum of a specified number
  # we borrow ideas from dynamic programming
  # where we store intermediate results in tables to be
  # looked up later on
  sum_trials = sum(n_trials)
  n_components = length(n_trials)
  # First, we store the mass of each binomial component
  # in a table
  # The number of row is the sum of all numbers of trials +1
  # representing the maximum that we can get
  # This would of course waste some space as the first
  # several components would not reach the maximum surely but it's simple
  individual_table = matrix(0, nrow=sum_trials+1,
                            ncol=n_components)
  for(component in 1 : n_components)
  {
    individual_table[,component] = dbinom(x=seq(0, sum_trials),
                                          size=n_trials[component],
                                          prob=success_probs[component])
    # individual_table[(n_trials[component]+2):(sum_trials+1), component]=0
  }

  # convolution_table now stores the mass functions
  convolution_table = matrix(0, nrow=sum_trials+1,
                             ncol=n_components)
  # We first fill in the first column which is identical
  # to a binomial distribution
  convolution_table[,1] = dbinom(x=seq(0, sum_trials),
                                 size=n_trials[1],
                                 prob=success_probs[1])
  if(n_components > 1)
  {
    # Then we add components one by one
    for(component in 2 : n_components)
    {
      for(s in 0 : sum_trials)
      {
        p = 0
        for(i in 0 : s)
        {
          p = p + convolution_table[i+1, component-1] * individual_table[s-i+1, component]
        }
        convolution_table[s+1, component] = p
      }
    }
  }
  # The last column is the pmf desired
  pmf_vector = convolution_table[, n_components]
  pmf = make_pmf(pmf_vector=pmf_vector)
  return(pmf)
}


#' Use recursion to compute the exact pmf of a binomial convolution distribution
#'
#' This function uses recursion to compute the pmf of a specified
#' binomial convolution distribution
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The pmf of the specified binomial convolution distribution
#' @export
exact_pmf_recursive <- function(n_trials,
                                success_probs)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)

  pmf_vector = numeric(sum_trials+1)
  # s=0
  # For the sum to be 0, every component is 0
  p = 1
  for(component in 1 : n_components)
  {
    p = p * (1 - success_probs[component])^(n_trials[component])
  }
  pmf_vector[1] = p
  # s > 0
  for(s in 1 : sum_trials)
  {
    p = 0
    for(j in 1 : s)
    {
      temp = 0
      for(component in 1 : n_components)
      {
        temp = temp + n_trials[component] * (success_probs[component]/(1-success_probs[component]))^j
      }
      temp = temp * pmf_vector[s-j+1] * (-1)^(j-1)
      p = p + temp
    }
    pmf_vector[s+1] = p/s
  }
  pmf = make_pmf(pmf_vector=pmf_vector)
  return(pmf)
}


#' Use saddlepoint to approximate the pmf of a binomial convolution distribution
#'
#' This function uses saddlepoint approximation to compute the pmf of a specified
#' binomial convolution distribution
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The pmf of the specified binomial convolution distribution
#' @export
approximate_pmf_saddlepoint <- function(n_trials,
                                        success_probs)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)

  pmf_vector = numeric(sum_trials+1)
  # At the boundaries, the saddlepoint equaitons can not
  # be solved. However, they are easy to compute
  p = 1
  for(component in 1 : n_components)
  {
    p = p * (1 - success_probs[component])^(n_trials[component])
  }
  pmf_vector[1] = p

  p = 1
  for(component in 1 : n_components)
  {
    p = p * (success_probs[component])^(n_trials[component])
  }
  pmf_vector[sum_trials+1] = p

  # For the rest of the values
  for(s in 1 : (sum_trials-1))
  {
    # First solve the saddlepoint equation
    cgf = theoretical_cgf(n_trials=n_trials,
                          success_probs=success_probs)
    dcgf = theoretical_dcgf(n_trials=n_trials,
                            success_probs=success_probs)
    mu_hat = uniroot(function(x){return(dcgf(x)-s)},
                     interval = c(-100,100))$root

    ddcgf = theoretical_ddcgf(n_trials=n_trials,
                              success_probs=success_probs)
    pmf_vector[s+1] = exp(cgf(mu_hat) - mu_hat * s)/(sqrt(2*pi * ddcgf(mu_hat)))

  }
  # As saddlepoint approximation does not guarantee a
  # pmf, we normalize the result keeping the boundaries as they are
  temp = pmf_vector[1] + pmf_vector[sum_trials+1]
  pmf_vector[2:sum_trials] = pmf_vector[2:sum_trials]/sum(pmf_vector[2:sum_trials]) * (1-temp)
  pmf = make_pmf(pmf_vector=pmf_vector)
  return(pmf)
}


#' Use Kolmogorov to approximate the pmf of a binomial convolution distribution
#'
#' This function uses Kolmogorov approximation to compute the pmf of a specified
#' binomial convolution distribution
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The pmf of the specified binomial convolution distribution
#' @export
approximate_pmf_kolmogorov <- function(n_trials,
                                       success_probs)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)

  p0 = dbinom(x=seq(0,  sum_trials),
              size=sum_trials,
              prob=sum(n_trials*success_probs)/sum_trials)

  delta_p0 = matrix(0, nrow=6,
                    ncol=length(p0))
  delta_p0[1,] = diff(c(0, p0))
  for(k in 2 : 6)
  {
    delta_p0[k,] = diff(c(0, delta_p0[k-1,]))
  }

  mu_1_0 = sum(n_trials*success_probs)

  vs = moments(n_trials=n_trials,
               success_probs=success_probs)
  e_vs = moments(n_trials=sum_trials,
                 success_probs=0.1)
  a_js = numeric(6)
  a_js[1] = (-1) * (vs[1] - mu_1_0)

  p_s = matrix(0, nrow=6,
               ncol=length(p0))
  p_s[1,] = p0 + a_js[1] * delta_p0[1,]

  for(k in 2 : 6)
  {
    temp_0 = 0
    for(j in 1 : (k-1))
    {
      temp = 0
      for(i in j : k)
      {
        c_kji = choose(k, i) * (-1)^j * factorial(j) * kStatistics::nStirling2(i,j)
        if(i == k)
        {
          temp = temp + c_kji
        }else{
          temp = temp + c_kji * e_vs[k-i]
        }
      }
      temp = temp * a_js[j]
      temp_0 = temp_0 + temp
    }
    mu_k_k_1 = e_vs[k] + temp_0

    a_js[k] = (-1)^k * (vs[k] - mu_k_k_1)/factorial(k)
    p_s[k,] = p_s[k-1,] + a_js[k] * delta_p0[k,]
  }
  pmf_vector = p_s[6,]/sum(p_s[6,])
  pmf = make_pmf(pmf_vector=pmf_vector)
  return(pmf)
}


#' Use discrete fourier transform to compute the exact pmf of a binomial convolution distribution
#'
#' This function uses dft to compute the pmf of a specified
#' binomial convolution distribution
#' @param n_trials A vector of positive integers each standing for the number of trials of one component
#' @param success_probs A vector of floats between 0 and 1 each standing for the probability of success of one component
#' @return The pmf of the specified binomial convolution distribution
#' @export
exact_pmf_dft <- function(n_trials,
                          success_probs)
{
  sum_trials = sum(n_trials)
  n_components = length(n_trials)

  omega = 2*pi / (sum_trials+1)

  x_list = dft_helper(n_trials=n_trials,
                      success_probs=success_probs)
  temp_matrix = x_list$x
  x = apply(temp_matrix, MARGIN=1, FUN=prod)
  pmf_vector = numeric(1+sum_trials)
  l = 0 : sum_trials
  for(s in 0 : sum_trials)
  {
    pmf_vector[s+1] = Re(mean(x*exp(-1i * omega * s * l)))
  }

  # pmf_vector = pmf_vector/sum(pmf_vector)
  pmf = make_pmf(pmf_vector=pmf_vector)
  return(pmf)
}



