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


#' Converts a vector into a pmf
#'
#' This function takes a vector and turns it into
#' a pmf whose support starts at 0
#' @param pmf_vector A vector of probabilities
#' @return A pmf whose support starts at 0
make_pmf <- function(pmf_vector)
{
  pmf <- function(s)
  {
    if((s < 0) | (s > (length(pmf_vector)-1)))
    {
      return(0)
    }else{
      return(pmf_vector[s+1])
    }
  }
  return(Vectorize(pmf, vectorize.args = "s"))
}

#' Truncate a vector so that everything falls between the boundaries
#'
#' @param vec A vector
#' @param UB The upper bound
#' @param LB The lower bound
#' @return The trunccated vector
prob_fence <- function(vec, UB=1, LB=0){
  return(pmax(LB, pmin(vec, UB)))
}


#' Transform a R2 vector from Cartisean coordinate to polar coordinate
#'
#' @param vec A vector
#' @return The radius and the angle
vector_to_theta_r <- function(vec)
{
  r = sqrt(sum(vec^2))
  theta = atan(vec[2]/vec[1])
  if((all(vec < 0)) | ((vec[2] > 0) & (vec[1] < 0)))
  {
    theta = theta + pi
  }

  return(list(r=r,
              theta=theta))
}


#' Returns the index of the first element that is greater than or equal to p
#'
#' @param x A vector
#' @param p A vector
#' @return The radius and the angle
v_first_which <- Vectorize(FUN=function(x, p) {
  ind = which(x >= p)
  if(length(ind) > 0)
  {
    result = ind[1]
  }else{
    result = length(x)
  }
  return(result)
}, vectorize.args = "p")


#' Returns the index of the last element that is greater than p
#'
#' @param x A vector
#' @param p A vector
#' @return The radius and the angle
v_last_which <- Vectorize(FUN=function(x, p) {
  ind = which(x > p)
  if(length(ind) > 0)
  {
    result = ind[length(ind)]
  }else{
    result = 0
  }
  return(result)
}, vectorize.args = "p")


#' Check if a vector contains all probabilities
#'
#' @param est A vector
#' @return If all elements are probabilities
is_legit_solution <- function(est)
{
  return(all((est >= 0) & (est <= 1)))
}
