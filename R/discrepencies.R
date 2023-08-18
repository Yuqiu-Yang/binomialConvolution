#' The integrand of the final discrepancy function
#'
#' This function returns the integrand of the discrepancy function
#' given the theoretical generating function, its empirical
#' estimator, and a weight function
#' @param theoretical_gf The theoretical generating function (MGF OR CGF)
#' @param empirical_gf The empirical estimator of the gf (MGF OR CGF)
#' @param weight_function A weighting function that ensures the convergence of the integral
#' @return The integrand of the discrepancy function
#' @export
discrepancy_integrand <- function(theoretical_gf,
                                  empirical_gf,
                                  weight_function)
{
  integrand <- function(s)
  {
    return((theoretical_gf(s) - empirical_gf(s))^2 * weight_function(s))
  }
  return(integrand)
}


#' Compute discrepancy
#'
#' This function performs numerical integration given
#' the integrand as well as the limits of the integral
#' @param integrand The integrand function
#' @param lower_limit The lower bound of the integral
#' @param upper_limit The upper bound of the integral
#' @return The discrepancy value
#' @export
compute_discrepancy <- function(integrand,
                                lower_limit,
                                upper_limit)
{
  discrepancy = integrate(integrand,
                          lower = lower_limit,
                          upper = upper_limit)$value
  return(discrepancy)
}


