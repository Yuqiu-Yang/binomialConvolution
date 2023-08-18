library(plotly)



#' @export
objective_function_1 <- function(success_probs,
                                 n_trials,
                                 samples,
                                 mode = "mgf",
                                 weight_function=function(s){return(exp(-s^2))},
                                 lower_limit=-5,
                                 upper_limit=5)
{
  if(mode == "mgf")
  {
    egf = empirical_mgf(samples = samples)
    gf = theoretical_mgf(n_trials = n_trials,
                         success_probs = success_probs)
  }else{
    egf = empirical_cgf(samples = samples)
    gf = theoretical_cgf(n_trials = n_trials,
                          success_probs = success_probs)
  }

  integrand = discrepancy_integrand(theoretical_gf=gf,
                                    empirical_gf=egf,
                                    weight_function=weight_function)
  discrepancy = compute_discrepancy(integrand=integrand,
                                    lower_limit=lower_limit,
                                    upper_limit=upper_limit)
  return(discrepancy)
}


#' @export
visualize_objective_function_1 <- function(success_probs,
                                           n_trials,
                                           samples,
                                           weight_function=function(s){return(exp(-s^2))},
                                           lower_limit=-5,
                                           upper_limit=5,
                                           grid_p1=seq(0, 1, length.out=50),
                                           grid_p2=seq(0, 1, length.out=50))
{
  n_p1 = length(grid_p1)
  n_p2 = length(grid_p2)
  dis_m = matrix(0, nrow=n_p1, ncol=n_p2)
  for(i in 1 : n_p1)
  {
    for(j in 1 : n_p2)
    {
      dis_m[i,j] = objective_function(ps=c(grid_p1[i], grid_p2[j]),
                                      n_trials=n_trials,
                                      samples=samples,
                                      weight_function=weight_function,
                                      lower_limit=lower_limit,
                                      upper_limit=upper_limit)
    }
  }
  plot_ly(x=grid_p1,
          y=grid_p2,
          z=dis_m,
          contours=list(z=list(show=TRUE,
                               start=0,
                               end=1,
                               size=0.02,
                               color="red"))) %>% add_surface()
}

