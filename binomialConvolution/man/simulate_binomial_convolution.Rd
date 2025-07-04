% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_simulation.R
\name{simulate_binomial_convolution}
\alias{simulate_binomial_convolution}
\title{Simulate Binomial Convolution random variables}
\usage{
simulate_binomial_convolution(n_samples, n_trials, success_probs)
}
\arguments{
\item{n_samples}{Number of samples to be generated}

\item{n_trials}{A vector of positive integers each standing for the number of trials of one component}

\item{success_probs}{A vector of floats between 0 and 1 each standing for the probability of success of one component}
}
\value{
A vector of length n_samples
}
\description{
This function will generate a vector of Binomial Convolution
random variables given the parameters of each binomial component:
Mj and Pj
}
