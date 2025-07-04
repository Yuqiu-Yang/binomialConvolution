% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{non_central_moments}
\alias{non_central_moments}
\title{The first six moments about 0 of a binomial convolution distribution}
\usage{
non_central_moments(n_trials, success_probs)
}
\arguments{
\item{n_trials}{A vector of positive integers each standing for the number of trials of one component}

\item{success_probs}{A vector of floats between 0 and 1 each standing for the probability of success of one component}
}
\value{
The first six moments about 0 of the specified binomial convolution distribution
}
\description{
This function returns the first six moments about 0 of a
binomial covolution distribution given the number of trials
and the success probabilities
}
