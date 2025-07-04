% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{pbinomconv}
\alias{pbinomconv}
\title{CDF of a Binomial Convolution}
\usage{
pbinomconv(q, n_trials, success_probs, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
\item{q}{Vector of quantiles}

\item{n_trials}{A vector of positive integers each standing for the number of trials of one component}

\item{success_probs}{A vector of floats between 0 and 1 each standing for the probability of success of one component}

\item{n_samples}{Number of samples to be generated}
}
\value{
A vector of cumulative probabilities
}
\description{
CDF of a Binomial Convolution
}
