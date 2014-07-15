\name{predPref}
\alias{predPref}
\title{estimates parameters of predator preferences model and calculates LRT}
\usage{
predPref(eaten, caught, alpha = 0.05, em_maxiter = 100)
}
\arguments{
  \item{eaten}{a df of eatings preferences; TxS}

  \item{caught}{a df of caught prey species; TxS}

  \item{alpha}{LRT level of significance}

  \item{em_maxiter}{maximum number of iterations allowed
  for EM algorithm}
}
\description{
estimates parameters of predator preferences model and
calculates LRT
}

