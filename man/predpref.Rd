\name{predPref}
\alias{predPref}
\title{estimates parameters of predator preferences model and calculates LRT}
\usage{
predPref(eaten, caught, alpha = 0.05, index_c = TRUE, em_maxiter = 100)
}
\arguments{
  \item{eaten}{a df of eatings preferences; TxS}

  \item{caught}{a df of caught prey species; TxS}

  \item{alpha}{LRT level of significance}

  \item{index_c}{boolean; TRUE indexes c in null hypothesis
  by t}

  \item{em_maxiter}{maximum number of iterations allowed
  for EM algorithm}
}
\description{
estimates parameters of predator preferences model and
calculates LRT
}

