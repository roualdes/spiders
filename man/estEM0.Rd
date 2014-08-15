\name{estEM0}
\alias{estEM0}
\title{estimates parameters from the null model with non-count data via EM}
\usage{
estEM0(Zdst, Ydst, J, I, indexC, em_maxiter)
}
\arguments{
  \item{Zdst}{matrix of sums of binary response, whether
  prey species s was eaten or not during occurrence t; rows
  indexed by time, and cols indexed by prey species, TxS}

  \item{Ydst}{matrix sum of number of caught prey species s
  during occurrence t; rows indexed by time, and cols
  indexed by prey species, TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}

  \item{indexC}{boolean; TRUE indexes c in null hypothesis
  by t}

  \item{em_maxiter}{maximum number of iterations allowed
  for EM algorithm}
}
\description{
estimates parameters from the null model with non-count
data via EM
}

