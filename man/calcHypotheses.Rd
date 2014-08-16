\name{calcHypotheses}
\alias{calcHypotheses}
\title{function to calculate hypotheses given a user specifed null and alternative}
\usage{
calcHypotheses(hyp, Xdst, Ydst, J, I, balanced, EM, em_maxiter)
}
\arguments{
  \item{hyp}{a 2-tuple specifying the null and alternative
  hypotheses, respectively}

  \item{Xdst}{matrix of sums of number of eaten prey
  species s during occurrence t; rows indexed by time, and
  cols indexed by prey species, TxS}

  \item{Ydst}{matrix sum of number of caught prey species s
  during occurrence t; rows indexed by time, and cols
  indexed by prey species, TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}

  \item{balanced}{boolean specifying balanced data or not}

  \item{EM}{boolean specifying if EM algorithm should be
  used}

  \item{em_maxiter}{maximum number of iterations allowed
  for EM algorithm}
}
\description{
function to calculate hypotheses given a user specifed null
and alternative
}

