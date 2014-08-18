\name{estCt}
\alias{estCt}
\title{estimates parameters from hypothesis lambda = gamma c'_t; S*T + T free parameters}
\usage{
estCt(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED)
}
\arguments{
  \item{Xdst}{matrix of sums of number of eaten prey
  species s during occurrence t; rows indexed by time, and
  cols indexed by prey species, TxS}

  \item{Ydst}{matrix sum of number of caught prey species s
  during occurrence t; rows indexed by time, and cols
  indexed by prey species, TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}

  \item{EM}{boolean; whether or not EM algorithm is used}

  \item{em_maxiter}{integer specifying max number of EM
  iterations}

  \item{BALANCED}{boolean; whether or not data are
  BALANCED}
}
\description{
estimates parameters from hypothesis lambda = gamma c'_t;
S*T + T free parameters
}

