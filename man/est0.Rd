\name{est0}
\alias{est0}
\title{estimates parameters from the null model when data unbalanced}
\usage{
est0(Xdst, Ydst, J, I)
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
}
\description{
estimates parameters from the null model when data
unbalanced
}

