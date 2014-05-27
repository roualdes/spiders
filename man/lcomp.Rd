\name{lcomp}
\alias{lcomp}
\title{complete data log-likelihood for null or simple alternative hypotheses}
\usage{
lcomp(Xdst, Ydst, lambda, gamma, J, I, c)
}
\arguments{
  \item{Xdst}{matrix of sums of number of eaten prey
  species s during occurrence t; rows indexed by time and
  cols indexed by prey species, TxS}

  \item{Ydst}{matrix sum of number of caught prey species s
  during occurrence t; rows indexed by time and cols
  indexed by prey species, TxS}

  \item{lambda}{matrix of parameters representing rates
  predator ate prey species s in time period t; TxS}

  \item{gamma}{matrix of parameters representing rates
  traps caught prey species s in time period t; TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}

  \item{c}{scalar or vector (indexed by time) in
  hypotheses}
}
\description{
complete data log-likelihood for null or simple alternative
hypotheses
}
\details{
this one should be used for optimization
}

