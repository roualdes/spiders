\name{se}
\alias{se}
\title{calculate standard error of estimates of predator preferences model}
\usage{
se(lambda, gamma, c, Xdst, Ydst, J, I)
}
\arguments{
  \item{lambda}{parameter estimates of lambda from EM}

  \item{gamma}{parameter estimates of gamma from EM}

  \item{c}{parameter estimates of c from EM}

  \item{Xdst}{matrix of sum of prey species s eaten during
  occurrence t; rows indexed by time, and cols indexed by
  prey species, TxS}

  \item{Ydst}{matrix of sum of number of caught prey
  species s during occurrence t; rows indexed by time, and
  cols indexed by prey species, TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}
}
\description{
calculate standard error of estimates of predator
preferences model
}
\details{
output is an array where the rownames are indexed by
{(T,S)}, T enumerated first
}

