\name{hess}
\alias{hess}
\title{calculate second derivatives of Q, H}
\usage{
hess(lambda, gamma, c, Zdst, Ydst, J, I)
}
\arguments{
  \item{lambda}{parameter estimates of lambda from EM}

  \item{gamma}{parameter estimates of gamma from EM}

  \item{c}{parameter estimates of c from EM}

  \item{Zdst}{matrix of sums of binary response, whether
  prey species s was eaten or not during occurrence t; rows
  indexed by time, and cols indexed by prey species, TxS}

  \item{Ydst}{matrix sum of number of caught prey species s
  during occurrence t; rows indexed by time, and cols
  indexed by prey species, TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}
}
\description{
calculate second derivatives of Q, H
}

