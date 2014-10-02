\name{llEM}
\alias{llEM}
\title{non-observed count log-likelihood of predators preferances}
\usage{
llEM(Zdst, Ydst, lambda, gamma, J, I, c = NULL)
}
\arguments{
  \item{Zdst}{matrix of sums of indicators whether or not
  predator ate prey species s during occurrence t; TxS}

  \item{Ydst}{matrix sum of number of caught prey species s
  during occurrence t; TxS}

  \item{lambda}{matrix of parameters representing rates
  predator ate prey species s in time period t; TxS}

  \item{gamma}{matrix of parameters representing rates
  traps caught prey species s in time period t; TxS}

  \item{J}{vector of predators caught in each time period}

  \item{I}{vector of number of days all traps were left out
  in a given time period}

  \item{c}{scalar in null hypotheses}
}
\description{
log-likelihood for non-observed count data; model with EM
}

