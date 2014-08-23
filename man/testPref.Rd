\name{testPref}
\alias{testPref}
\title{test predPref function by simulating data and fitting model}
\usage{
testPref(J, I, lambda, gamma, M = 100, hyp = c("C", "Cst"), EM = FALSE,
  em_maxiter = 100)
}
\arguments{
  \item{J}{number of predators caught at each time}

  \item{I}{effective number of traps at each time}

  \item{lambda}{matrix of rates at which predator eats prey
  species; TxS}

  \item{gamma}{matrix of rates at which prey species is
  seen in habitat; TxS}

  \item{M}{number of simulated datasets}

  \item{hyp}{a 2-tuple specifying the null and alternative
  hypotheses, respectively}

  \item{EM}{boolean specifying test of EM algorithm}

  \item{em_maxiter}{maximum number of iterations allowed
  for EM algorithm}
}
\description{
test predPref function by simulating data and fitting model
}

