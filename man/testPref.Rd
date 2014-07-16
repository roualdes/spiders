\name{testPref}
\alias{testPref}
\title{test predPref function by simulating data and fitting model}
\usage{
testPref(S, T, J, I, lambda, gamma, M = 100, EM = F, n = 4)
}
\arguments{
  \item{S}{number of prey species}

  \item{T}{number of time periods}

  \item{J}{number of predators caught at each time}

  \item{I}{effective number of traps at each time}

  \item{lambda}{matrix of rates at which predator eats prey
  species; TxS}

  \item{gamma}{matrix of rates at which prey species is
  seen in habitat; TxS}

  \item{M}{number of simulated datasets}

  \item{EM}{boolean specifying test of EM algorithm}

  \item{n}{number of parameters to randomly sample; max
  allowed S*T}
}
\description{
test predPref function by simulating data and fitting model
}

