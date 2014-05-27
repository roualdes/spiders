\name{simPref}
\alias{simPref}
\title{simulate data for predator preferences model}
\usage{
simPref(S, T, J, I, lambda, gamma)
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
}
\description{
simulate data for predator preferences model
}

