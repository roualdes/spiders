\name{simulate}
\alias{simulate}
\title{simulate data for predator preferences model}
\usage{
simulate(S, T, J, I, lambda, gamma)
}
\arguments{
  \item{S}{number of prey species}

  \item{T}{number of times}

  \item{J}{number of predators caught at each time; assumed
  constant for now}

  \item{I}{number of traps at each time; assumed constant
  for now}

  \item{lambda}{vector of rates at which predator eats prey
  species, indexed by s and t}

  \item{gamma}{vector of rates at which prey species is
  seen in habitat, indexed by s and t}
}
\description{
simulate data for predator preferences model
}

