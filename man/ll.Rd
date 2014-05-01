\name{ll}
\alias{ll}
\title{log-likelihood of predator preferances model}
\usage{
ll(Xdst, Ydst, lambda, gamma, J, I, c = NULL)
}
\arguments{
  \item{Xdst}{sum of number of eaten prey species s during
  occurrence t}

  \item{Ydst}{sum of number of caught prey species s during
  occurrence t}

  \item{lambda}{vector of parameters relative to eaten
  prey, indexed by s and t}

  \item{gamma}{vector of parameters relative to caught
  prey, indexed by s and t}

  \item{J}{number of predators caught at each time; assumed
  constant for now}

  \item{I}{number of traps at each time; assumed constant
  for now}

  \item{c}{scalar in our base null hypothesis}
}
\description{
log-likelihood of predator preferances model
}

