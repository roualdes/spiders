\name{testPref}
\alias{testPref}
\title{simulate and test predPref}
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
\value{
A list of two elements, one for each hypothesis. Each
element contains parameter estimates as calculated by the
respective hypothesis, and the number of iterations used to
fit said model.
}
\description{
simulate data and test function \code{predPref} on each
simulated dataset
}
\examples{
# set parameters
Predators <- Traps <- 100
PreySpecies <- 2
Times <- 5
g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)     # gamma
l <- matrix(seq(0.4,1.8,length.out=5)*sqrt(2), nrow=Times, ncol=PreySpecies) # ct

# test functions
\dontrun{
testPref(Predators, Traps, l, g, M=10, hyp=c('c', 'cst'), EM=F)
}
}
\seealso{
\code{\link{simPref}} \code{\link{plotTestPref}}
\code{\link{mean.testPref}}
}

