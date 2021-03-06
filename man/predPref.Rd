% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/predPref.R
\name{predPref}
\alias{predPref}
\title{parameter estimation}
\usage{
predPref(eaten, caught, hypotheses = c("c", "Ct"), alpha = 0.05,
  em_maxiter = 1000)
}
\arguments{
\item{eaten}{a dataframes of eatings preferences; TxS}

\item{caught}{a dataframes of caught prey species; TxS}

\item{hypotheses}{a 2-tuple specifying the null and alternative hypotheses, respectively}

\item{alpha}{LRT level of significance}

\item{em_maxiter}{maximum number of iterations allowed for EM algorithm}
}
\value{
A list of class 'predPref' with the following elements:

null: parameters as estimated under the specified null hypothesis.

alt: parameters as estimated under the specified alternative hypothesis.

loglikH0: the null hypothesis log-likelihood, with constants not accounted for.

loglikH1: the alternative hypothesis log-likelihood, with constants not accounted for.

J: a column vector of dimension T containing the number of predators in
each time period.

I: a column vector of dimension T containing the number of traps in each time period.

LRT: the likelihood ratio test statistics.

hypotheses: a 2-tuple of the user specified hypotheses.

data.name: a character string giving the names of the data.
}
\description{
Estimates parameters of predator preferences model and calculates LRT.
Eaten and caught dataframes are indexed with rows across time points
and columns of prey species.
}
\examples{
# set parameters
Predators <- Traps <- 100
PreySpecies <- 2
Times <- 5
g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)     # gamma
l <- matrix(seq(0.4,1.8,length.out=5)*sqrt(2), nrow=Times, ncol=PreySpecies) # ct

# fit model
\dontrun{
fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=FALSE)
predPref(fdata$eaten, fdata$caught, hypotheses=c('ct', 'cst'))
}
}
\seealso{
\code{\link{simPref}} \code{\link{summary.predPref}}
}

