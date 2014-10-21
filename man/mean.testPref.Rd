\name{mean.testPref}
\alias{mean.testPref}
\title{means of \code{testPref}}
\usage{
\method{mean}{testPref}(x, ..., hypothesis = c("null", "alt"))
}
\arguments{
  \item{x}{a testPref object as returned by the eponymous
  function}

  \item{...}{additional arguments}

  \item{hypothesis}{specify which hypothesis}
}
\description{
summarize output of \code{testPref} by calculating means,
across all simulations, for a specified hypothesis
}
\details{
This function essentially calculates means over multiple
replications and fits of a given set of hypotheses and
parameter values, as is output from the function
\code{testPref}.
}
\seealso{
\code{\link{simPref}} \code{\link{testPref}}
\code{\link{plotTestPref}}
}

