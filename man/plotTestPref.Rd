\name{plotTestPref}
\alias{plotTestPref}
\title{plot testPref}
\usage{
plotTestPref(x, hypothesis = c("null", "alt"), lambda = NULL,
  gamma = NULL)
}
\arguments{
  \item{x}{a testPref object as returned by the eponymous
  function}

  \item{hypothesis}{specify which hypothesis to plot}

  \item{lambda}{a matrix of true values of the parameter
  lambda; TxS}

  \item{gamma}{a matrix of true values of the parameter
  gamma; TxS}
}
\description{
plot the output from \code{testPref}
}
\details{
If true values of lambda and gamma are given, then tick
(rugs) marks will be placed at those values along the
x-axis.
}
\seealso{
\code{\link{simPref}} \code{\link{testPref}}
\code{\link{mean.testPref}}
}

