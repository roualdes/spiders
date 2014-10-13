\name{testC}
\alias{testC}
\title{linear contrast of c_st}
\usage{
testC(x, b, mu = 0, alternative = c("two.sided", "less", "greater"),
  conf.level = 0.95)
}
\arguments{
  \item{x}{a prefPref object as fit by the eponymous
  function}

  \item{b}{a vector to linearly transform c_st}

  \item{mu}{a number to test the linear contrast against in
  the null}

  \item{alternative}{string to specify alternative
  hypothesis}

  \item{conf.level}{confidence level of the interval}
}
\value{
A list with class '"htest"' containing the following
components:

statistic: the value of the t-statistic.

parameter: the degrees of freedom for the t-statistic.

p.value: the p-value for the test.

conf.int: a confidence interval for the mean appropriate to
the specified alternative hypothesis.

estimate: the estimated mean or difference in means
depending on whether it was a one-sample test or a
two-sample test.

null.value: the specified hypothesized value of the mean or
mean difference depending on whether it was a one-sample
test or a two-sample test.

alternative: a character string describing the alternative
hypothesis.

method: a character string indicating what type of t-test
was performed.

data.name: a character string giving the names of the data.
}
\description{
estimates linear contrasts of the elements of c, c_s, c_t,
c_st
}
\details{
The input vector b performs the linear transformation t(b)
\%*\% matrix(c_st), so that c_st becomes a column vector by
indexing t first and then s.  Hence there is no requirement
of a linear contrast, a linear transformation such that
t(b) \%*\% matrix(1, nrow=length(b)) != 0 is allowed.
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
pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('ct', 'cst'))
testC(pref, b = c(0,1, -1, 0, 0))
}
}

