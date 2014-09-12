\name{testC}
\alias{testC}
\title{estimates linear contrasts of c, c_s, c_t, c_st}
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
\description{
estimates linear contrasts of c, c_s, c_t, c_st
}

