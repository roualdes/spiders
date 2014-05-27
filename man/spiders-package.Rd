\name{spiders-package}
\alias{spiders-package}
\alias{spiders}
\docType{package}
\title{ Fits predator preferences model. }
\description{
Models predator preferences over an array of time and prey species.
}
\details{
\tabular{ll}{
Package: \tab spiders\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-30\cr
License: \tab GPL \cr
}

}
\author{
Edward A. Roualdes <edward.roualdes@uky.edu>
}
\references{
Reliability Estimates for Ivlev's Electivity Index, the Forage Ratio, and a Proposed Linear Index of Food Selection
Richard E. Strauss 
Transactions of the American Fisheries Society 
Vol. 108, Iss. 4, 1979
}

\keyword{ spiders }
\seealso{
}
\examples{
## make up some numbers
Predators <- 20
Traps <- 20                              
PreySpecies <- 3
Times <- 12
ST <- Times*PreySpecies
l <- matrix(1:ST, nrow=Times, ncol=PreySpecies)
g <- matrix(2*(1:ST), nrow=Times, ncol=PreySpecies)

##  simulate data
fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g)

## calculate model
(prefs <- predPref(fdata$eaten, fdata$caught))

}
