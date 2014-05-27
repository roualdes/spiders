##' estimates parameters from the null model when balanced data
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est0b <- function(Xdst, Ydst, J, I) {
    i <- I[1]
    iXdY <- i*sumST(Xdst)/sumST(Ydst)
    gammaHat <- (Xdst + Ydst) / (iXdY + i)
    cHat <- iXdY/J[1]
    return(list('gamma' = gammaHat, 'c' = cHat))
}

##' estimates parameters from the general alternative model when balanced data
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est1b <- function(Xdst, Ydst, J, I) {
    gammaHat <- Ydst/I[1]
    lambdaHat <- Xdst/J[1]
    return(list('gamma' = gammaHat, 'lambda' = lambdaHat))
}
