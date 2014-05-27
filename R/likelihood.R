##' log-likelihood of predator preferances model for null, simple, or general hypotheses with balanced data
##'
##' @details this one should be used for asymptotics
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param lambda matrix of parameters representing rates predator ate prey species s in time period t; TxS
##' @param gamma matrix of parameters representing rates traps caught prey species s in time period t; TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param c scalar or vector (indexed by time) in hypotheses
llb <- function(Xdst, Ydst, lambda, gamma, J, I, c=NULL) {
    if (missing(c)) {
        ## general alternative
        out <- -J*sumST(lambda) + log(J)*sumST(Xdst) + sumST(Xdst*log(lambda)) - I*sumST(gamma) + log(I)*sumST(Ydst) + sumST(Ydst*log(gamma))
    } else {
        ## null & simple alternative
        out <- -c*J*sumST(gamma) + log(c*J)*sumST(Xdst) + sumST(Xdst*log(gamma)) - I*sumST(gamma) + sumST(Ydst*log(gamma)) + log(I)*sumST(Ydst)
    }
    out
}

##' complete data log-likelihood for null or simple alternative hypotheses
##'
##' @details this one should be used for optimization
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param lambda matrix of parameters representing rates predator ate prey species s in time period t; TxS
##' @param gamma matrix of parameters representing rates traps caught prey species s in time period t; TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param c scalar or vector (indexed by time) in hypotheses
lcomp <- function(Xdst, Ydst, lambda, gamma, J, I, c) {
    - sumT(c*J*sumSp(gamma)) + sumT(log(c)*sumSp(Xdst)) +
        sumST(Xdst*log(gamma)) - sumT(I*sumSp(gamma)) + sumST(Ydst*log(gamma))        
}

## ## null likelihood used in optimization
## l0 <- function(theta) {
##     c <- theta[1]
##     gamma <- matrix(theta[-1], nrow=T, ncol=S)
##     lcomp(c, gamma, J, ND, Xdst, Ydst)
## }

## ## simple alternative likelihood used in optimization
## l1 <- function(theta) {
##     idx <- 1:T
##     c <- theta[idx]
##     gamma <- matrix(theta[-idx], nrow=T, ncol=S)
##     lcomp(c, gamma, J, ND, Xdst, Ydst)
## }

## ## gradient of null likelihood
## grl0 <- function(theta) {
##     c(-sumT(J*sumSp(gamma)) + sumST(Xdst)/theta[1], # d/dc l0
##       )
## }
