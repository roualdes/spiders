##' log-likelihood of predator preferances model with balanced data
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param lambda matrix of parameters representing rates predator ate prey species s in time period t; TxS
##' @param gamma matrix of parameters representing rates traps caught prey species s in time period t; TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param c scalar in null hypotheses
llb <- function(Xdst, Ydst, lambda, gamma, J, I, c=NULL) {
    if ( missing(c) ) {
        ## general alternative
        out <- -J*sumST(lambda) + log(J)*sumST(Xdst) + sumST(Xdst*log(lambda)) 
    } else {
        ## null
        out <- -c*J*sumST(gamma) + log(c*J)*sumST(Xdst) + sumST(Xdst*log(gamma))
    }
    out - I*sumST(gamma) + log(I)*sumST(Ydst) + sumST(Ydst*log(gamma))
}

##' log-likelihood of predator preferances model with unbalanced data
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param lambda matrix of parameters representing rates predator ate prey species s in time period t; TxS
##' @param gamma matrix of parameters representing rates traps caught prey species s in time period t; TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param c scalar in null hypotheses
ll <- function(Xdst, Ydst, lambda, gamma, J, I, c) {
    if ( missing(c) ) {
        ## general alternative
        out <- -sumT(J*sumSp(lambda)) + sumT(log(J)*sumSp(Xdst)) + sumST(Xdst*log(lambda))
    } else {
        ## null
        out <- -c*sumT(J*sumSp(gamma)) + sumT(log(c*J)*sumSp(Xdst)) + sumST(Xdst*log(gamma)) 
    }
    out - sumT(I*sumSp(gamma)) + sumT(log(I)*sumSp(Ydst)) + sumST(Ydst*log(gamma))
}



##' log-likelihood of predators preferances model with EM
##'
##' @param Zdst matrix of sums of indicators whether or not predator ate prey species s during occurrence t; TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; TxS
##' @param lambda matrix of parameters representing rates predator ate prey species s in time period t; TxS
##' @param gamma matrix of parameters representing rates traps caught prey species s in time period t; TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param c scalar in null hypotheses
llEM <- function(Zdst, Ydst, lambda, gamma, J, I, c) {
    if ( missing(c) ) {
        ## avoid log(0) when hat{lambda} = 0
        ## these values contribute nothing, so ignore them
        ## TBD might also need to watch for any hat{gamma} = 0
        w <- which(Zdst != 0, arr.ind=T)
        ## general alternative
        elambda <- exp(lambda[w])
        EX <- lambda[w]*elambda / (elambda - 1)
        out <- -sumT(J*sumSp(lambda)) + sumST(Zdst[w]*log(lambda[w])*EX)
    } else {
        ## null
        cg <- c*gamma
        elambda <- exp(cg)
        EX <- cg*elambda / (elambda - 1)
        out <- -c*sumT(J*sumSp(gamma)) + sumST(Zdst*log(cg)*EX)
    }
    out - sumT(I*sumSp(gamma)) + sumT(log(I)*sumSp(Ydst)) + sumST(Ydst*log(gamma))
}
