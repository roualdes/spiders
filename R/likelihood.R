##' log-likelihood of predator preferances model with unbalanced data
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time and cols indexed by prey species, TxS
##' @param lambda matrix of parameters representing rates predator ate prey species s in time period t; TxS
##' @param gamma matrix of parameters representing rates traps caught prey species s in time period t; TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param c scalar in null hypotheses
ll <- function(Xdst, Ydst, lambda, gamma, J, I, c=NULL) {
    if ( is.null(c) ) {
        ## general alternative
        out <- -sumT(J*sumSp(lambda)) + sumT(log(J)*sumSp(Xdst)) + sumST(Xdst*log(lambda))
    } else {
        lc <- length(c)
        ## null
        if ( lc == nrow(Xdst) || lc == 1 ) {
            out <- -1*sumT(c*J*sumSp(gamma)) + sumT(log(c*J)*sumSp(Xdst)) + sumST(Xdst*log(gamma)) 
        } else {
            s <- seq_len(ncol(Xdst))
            cGamma <- sapply(s, function(j) gamma[,j]*c[j])
            out <- -1*sumT(J*sumSp(cGamma)) + sumT(log(J)*sumSp(Xdst)) + sumST(Xdst*log(cGamma))           
        }
        
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
llEM <- function(Zdst, Ydst, lambda, gamma, J, I, c=NULL) {
    if ( is.null(c) ) {
        ## avoid log(0) when hat{lambda} = 0
        ## these values contribute nothing, so ignore them
        ## TBD might also need to watch for any hat{gamma} = 0
        w <- which(Zdst != 0, arr.ind=T)
        ## general alternative
        out <- sumST(Zdst[w]*log(1-exp(-lambda[w])) - (J[w[,1]]-Zdst[w])*lambda[w])
    } else {
        lc <- length(c)
        ## null
        if ( lc == nrow(Zdst) || lc == 1 ) {
            cGamma <- c*gamma
        } else {
            s <- seq_len(ncol(Zdst))
            cGamma <- sapply(s, function(j) gamma[,j]*c[j])
        }
        out <- sumST(Zdst*log(1-exp(-cGamma)) - (J-Zdst)*cGamma)
    }
    out - sumT(I*sumSp(gamma)) + sumT(log(I)*sumSp(Ydst)) + sumST(Ydst*log(gamma))
}
