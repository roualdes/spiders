##' log-likelihood of predator preferances model
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
        out <- -sumT(J*sumSp(lambda)) + sumT(log(J)*sumSp(Xdst)) + sumST(Xdst*log(lambda))
    } else {
        lc <- length(c)
        if ( lc == nrow(Xdst) || lc == 1 ) {
            out <- -sumT(c*J*sumSp(gamma)) + sumT(log(c*J)*sumSp(Xdst)) + sumST(Xdst*log(gamma)) 
        } else if ( lc == ncol(Xdst) ) {
            s <- seq_len(ncol(Xdst))
            cGamma <- sapply(s, function(j) gamma[,j]*c[j])
            out <- -sumT(J*sumSp(cGamma)) + sumT(log(J)*sumSp(Xdst)) + sumST(Xdst*log(cGamma))
        } else {
            out <- -sumT(J*sumSp(c*gamma)) + sumT(log(J)*sumSp(Xdst)) + sumST(Xdst*log(c*gamma))
        }
    }
    out - sumT(I*sumSp(gamma)) + sumST(Ydst*log(gamma))
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
## llEM <- function(Zdst, Ydst, lambda, gamma, J, I, c=NULL) {
##     if ( is.null(c) ) {
##         ## avoid log(0) when hat{lambda} = 0
##         ## these values contribute nothing, so ignore them
##         ## TBD might also need to watch for any hat{gamma} = 0
##         w <- which(Zdst != 0, arr.ind=T)
##         out <- sumST(Zdst[w]*log(1-exp(-lambda[w])) - (J[w[,1]] - Zdst[w])*lambda[w])
##     } else {
##         lc <- length(c)
##         if ( lc == 1 ) {
##             out <- sumST(Zdst*log(1-exp(-c*gamma))) - c*sumST((J-Zdst)*gamma)
##         } else if ( lc == nrow(Zdst)  ) {
##             out <- sumST(Zdst*log(1-exp(-c*gamma))) - sumT(c*sumSp((J-Zdst)*gamma))
##         } else if ( lc == ncol(Zdst) ){
##             s <- seq_len(ncol(Zdst))
##             cGamma <- sapply(s, function(j) gamma[,j]*c[j])
##             out <- sumST(Zdst*log(1-exp(-cGamma))) - sumST((J-Zdst)*cGamma)
##         } else {
##             out <- sumST(Zdst*log(1-exp(-c*gamma))) - sumST((J-Zdst)*c*gamma)
##         }
##     }
##     out - sumT(I*sumSp(gamma)) + sumST(Ydst*log(gamma))
## }
llEM <- function(Zdst, Ydst, lambda, gamma, J, I, c=NULL) {
    if ( is.null(c) ) {
        l <- lambda
    } else {
        lc <- length(c)
        if ( lc ==  ncol(Zdst) ) {
            l <- sapply(seq_len(lc), function(j) gamma[,j]*c[j])            
        } else {
            l <- c*gamma
        }
    }
    g <- gamma*I
    out <- 0
    for ( i in seq_len(nrow(Zdst)) ) {
        for ( j in seq_len(ncol(Zdst)) ) {
            out <- out + dpois(Ydst[i,j], lambda=g[i,j], log=TRUE) +
                dbinom(Zdst[i,j], size=J[i], prob=1-exp(-l[i,j]), log=TRUE)
        }
    }
    out
}
