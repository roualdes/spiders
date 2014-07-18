##' calculate second derivatives of Q, H
##'
##' @param lambda parameter estimates of lambda from EM
##' @param gamma parameter estimates of gamma from EM
##' @param c parameter estimates of c from EM
##' @param Zdst matrix of sums of binary response, whether prey species s was eaten or not during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
hess <- function(lambda, gamma, c, Zdst, Ydst, J, I) {
    
    ## some numbers
    S <- ncol(lambda); s <- seq_len(S)
    T <- nrow(lambda); t <- seq_len(T)
    ST <- S*T; st <- seq_len(st)
    elambda <- exp(lambda)
    EX <- lambda*elambda / (elambda - 1)
    ZEX <- Zdst*EX

    if (missing(c)) {
        ## fill negative Q
        Q <- rep(0, 2*ST)
        Q[st] <- unlist(ZEX/lambda^2)
        Q[-st] <- unlist(Ydst/gamma^2)

        ## fill H
        H <- rep(0, 2*ST)
        enlam <- exp(-lambda)
        H[st] <- unlist(J*(enlam/(1-enlam) + enlam*enlam / (1-enlam)^2))
    } else {
        
    }
}
