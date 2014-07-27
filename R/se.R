##' calculate standard error of EM estimates
##'
##' @param lambda parameter estimates of lambda from EM
##' @param gamma parameter estimates of gamma from EM
##' @param c parameter estimates of c from EM
##' @param Zdst matrix of sums of binary response, whether prey species s was eaten or not during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix of sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @details output is an array where the rownames are indexed by {(T,S)}, T enumerated first
seEM <- function(lambda, gamma, c, Zdst, Ydst, J, I) {

    ## some numbers
    S <- ncol(gamma)
    s <- seq_len(S)
    T <- nrow(gamma)
    t <- seq_len(T)
    ST <- S*T
    indices <- unlist(lapply(s,
                             function(y) sapply(t,
                                                function(x) paste(x, y, sep=''))))    
    if ( is.null(c) ) {                 # H1
        out <- matrix(0, nrow=ST, ncol=2)
        colnames(out) <- c('lambda', 'gamma')
        out[,'gamma'] <- unlist(Ydst/gamma^2)
        elambda <- exp(lambda)
        out[,'lambda'] <- unlist(J*elambda/(1-elambda)^2)
        rownames(out) <- indices
    } else {                            # H0
        out <- rep(0, ST+1)
        l <- c*gamma
        el <- exp(l)
        tmp <- J*el/(1-el)^2
        g2 <- gamma^2
        out[1] <- sumST(tmp*g2)
        out[-1] <- unlist(tmp*c^2 + Ydst/g2)
        names(out) <- c('c', paste('gamma', indices, sep=''))
    }
    sqrt(1/out)
}


##' calculate standard error of estimates of predator preferences model
##'
##' @param lambda parameter estimates of lambda from EM
##' @param gamma parameter estimates of gamma from EM
##' @param c parameter estimates of c from EM
##' @param Xdst matrix of sum of prey species s eaten during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix of sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @details output is an array where the rownames are indexed by {(T,S)}, T enumerated first
se <- function(lambda, gamma, c, Xdst, Ydst, J, I) {

    ## some numbers
    S <- ncol(gamma)
    s <- seq_len(S)
    T <- nrow(gamma)
    t <- seq_len(T)
    ST <- S*T
    indices <- unlist(lapply(s,
                             function(y) sapply(t,
                                                function(x) paste(x, y, sep=''))))
    if ( is.null(c) ) {                 # H1
        out <- matrix(0, nrow=ST, ncol=2)
        colnames(out) <- c('lambda', 'gamma')
        out[,'gamma'] <- unlist(Ydst/gamma^2)
        out[,'lambda'] <- unlist(Xdst/lambda^2)
        rownames(out) <- indices
    } else {
        out <- rep(0, ST+1)
        out[1] <- sumST(Xdst/c^2)
        out[-1] <- unlist((Xdst+Ydst)/gamma^2)
        names(out) <- c('c', paste('gamma', indices, sep=''))
    }
    sqrt(1/out)
}
