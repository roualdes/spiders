##' function to calculate hypotheses given a user specifed null and alternative
##'
##' @param hyp a 2-tuple specifying the null and alternative hypotheses, respectively
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param balanced boolean specifying balanced data or not
##' @param EM boolean specifying if EM algorithm should be used
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
calcHypotheses <- function(hyp, Xdst, Ydst, J, I, balanced, EM, em_maxiter) {

    ## EM?
    if ( EM ) {
        ## null
        null <- estEM0(Xdst, Ydst, J, I, ifelse(hyp[1]=='index_c', TRUE, FALSE), em_maxiter)
        null$SE <- seEM(NULL, null$gamma, null$c, Xdst, Ydst, J, I)
        llH0 <- llEM(Xdst, Ydst, NA, null$gamma, J, I, null$c)

        ## alt
        if ( hyp[2] == 'general' ) {
            alt <- estEM1(Xdst, Ydst, J, I, em_maxiter)
            alt$SE <- seEM(alt$lambda, alt$gamma, NULL, Xdst, Ydst, J, I)
            llH1 <- llEM(Xdst, Ydst, alt$lambda, alt$gamma, J, I)
        } else {
            alt <- estEM0(Xdst, Ydst, J, I, TRUE, em_maxiter)
            alt$SE <- seEM(NULL, alt$gamma, alt$c, Xdst, Ydst, J, I)
            llH1 <- llEM(Xdst, Ydst, NA, alt$gamma, J, I, alt$c)
        }
    } else {

        ## null
        if ( balanced && hyp[1] == 'c') {
            ## solve analytically
            null <- est0b(Xdst, Ydst, J[1], I[1])
            null$SE <- se(NULL, null$gamma, null$c, Xdst, Ydst, J, I)
            llH0 <- llb(Xdst, Ydst, NA, null$gamma, J[1], I[1], null$c)
        } else {
            ## even if c is not indexed, data not balanced => need iterative solution
            null <- est0(Xdst, Ydst, J, I, ifelse(hyp[1] == 'index_c', TRUE, FALSE))
            null$SE <- se(NULL, null$gamma, null$c, Xdst, Ydst, J, I)
            llH0 <- ll(Xdst, Ydst, NA, null$gamma, J, I, null$c)
        }
        
        
        ## alt
        if ( hyp[2] == 'general' ) {
            ## general alternative
            alt <- est1(Xdst, Ydst, J, I)
            alt$SE <- se(alt$lambda, alt$gamma, NULL, Xdst, Ydst, J, I)
            llH1 <- ll(Xdst, Ydst, alt$lambda, alt$gamma, J, I)            
        } else {
            ## c is indexed
            alt <- est0(Xdst, Ydst, J, I, TRUE)
            alt$SE <- se(NULL, alt$gamma, alt$c, Xdst, Ydst, J, I)
            llH1 <- ll(Xdst, Ydst, NA, alt$gamma, J, I, alt$c)
        }
    }
    
    ## calculate degrees of freedom
    S <- ncol(Xdst); T <- nrow(Xdst); ST <- S*T
    nullDF <- ST + ifelse( hyp[1] == 'index_c', T, 1)
    altDF <- ST + ifelse( hyp[2] == 'general', ST, T)
    df <- altDF - nullDF
    
    list('llH0' = llH0, 'llH1' = llH1, 'null' = null, 'alt' = alt, 'df' = df)
}

##' function to check user specified hypotheses
##'
##' @param hyp a 2-tuple specifying the null and alternative hypotheses, respectively
checkHypotheses <- function(hyp) {
    
    ## initialize output
    H <- rep(0, 2)
    
    ## H0
    if ( grepl('index', hyp[1]) ) {
        H[1] <- 'index_c'
    } else if ( grepl('c', hyp[1]) ) {
        H[1] <- 'c'
    } else {
        stop('Hypotheses specified incorrectly; please check documentation.')
    }

    ## H1
    if ( grepl('gen', hyp[2]) ) {
        H[2] <- 'general'
    } else if ( grepl('index', hyp[2]) ) {
        H[2] <- 'index_c'
    }

    ## both
    if ( H[1] == H[2] ) stop('Null and Alternative hypotheses must be different.')

    H
}
