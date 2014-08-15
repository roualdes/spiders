##' estimates parameters of predator preferences model and calculates LRT
##'
##' @param eaten a df of eatings preferences; TxS
##' @param caught a df of caught prey species; TxS
##' @param alpha LRT level of significance
##' @param index_c boolean; TRUE indexes c in null hypothesis by t
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##' @export
predPref <- function(eaten, caught, alpha=0.05, index_c = TRUE, em_maxiter=100) {

    xNames <- colnames(eaten)
    yNames <- colnames(caught)

    ## traps left out for different numbers of days?
    if ( !('adj' %in% xNames) ) eaten$adj <- 1
    if ( !('adj' %in% yNames) ) caught$adj <- 1

    ## prey names
    extraVars <- c('time', 'adj')
    preyNames <- setdiff(xNames, extraVars)

    ## data errors
    X <- eaten[,preyNames]
    Y <- caught[,preyNames]
    if ( any(X < 0) || any(Y < 0) ) stop("Count data can not be less than 0.")

    ## do we run EM?
    EM <- ifelse(!any(X > 1), TRUE, FALSE)

    ## data for calculations
    Xdst <- getTimeCounts(eaten, preyNames)[,preyNames, drop=F]
    Ydst <- getTimeCounts(caught, preyNames)[,preyNames, drop=F]

    ## errors with time points
    if ( nrow(Xdst) != nrow(Ydst) ) stop("Differing number of time points in eaten/caught data.")

    ## predators (J), traps (I), prey species (S), times (T)
    J <- getTimeCounts(eaten, 'adj')[,2]
    I <- getTimeCounts(caught, 'adj')[,2] # total days traps were out each t
    S <- length(preyNames)
    T <- nrow(Xdst)                  # assuming same times in both X,Y

    ## EM?
    if ( EM ) {
        ## estimate parameters
        null <- estEM0(Xdst, Ydst, J, I, index_c, em_maxiter)
        alt <- estEM1(Xdst, Ydst, J, I, em_maxiter)

        ## standard errors
        null$SE <- seEM(NULL, null$gamma, null$c, Xdst, Ydst, J, I)
        alt$SE <- seEM(alt$lambda, alt$gamma, NULL, Xdst, Ydst, J, I)
        
        ## calc likelihoods
        llH0 <- llEM(Xdst, Ydst, NA, null$gamma, J, I, null$c)
        llH1 <- llEM(Xdst, Ydst, alt$lambda, alt$gamma, J, I)
    } else {

        ## balanced data
        if ( length(unique(J)) == 1 && length(unique(I)) == 1 && index_c == FALSE ) {
            ## estimate parameters
            null <- est0b(Xdst, Ydst, J[1], I[1])
            alt <- est1b(Xdst, Ydst, J[1], I[1])

            ## standard errors
            null$SE <- se(NULL, null$gamma, null$c, Xdst, Ydst, J, I)
            alt$SE <- se(alt$lambda, alt$gamma, NULL, Xdst, Ydst, J, I)

            ## calc likelihoods
            llH0 <- llb(Xdst, Ydst, NA, null$gamma, J[1], I[1], null$c)
            llH1 <- llb(Xdst, Ydst, alt$lambda, alt$gamma, J[1], I[1])
            
        } else {                        # not balanced
            ## estimate parameters
            null <- est0(Xdst, Ydst, J, I, index_c)
            alt <- est1(Xdst, Ydst, J, I)

            ## standard errors
            null$SE <- se(NULL, null$gamma, null$c, Xdst, Ydst, J, I)
            alt$SE <- se(alt$lambda, alt$gamma, NULL, Xdst, Ydst, J, I)

            ## calc likelihoods
            llH0 <- ll(Xdst, Ydst, NA, null$gamma, J, I, null$c)
            llH1 <- ll(Xdst, Ydst, alt$lambda, alt$gamma, J, I)
        }
    }

    ## LRT stats
    ## calculate degrees of freedom on asymptotic chi-squared
    df <- S*T-length(c)
    Lambda <- -2*(llH0 - llH1)
    
    out <- list('alt' = alt, 'null' = null,
                'loglikH1' = llH1, 'loglikH0' = llH0,
                'numPredators' = J, 'numTraps' = I,
                'Lambda' = Lambda, 'df' = df,
                'p.value' = pchisq(Lambda, df=df, lower.tail=F))
    class(out) <- 'predPref'
    out
}
