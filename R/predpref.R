##' estimates parameters of predator preferences model and calculates LRT
##'
##' @param eaten a df of eatings preferences; TxS
##' @param caught a df of caught prey species; TxS
##' @param alpha LRT level of significance
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##' @export
predPref <- function(eaten, caught, alpha=0.05, em_maxiter=100) {

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
    EM <- FALSE                         # assume no
    if ( !any(X > 1) ) EM <- TRUE    # double check

    ## predators (J), traps (I), prey species (S), times (T)
    S <- length(preyNames)
    J <- getMonthCounts(eaten, 'adj')[,2]
    I <- getMonthCounts(caught, 'adj')[,2] # total days traps were out each t
    Xdst <- getMonthCounts(eaten, preyNames)[,preyNames]
    Ydst <- getMonthCounts(caught, preyNames)[,preyNames]
    T <- nrow(Xdst)                  # assuming same times in both X,Y

    ## errors with time points
    if ( nrow(Xdst) != nrow(Ydst) ) stop("Differing number of time points in eaten/caught data.")
    
    ## estimate parameters
    ## balanced data
    if ( length(unique(J)) == 1 && length(unique(I)) == 1 ) {

        ## balanced data
        null <- est0b(Xdst, Ydst, J[1], I[1])
        gAlt <- est1b(Xdst, Ydst, J[1], I[1])

        ## calc likelihoods
        llH0 <- llb(Xdst, Ydst, NA, null$gamma, J[1], I[1], null$c)
        llH1 <- llb(Xdst, Ydst, gAlt$lambda, gAlt$gamma, J[1], I[1])
        
    } else {                            # not balanced
        if ( EM ) {                     # EM?
            null <- estEM0(Xdst, Ydst, J, I, em_maxiter)
            gAlt <- estEM1(Xdst, Ydst, J, I, em_maxiter)

            ## calc likelihoods
            llH0 <- llEM(Xdst, Ydst, NA, null$gamma, J, I, null$c)
            llH1 <- llEM(Xdst, Ydst, gAlt$lambda, gAlt$gamma, J, I)
            
        } else {                        # unbalanced
            null <- est0(Xdst, Ydst, J, I)
            gAlt <- est1(Xdst, Ydst, J, I)

            ## calc likelihoods
            llH0 <- ll(Xdst, Ydst, NA, null$gamma, J, I, null$c)
            llH1 <- ll(Xdst, Ydst, gAlt$lambda, gAlt$gamma, J, I)
        }
    }

    ## LRT stats
    ## calculate degrees of freedom on asymptotic chi-squared
    df <- S*T-1
    LRT <- -2*(llH0 - llH1)
    
    list('gAlt' = gAlt, 'null' = null,
        'loglikH1' = llH1, 'loglikH0' = llH0,
         '-2log(LRT)' = LRT, 'df' = df, 'p.value' = pchisq(LRT, df=df, lower.tail=F))
}
