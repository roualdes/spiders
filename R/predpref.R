##' estimates parameters of predator preferences model and calculates LRT
##'
##' @param eaten a df of eatings preferences; TxS
##' @param caught a df of caught prey species; TxS
##' @param alpha LRT level of significance
##' @export
predPref <- function(eaten, caught, alpha=0.05) {

    ## preliminary
    xNames <- colnames(eaten)
    yNames <- colnames(caught)
    if ( !('adj' %in% xNames) ) eaten$adj <- 1
    if ( !('adj' %in% yNames) ) caught$adj <- 1

    extraVars <- c('time', 'adj')
    preyNames <- setdiff(xNames, extraVars)

    ## predators (J), traps (I), prey species (S), times (T)
    S <- length(preyNames)
    J <- getMonthCounts(eaten, 'adj')[,2]
    I <- getMonthCounts(caught, 'adj')[,2] # total days traps were out each t
    Xdst <- getMonthCounts(eaten, preyNames)[,preyNames]
    Ydst <- getMonthCounts(caught, preyNames)[,preyNames]
    T <- nrow(Xdst)                  # assuming same times in both X,y

    ## estimate parameters
    if ( length(unique(J)) == 1 && length(unique(I)) == 1 ) {

        ## balanced data
        null <- est0b(Xdst, Ydst, J, I)
        gAlt <- est1b(Xdst, Ydst, J, I)

        ## calc likelihoods
        llH0 <- llb(Xdst, Ydst, NA, null$gamma, J[1], I[1], null$c)
        llH1 <- llb(Xdst, Ydst, gAlt$lambda, gAlt$gamma, J[1], I[1])

        ## calculate degrees of freedom on asymptotic chi-squared
        df <- S*T-1
        
    } else {
        stop('unbalanced data is not yet implemented.')
    }

    ## lrt stat
    LRT <- -2*(llH0 - llH1)
    
    list('loglikH1' = llH1, 'loglikH0' = llH0,
         '-2log(LRT)' = LRT, 'df' = df, 'p.value' = pchisq(LRT, df=df, lower.tail=F))
}
