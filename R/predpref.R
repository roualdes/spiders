##' estimates parameters of predator preferences model and calculates LRT
##'
##' @param eaten a df of eatings preferences; TxS
##' @param caught a df of caught prey species; TxS
##' @param hypotheses a 2-tuple specifying the null and alternative hypotheses, respectively
##' @param alpha LRT level of significance
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##' @export
predPref <- function(eaten, caught, hypotheses = c('c', 'Ct'), alpha=0.05, em_maxiter=1000) {

    ## check hypotheses specification
    hypotheses <- checkHypotheses(hypotheses)

    ## data
    dname <- paste(deparse(substitute(eaten)), 'and', deparse(substitute(caught)))
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

    ## predators (J), traps (I)
    J <- getTimeCounts(eaten, 'adj')[,2]
    I <- getTimeCounts(caught, 'adj')[,2] # total days traps were out each t
    
    ## data for calculations
    Xdst <- getTimeCounts(eaten, preyNames)[,preyNames, drop=F]
    Ydst <- getTimeCounts(caught, preyNames)[,preyNames, drop=F]

    ## do we run EM?
    EM <- ifelse(!any(X > 1), TRUE, FALSE)
    
    ## are data balanced
    BAL <- length(unique(J)) == 1 && length(unique(I)) == 1

    ## errors with time points
    if ( nrow(Xdst) != nrow(Ydst) ) stop("Differing number of time points in eaten/caught data.")

    calcs <- calcHypotheses(hyp = hypotheses,
                      Xdst = Xdst, Ydst = Ydst, J=J, I=I,
                      balanced = BAL, EM=EM, em_maxiter = em_maxiter)
    llH0 <- calcs$llH0; llH1 <- calcs$llH1
    null <- calcs$null; alt <- calcs$alt
    df <- calcs$df

    ## LRT stat
    Lambda <- -2*(llH0 - llH1)
    lrt <- list(Lambda = Lambda, df = df,
                p.value = pchisq(Lambda, df=df, lower.tail=F))
    
    out <- list('alt' = alt, 'null' = null,
                'loglikH1' = llH1, 'loglikH0' = llH0,
                'numPredators' = J, 'numTraps' = I,
                LRT = lrt, hypotheses = hypotheses, data.name = dname)
    class(out) <- 'predPref'
    out
}
