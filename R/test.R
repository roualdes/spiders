##' test predPref function by simulating data and fitting model
##'
##' @param J number of predators caught at each time
##' @param I effective number of traps at each time
##' @param lambda matrix of rates at which predator eats prey species; TxS
##' @param gamma matrix of rates at which prey species is seen in habitat; TxS
##' @param M number of simulated datasets
##' @param hyp a 2-tuple specifying the null and alternative hypotheses, respectively
##' @param EM boolean specifying test of EM algorithm
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##' @export
testPref <- function(J, I, lambda, gamma, M=100, hyp=c('C', 'Cst'), EM = FALSE, em_maxiter = 100) {

    ## initialize output structure
    out <- vector('list', 2)
    names(out) <- c('null', 'alt')
    need_set_output <- TRUE             # need initialize storage within output?

    ## some numbers 
    S <- ncol(lambda)
    T <- nrow(lambda)
    
    for ( m in seq_len(M) ) {
        
        ## simulate data and fit model
        fdata <- simPref(S, T, J, I, lambda, gamma, EM=EM)
        prefs <- predPref(fdata$eaten, fdata$caught, hypotheses = hyp, em_maxiter = em_maxiter)

        ## initialize storage within output structure
        if ( need_set_output ) {
            ST <- length(unlist(prefs$null$gamma))
            out$null[['gamma']] <- out$alt[['gamma']] <- matrix(0, M, ST)
            
            ## null model storage
            if ( !is.null(prefs$null$c) ) 
                out$null[['c']] <- matrix(0, M, length(prefs$null$c))
            if ( !is.null(prefs$null$lambda) )
                out$null[['lambda']] <- matrix(0, M, ST)
            if ( !is.null(prefs$null$em_iters) )
                out$null[['em_iters']] <- rep(0, M)
            if ( !is.null(prefs$null$iters) )
                out$null[['iters']] <- rep(0, M)

            ## alt model storage
            if ( !is.null(prefs$alt$c) )
                out$alt[['c']] <- matrix(0, M, length(prefs$alt$c))
            if ( !is.null(prefs$alt$lambda) )
                out$alt[['lambda']] <- matrix(0, M, ST)            
            if ( !is.null(prefs$alt$em_iters) )
                out$alt[['em_iters']] <- rep(0, M)
            if ( !is.null(prefs$alt$iters) )
                out$alt[['iters']] <- rep(0, M)
            
            need_set_output <- FALSE
        }

        ## store output
        out$null$gamma[m,] <- unlist(prefs$null$gamma)
        out$alt$gamma[m,] <- unlist(prefs$alt$gamma)

        ## null
        if ( !is.null(prefs$null$c) )
            out$null$c[m,] <- prefs$null$c
        if ( !is.null(prefs$null$lambda) )
            out$null$lambda[m,] <- unlist(prefs$null$lambda)
        if ( !is.null(prefs$null$em_iters) )
                out$null$em_iters <- prefs$null$em_iters
        if ( !is.null(prefs$null$iters) )
            out$null[['iters']] <- prefs$null$iters
        
        ## alt
        if ( !is.null(prefs$alt$c) )
            out$alt$c[m,] <- prefs$alt$c
        if ( !is.null(prefs$alt$lambda) )
            out$alt$lambda[m,] <- unlist(prefs$alt$lambda)
        if ( !is.null(prefs$alt$iters) )
            out$alt$iters <- prefs$alt$iters
    }
    class(out) <- 'testPref'
    out
}

##' plot the output of testPref
##'
##' @param x a testPref object as returned by the eponymous function
##' @param hypothesis specify which hypothesis to plot
##' @export
plotTestPref <- function(x, hypothesis = 'null') {
    
    ## some numbers
    ST <- ncol(x$null$gamma); st <- seq_len(ST)

    ## data
    nullGamma <- data.frame('gamma' = as.vector(x$null$gamma), 'index' = st)
    altGamma <- data.frame('gamma' = as.vector(x$alt$gamma), 'index' = st)

    ## more null data
    if ( !is.null(x$null$c) )
        nullC <- data.frame('c' = as.vector(x$null$c), 'index' = seq_len(ncol(x$null$c)))
    if ( !is.null(x$null$lambda) )
        nullLambda <- data.frame('lambda' = as.vector(x$null$lambda), 'index' = st)

    ## more alt data
    if ( !is.null(x$alt$c) )
        altC <- data.frame('c' = as.vector(x$alt$c), 'index' = seq_len(ncol(x$alt$c)))
    if ( !is.null(x$alt$lambda) )
        altLambda <- data.frame('lambda' = as.vector(x$alt$lambda), 'index' = st)

    ## plot data
    require(lattice)
    require(gridExtra)
    
    two <- densityplot(~gamma, data=nullGamma, groups=index)
    four <- densityplot(~gamma, data=altGamma, groups=index)
    if ( !is.null(x$null$c) )
        one <- densityplot(~c, data=nullC, groups=index)
    if ( !is.null(x$null$lambda) )
        one <- densityplot(~lambda, data=nullLambda, groups=index)

    if ( !is.null(x$alt$c) )
        three <- densityplot(~c, data=altC, groups=index)
    if ( !is.null(x$alt$lambda) )
        three <- densityplot(~lambda, data=altLambda, groups=index)

    if ( hypothesis == 'null' ) {
        grid.arrange(one, two, nrow=2, main = 'Null Hypothesis')
    } else if ( hypothesis == 'alt' ) {
        grid.arrange(three, four, nrow=2, main = 'Alternative Hypothesis')
    } else {
        grid.arrange(one, two, three, four, nrow=2, main = '')
    }
}

##' summarize output from testPref()
##' 
##' @param object a testPref object as returned by the eponymous function
##' @param ... additional arguments
##' @param hypothesis specify which hypothesis to plot
##' @S3method 
mean.testPref <- function(object, ..., hypothesis = c('null', 'alt', 'both')) {

    hypothesis <- match.arg(hypothesis)
    meanFn <- function(x, ...) mean(x, ...)
    
    ## null
    if ( !is.null(object$null$c) ) {
        null <- list('gamma' = apply(object$null$gamma, 2, meanFn), 'c' = apply(object$null$c, 2, meanFn))
    } else {
        null <- list('gamma' = apply(object$null$gamma, 2, meanFn), 'lambda' = apply(object$null$lambda, 2, meanFn))
    }

    ## alt
    if ( !is.null(object$alt$c) ) {
        alt <- list('gamma' = apply(object$alt$gamma, 2, meanFn), 'c' = apply(object$alt$c, 2, meanFn))
    } else {
        alt <- list('gamma' = apply(object$alt$gamma, 2, meanFn), 'lambda' = apply(object$alt$lambda, 2, meanFn))
    }

    if ( isTRUE(hypothesis == 'null') ) {
        null
    } else if ( isTRUE(hypothesis == 'alt') ) {
        alt
    } else {
        list('null' = null, 'alt' = alt)
    }
}
