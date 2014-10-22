##' @title simulate and test predPref
##'
##' @description simulate data and test function \code{predPref} on each simulated dataset
##'
##' @return A list of two elements, one for each hypothesis. Each element contains
##' parameter estimates as calculated by the respective hypothesis, and the number of
##' iterations used to fit said model.
##'
##' @param J number of predators caught at each time
##' @param I effective number of traps at each time
##' @param lambda matrix of rates at which predator eats prey species; TxS
##' @param gamma matrix of rates at which prey species is seen in habitat; TxS
##' @param M number of simulated datasets
##' @param hyp a 2-tuple specifying the null and alternative hypotheses, respectively
##' @param EM boolean specifying test of EM algorithm
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##' @param storeSeed boolean to store random seed(s) or not
##'
##' @seealso \code{\link{simPref}} \code{\link{plotTestPref}} \code{\link{mean.testPref}}
##'
##' @examples
##'
##' # set parameters
##' Predators <- Traps <- 100
##' PreySpecies <- 2
##' Times <- 5
##' g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)     # gamma
##' l <- matrix(seq(0.4,1.8,length.out=5)*sqrt(2), nrow=Times, ncol=PreySpecies) # ct
##'
##' # test functions
##' \dontrun{
##' testPref(Predators, Traps, l, g, M=10, hyp=c('c', 'cst'), EM=F)
##' }
##'
##' @export
testPref <- function(J, I, lambda, gamma, M=100, hyp=c('C', 'Cst'), EM = FALSE, em_maxiter = 100, storeSeed = FALSE) {

    ## initialize output structure
    out <- vector('list', 2)
    names(out) <- c('null', 'alt')
    needSetOutput <- TRUE             # need initialize storage within output?
    if ( storeSeed ) {
        rseed <- vector('list', M)
        invisible(runif(1))             # ensure random seed exists
    }

    ## some numbers 
    S <- ncol(lambda)
    T <- nrow(lambda)
    
    for ( m in seq_len(M) ) {

        if ( storeSeed ) {
            rseed[[m]] <- .Random.seed
        }
        
        ## simulate data and fit model
        fdata <- simPref(S, T, J, I, lambda, gamma, EM=EM)
        prefs <- predPref(fdata$eaten, fdata$caught, hypotheses = hyp, em_maxiter = em_maxiter)

        ## initialize storage within output structure
        if ( needSetOutput ) {
            ST <- length(unlist(prefs$null$gamma))
            out$null[['gamma']] <- out$alt[['gamma']] <- matrix(0, M, ST)
            
            ## null model storage
            if ( !is.null(prefs$null$c) ) {
                out$null[['c']] <- matrix(0, M, length(prefs$null$c))
            }
            if ( !is.null(prefs$null$lambda) ) {
                out$null[['lambda']] <- matrix(0, M, ST)
            }
            if ( !is.null(prefs$null$em_iters) ) {
                out$null[['em_iters']] <- rep(0, M)
            }
            if ( !is.null(prefs$null$iters) ) {
                out$null[['iters']] <- rep(0, M)
            }
                

            ## alt model storage
            if ( !is.null(prefs$alt$c) ) {
                out$alt[['c']] <- matrix(0, M, length(prefs$alt$c))
            }
            if ( !is.null(prefs$alt$lambda) ) {
                out$alt[['lambda']] <- matrix(0, M, ST)            
            }
            if ( !is.null(prefs$alt$em_iters) ) {
                out$alt[['em_iters']] <- rep(0, M)
            }
            if ( !is.null(prefs$alt$iters) ) {
                out$alt[['iters']] <- rep(0, M)
            }
            
            needSetOutput <- FALSE
        }

        ## store output
        out$null$gamma[m,] <- unlist(prefs$null$gamma)
        out$alt$gamma[m,] <- unlist(prefs$alt$gamma)

        ## null
        if ( !is.null(prefs$null$c) ) {
            out$null$c[m,] <- prefs$null$c
         }
        if ( !is.null(prefs$null$lambda) ) {
            out$null$lambda[m,] <- unlist(prefs$null$lambda)
         }
        if ( !is.null(prefs$null$em_iters) ) {
            out$null$em_iters <- prefs$null$em_iters
        }
        if ( !is.null(prefs$null$iters) ) {
            out$null[['iters']] <- prefs$null$iters
        }
        
        ## alt
        if ( !is.null(prefs$alt$c) ) {
            out$alt$c[m,] <- prefs$alt$c
        }
        if ( !is.null(prefs$alt$lambda) ) {
            out$alt$lambda[m,] <- unlist(prefs$alt$lambda)
        }
        if ( !is.null(prefs$alt$iters) ) {
            out$alt$iters <- prefs$alt$iters
        }
            
    }
    class(out) <- 'testPref'
    attr(out, 'ST') <- S*T
    if ( storeSeed ) {
        attr(out, "seed") <- rseed        
    }
    
    out
}

##' @title plot testPref
##'
##' @description plot the output from \code{testPref}
##'
##' @details Function relies on packages \code{lattice} and \code{gridExtra}.  If
##' true values of lambda and gamma are given, then tick marks will be placed
##' at those values along the x-axis.
##'
##' @param x a testPref object as returned by the eponymous function
##' @param hypothesis specify which hypothesis to plot
##' @param lambda a matrix of true values of the parameter lambda; TxS
##' @param gamma a matrix of true values of the parameter gamma; TxS
##'
##' @seealso \code{\link{simPref}} \code{\link{testPref}} \code{\link{mean.testPref}}
##' 
##' @export
plotTestPref <- function(x, hypothesis = c('null', 'alt'), lambda = NULL, gamma = NULL) {
    
    ## some numbers
    hypothesis <- match.arg(hypothesis)
    M <- nrow(x$null$gamma)                 # number of replications run
    ST <- attr(x, 'ST'); st <- seq_len(ST)
    params <- ifelse( !is.null(lambda) && !is.null(gamma), TRUE, FALSE )

    ## data
    nullGamma <- data.frame('gamma' = as.vector(x$null$gamma), 'index' = st)
    altGamma <- data.frame('gamma' = as.vector(x$alt$gamma), 'index' = st)

    ## more null data
    if ( !is.null(x$null$c) ) {
        nullC <- data.frame('c' = as.vector(x$null$c),
                            'index' = rep(seq_len(ncol(x$null$c)), each=M))        
    }
    if ( !is.null(x$null$lambda) ) {
        nullLambda <- data.frame('lambda' = as.vector(x$null$lambda), 'index' = st)
    }
        

    ## more alt data
    if ( !is.null(x$alt$c) ) {
        altC <- data.frame('c' = as.vector(x$alt$c),
                           'index' = rep(seq_len(ncol(x$alt$c)), each=M))
    }
    if ( !is.null(x$alt$lambda) ) {
        altLambda <- data.frame('lambda' = as.vector(x$alt$lambda), 'index' = st)        
    }

    ## plot data
    require(lattice)

    if ( params ) {
        two <- densityplot(~gamma, data=nullGamma, groups=index,
                           panel = function(...) {
                               panel.densityplot(...)
                               panel.rug(x = unique(as.vector(gamma)), lwd=3)})
    } else {
        two <- densityplot(~gamma, data=nullGamma, groups=index)
    }

    if ( params ) {
        four <- densityplot(~gamma, data=altGamma, groups=index,
                           panel = function(...) {
                               panel.densityplot(...)
                               panel.rug(x = unique(as.vector(gamma)), lwd=3)})
    } else {
        four <- densityplot(~gamma, data=altGamma, groups=index)
    }

    if ( !is.null(x$null$c) ) {
        if ( params ) {
            trueC <- unique(apply(lambda/gamma, 1, unique))
            one <- densityplot(~c, data=nullC, groups=index,
                               main='Null Hypothesis',
                               panel = function(...) {
                                   panel.densityplot(...)
                                   panel.rug(x = trueC, lwd=3)})
        } else {
            one <- densityplot(~c, data=nullC, groups=index,
                               main='Null Hypothesis')
        }
    }
        
    if ( !is.null(x$null$lambda) )
        one <- densityplot(~lambda, data=nullLambda, groups=index)

    if ( !is.null(x$alt$c) ) {
        if ( params ) {
            trueC <- unique(apply(lambda/gamma, 1, unique))
            three <- densityplot(~c, data=altC, groups=index,
                                 main='Alternative Hypothesis',
                                 panel = function(...) {
                                     panel.densityplot(...)
                                     panel.rug(x = trueC, lwd=3)})
        } else {
            three <- densityplot(~c, data=altC, groups=index,
                                 main='Alternative Hypothesis')            
        }
    }

    if ( !is.null(x$alt$lambda) )
        three <- densityplot(~lambda, data=altLambda, groups=index)

    if ( hypothesis == 'null' ) {
        print(one, split=c(1,1,1,2), more=TRUE)
        print(two, split=c(1,2,1,2))
    } else {
        print(three, split=c(1,1,1,2), more=TRUE)
        print(four, split=c(1,2,1,2))
    }
}

##' @title means of \code{testPref}
##'
##' @description summarize output of \code{testPref} by calculating means, across all
##' simulations, for a specified hypothesis
##'
##' @details This function essentially calculates means over multiple replications and
##' fits of a given set of hypotheses and parameter values, as is output from the
##' function \code{testPref}.
##' 
##' @param x a testPref object as returned by the eponymous function
##' @param ... additional arguments
##' @param hypothesis specify which hypothesis
##'
##' @seealso \code{\link{simPref}} \code{\link{testPref}} \code{\link{plotTestPref}}
##' 
##' @S3method 
mean.testPref <- function(x, ..., hypothesis = c('null', 'alt')) {

    hypothesis <- match.arg(hypothesis)
       
    ## null
    if ( !is.null(x$null$c) ) {
        null <- list(gamma = colMeans(x$null$gamma),
                     c = colMeans(x$null$c))
    } else {
        null <- list(gamma = colMeans(x$null$gamma),
                     lambda = colMeans(x$null$lambda))
    }

    ## alt
    if ( !is.null(x$alt$c) ) {
        alt <- list(gamma = colMeans(x$alt$gamma),
                    c = colMeans(x$alt$c))
    } else {
        alt <- list(gamma = colMeans(x$alt$gamma),
                    lambda = colMeans(x$alt$lambda))
    }

    if ( isTRUE(hypothesis == 'null') ) {
        null
    } else {
        alt
    }
}
