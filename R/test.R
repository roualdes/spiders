##' test predPref function by simulating data and fitting model
##'
##' @param S number of prey species
##' @param T number of time periods
##' @param J number of predators caught at each time
##' @param I effective number of traps at each time
##' @param lambda matrix of rates at which predator eats prey species; TxS
##' @param gamma matrix of rates at which prey species is seen in habitat; TxS
##' @param M number of simulated datasets
##' @param EM boolean specifying test of EM algorithm
##' @param hyp a 2-tuple with names 'null' and 'alt' specifying the null and alternative hypotheses
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
##' @param n number of parameters to randomly sample; max allowed S*T
##' @export
testPref <- function(S, T, J, I, lambda, gamma, M=100, EM=F, hyp = c('c', 'gen'), em_maxiter = 100, n=4) {


    ## initialize output strucutres
    alt <- null <- as.data.frame(matrix(NA, nrow=M*n, ncol=3))
    colnames(alt) <- c('f', 'lambda', 'gamma')
    colnames(null) <- c('f', 'c', 'gamma')
    iters <- matrix(1, nrow=M, ncol=2)
    colnames(iters) <- c('null', 'gAtl')
    jdx <- 2:(S+1)                      # selects only data columns
    ST <- S*T
    if ( n <= ST ) {
        s <- sort(sample(1:ST, n))            # randomly sample n parameters
    } else {
        stop(sprintf('n must be smaller than S*T = %d', ST))
    }


    ## simulations
    for (m in seq_len(M)) {
        
        ## simulate data
        fdata <- simPref(S, T, J, I, lambda, gamma, EM=EM)
        
        ## fit model
        prefs <- predPref(fdata$eaten, fdata$caught, hypotheses = hyp, em_maxiter = em_maxiter)
        
        index_c <- ifelse(length(prefs$null$c)>1, TRUE, FALSE) # might need to fix this line in future
        idx <- (n*(m-1)+1); idxs <- idx:(idx+(n-1))
        
        ## store estimates
        if (EM) {
            iters[m,] <- c(prefs$null$em_iters, prefs$alt$em_iters)
        }
        alt[idxs,] <- cbind(s, unlist(prefs$alt$lambda)[s],
                            unlist(prefs$alt$gamma)[s])
        if (index_c) {
            r <- s[which(s<=T)]
            if (any(s>T)) r <- c(r, s[which(s>T)]-T)
            null[idxs,] <- cbind(s, prefs$null$c[r], unlist(prefs$null$gamma)[s])
        } else {
            null[idxs,] <- cbind(s, rep(prefs$null$c, n), unlist(prefs$null$gamma)[s])            
        }

    }
    out <- list('null' = null, 'alt' = alt, 'iters' = iters)
    class(out) <- 'testPref'
    out
}

##' plot the output of testPref
##'
##' @param x a testPref object as returned by the eponymous function
##' @export
plotTestPref <- function(x) {
    require(lattice)
    require(gridExtra)
    null <- x$null
    alt <- x$alt
    grid.arrange(densityplot(~c+gamma, data=null, groups=f,
                             main='Null Hypothesis', xlab='',
                             scales = list(y = list(relation = "free"),
                                 x = list(relation='free')),
                             ),
                 densityplot(~lambda+gamma, data=alt, group=f,
                             xlab='', main='Alternative Hypothesis',
                             scales = list(y = list(relation = "free"),
                                 x = list(relation='free'))),
                 nrow=2)
}

##' summarize output from testPref()
##' 
##' @param object a testPref object as returned by the eponymous function
##' @param ... additional arguments
##' @param fun a vectorized function to summarize each column by; defaults to mean()
##' @S3method 
summary.testPref <- function(object, ..., fun = mean) {
    fun <- match.fun(fun)
    list('null' = ddply(object$null, .(f), colwise(fun)),
         'alt' = ddply(object$alt, .(f), colwise(fun)))
}

##' calculate bias from output of testPref
##'
##' @param x a testPref object as returned by the eponymous function
##' @param lambda true values of lambda; TxS numbers in a data.frame
##' @param gamma true values of gamma; TxS numbers in a data.frame
##' @export
calcBias <- function(x, lambda, gamma) {
    
    ## ensure proper structure and calculate means of each randomly sampled element
    if ( class(x) != 'testPref' ) stop('argument out is not of correct class.')
    summOut <- summary(x)
    null <- summOut$null                # H0 means
    alt <- summOut$alt                  # H1 means
    f <- unique(null$f)                 # randomly sampled indices
    T <- nrow(lambda)                   # units of time
    c <- lambda/gamma                   # true value of c
        
    ## match randomly sampled indices with parameter estimates for c
    cidx <- f[which(f<=T)]
    cidx <- c(cidx,f[which(f>T)]-T)

    ## output parameter - mean(estimates)
    list('null' = data.frame('gamma' = gamma[f] - null$gamma,
             'c' = c[cidx] - null$c),
         'alt' = data.frame('lambda' = lambda[f] - alt$lambda,
             'gamma' = gamma[f] - alt$gamma))
}
