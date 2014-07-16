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
##' @param n number of parameters to randomly sample; max allowed S*T
testPref <- function(S, T, J, I, lambda, gamma, M=100, EM=F, n=4) {


    ## initialize output strucutres
    alt <- null <- as.data.frame(matrix(NA, nrow=M*n, ncol=3))
    colnames(alt) <- c('f', 'lambda', 'gamma')
    colnames(null) <- c('f', 'c', 'gamma')
    iters <- matrix(1, nrow=M, ncol=2)
    colnames(iters) <- c('null', 'gAtl')
    jdx <- 2:(S+1)                      # selects only data columns
    s <- sample(1:(S*T), n)             # randomly sample n parameters

    ## simulations
    for (m in seq_len(M)) {
        
        ## simulate data
        fdata <- simPref(S, T, J, I, lambda, gamma)
        if (EM) fdata$eaten[,jdx][which(fdata$eaten[,jdx]>0, arr.ind=T)] <- 1
        
        ## fit model
        prefs <- predPref(fdata$eaten, fdata$caught)
        idx <- (n*(m-1)+1); idxs <- idx:(idx+(n-1))
        
        ## store estimates
        if (EM) {
            iters[m,] <- c(prefs$null$em_iters, prefs$gAlt$em_iters)
        }
        alt[idxs,] <- cbind(s, unlist(prefs$gAlt$lambda)[s],
                            unlist(prefs$gAlt$gamma)[s])
        null[idxs,] <- cbind(s, rep(prefs$null$c, n), unlist(prefs$null$gamma)[s])
    }
    return( list('null' = null, 'alt' = alt, 'iters' = iters))
}

##' plot the output of testPref
##'
##' @param H0 null hypothesis dataset
##' @param H1 alternative hypothesis dataset
plotTestPref <- function(H0, H1) {
    require(lattice)
    require(gridExtra)
    grid.arrange(densityplot(~c+gamma, data=H0, groups=f,
                             main='Null Hypothesis', xlab='',
                             scales = list(y = list(relation = "free"),
                                 x = list(relation='free'))),
                 densityplot(~lambda+gamma, data=H1, group=f,
                             xlab='', main='Alternative Hypothesis',
                             scales = list(y = list(relation = "free"),
                                 x = list(relation='free'))),
                 nrow=2)
}
