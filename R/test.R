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
    ST <- S*T
    if ( n <= ST ) {
        s <- sample(1:ST, n)            # randomly sample n parameters
    } else {
        stop(sprintf('n must be smaller than S*T = %d', ST))
    }


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
##' @param null null hypothesis dataset
##' @param alt alternative hypothesis dataset
plotTestPref <- function(null, alt) {
    require(lattice)
    require(gridExtra)
    grid.arrange(densityplot(~c+gamma, data=null, groups=f,
                             main='Null Hypothesis', xlab='',
                             scales = list(y = list(relation = "free"),
                                 x = list(relation='free'))),
                 densityplot(~lambda+gamma, data=alt, group=f,
                             xlab='', main='Alternative Hypothesis',
                             scales = list(y = list(relation = "free"),
                                 x = list(relation='free'))),
                 nrow=2)
}

##' calculate bias from output of testPref
##'
##' @param null null hypothesis dataset
##' @param alt alternative hypothesis dataset
##' @param lambda true values of lambda; TxS numbers
##' @param gamma true values of gamma; TxS numbers
calcBias <- function(null, alt, lambda, gamma) {
    ## get bias = param - est for each simulation
    s <- unique(null$f)
    H0bias <- as.data.frame(sapply(s, function(x) gamma[x] - null$gamma[which(null$f == x)]))
    H1biasG <- as.data.frame(sapply(s, function(x) gamma[x] - alt$gamma[which(alt$f == x)]))
    H1biasL <- as.data.frame(sapply(s, function(x) lambda[x] - alt$lambda[which(alt$f == x)]))
    colnames(H0bias) <- colnames(H1biasG) <- colnames(H1biasL) <- as.character(s)

    ## format alternative hypothesis data
    sL <- stack(H1biasL); sG <- stack(H1biasG)
    oL <- order(sL$ind); oG <- order(sG$ind)
    altDat <- as.data.frame(cbind(sL[oL,'values'], sG[oG,]))
    colnames(altDat) <- c('lambda', 'gamma', 'ind')

    ## format null data
    nullDat <- stack(H0bias)
    nullDat$variable <- 'gamma'
    
    return( list('null' = nullDat,
                 'alt' = melt(altDat, id.vars='ind')) )
}

##' plot bias calculations from calcBias
##'
##' @param Bias output from calcBias
##' @details Means and medians are represented by dots and bars, respectively.  
plotBias <- function(Bias) {
    require(lattice)
    require(gridExtra)
    nullM <- ddply(Bias$null, .(ind, variable), summarize, means=mean(values))
    altM <- ddply(Bias$alt, .(ind, variable), summarize, means=mean(value))
    grid.arrange(bwplot(values~ind|variable, data=Bias$null, main='Null Hypothesis',
                        strip=strip.custom(var.name='gamma'), pch='|', ylab='Bias',
                        xlab='',
                        panel = function(...) {
                            panel.bwplot(...)
                            panel.points(x=nullM$means, col='black', pch=20)
                            ## ensures order of means is correct
                            #panel.text(x=nullM$ind, y=0.5, labels=as.character(nullM$ind))
                            panel.abline(h=0, col='black', lty='dotted')
                        }),
                 bwplot(value~ind|variable, data=Bias$alt, main='Alternative Hypothesis',
                        pch='|', ylab='Bias', xlab='',
                        panel = function(...) {
                            panel.bwplot(...)
                            panel.points(x=altM$means, col='black', pch=20)
                            ## ensures order of means is correct
                            #panel.text(x=altM$ind, y=0.5, labels=as.character(altM$ind))
                            panel.abline(h=0, col='black', lty='dotted')
                        }),
                 nrow=2)
}
