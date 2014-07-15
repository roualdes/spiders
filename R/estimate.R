## balanced, completely observed data

##' estimates parameters from the null model when data balanced
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est0b <- function(Xdst, Ydst, J, I) {
    i <- I[1]
    iXdY <- i*sumST(Xdst)/sumST(Ydst)
    gammaHat <- (Xdst + Ydst) / (iXdY + i)
    cHat <- iXdY/J[1]
    return(list('gamma' = gammaHat, 'c' = cHat))
}

##' estimates parameters from the general alternative model when data balanced
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period 
est1b <- function(Xdst, Ydst, J, I) {
    gammaHat <- Ydst/I
    lambdaHat <- Xdst/J
    return(list('gamma' = gammaHat, 'lambda' = lambdaHat))
}


## unbalanced, completely observed data

##' estimates parameters from the null model when data unbalanced
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est0 <- function(Xdst, Ydst, J, I) {
    eps <- 1e-4
    
    ## some numbers
    S <- ncol(Xdst)
    T <- nrow(Xdst)
    XYdst <- Xdst + Ydst
    stXdst <- sumST(Xdst)
    iter <- 1; maxiter <- 20

    ## not sure this is the right spot for these checks
    ## ensure J & I have dimension T or 1
    lJ <- length(J)
    lI <- length(I)
    if ( lJ != T || lJ != 1 ) stop("J indexed oddly says est0.")
    if ( lI != T || lI != 1 ) stop("I indexed oddly says est0.")

    ## initial estimates; only need c
    cHat <- cHat_old <- 1

    ## iteratively update; relies on concavity of log-lik
    while ( TRUE ) {

        ## update parameters
        gammaHat <- XYdst / (J*cHat + I) # row-wise division
        cHat <- stXdst / sumT(J*sumSp(gammaHat))

        ## check convergence
        if ( converged(gammaHat, gammaHat_old) &&
            converged(cHat, cHat_old) ) break

        ## if not converged, update estimates for next iteration
        gammaHat_old <- gammaHat
        cHat_old <- cHat
        
        if ( iter > maxiter ) break
        iter <- iter+1
    }
    return( list('gamma' = gammaHat, 'c' = cHat) )
}

##' estimates parameters from the general alternative model when data unbalanced
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est1 <- function(Xdst, Ydst, J, I) {

    ## some numbers
    T <- nrow(Xdst)
    if ( length(J) != T ) stop("J indexed oddly says est1")
    if ( length(I) != T ) stop("I indexed oddly says est1")

    gammaHat <- Ydst / I                # row-wise division
    lambdaHat <- Xdst / J
    return( list('lambda' = lambdaHat, 'gamma' = gammaHat) )
}

## non-count data

##' estimates parameters from the null model with non-count data via EM
##'
##' @param Zdst matrix of sums of binary response, whether prey species s was eaten or not during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
estEM0 <- function(Zdst, Ydst, J, I, em_maxiter){

    ## some numbers
    S <- ncol(Zdst)
    T <- nrow(Zdst)
    em_iter <- inner_iter <- 1;
    inner_maxiter <- 20

    ## initialize some values
    cHat <- cHat_old <- runif(1)
    init <- est1(Zdst, Ydst, J, I)
    gammaHat <- gammaHat_old <- init$gamma 
    lambda <- elambda <- init$lambda       

    ## iterate EM
    while ( TRUE ) {
        
        ## expected value of Xjst
        lambda <- cHat*gammaHat
        elambda <- exp(lambda)
        EX <- lambda*elambda / (elambda - 1)

        ## iterate system of eqs
        ## while ( TRUE ) {

            ## previous estimates
            ## cHat_old2 <- cHat
            ## gammaHat_old2 <- gammaHat
            
            ## convenience
            ZEX <- Zdst*EX

            ## iterate simultaneous eqs
            gammaHat <- (ZEX + Ydst) / (cHat*J + I)
            cHat <- sumST(ZEX) / sumT(J*sumSp(gammaHat))

            ## ## check convergence of system of eqs
            ## if ( converged(cHat, cHat_old2) &&
            ##     converged(gammaHat, gammaHat_old2) ) break            
            ## inner_iter <- inner_iter+1
            ## #print(sprintf('inner iteration %d', inner_iter))
            
        ##     ## limit iterations
        ##     if ( inner_iter > inner_maxiter ) break
        ## }

        ## check convergence of EM
        if ( converged(cHat, cHat_old) &&
            converged(gammaHat, gammaHat_old) ) break
        
        ## if not converged, store updated estimates
        cHat_old <- cHat
        gammaHat_old <- gammaHat
        em_iter <- em_iter+1
        # print(sprintf('EM iteration %d found values: c = %f', em_iter, cHat))
        
        ## limit iterations
        if ( em_iter > em_maxiter )
            stop(sprintf('max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
        inner_iter <- 1
    }
    return( list('c' = cHat, 'gamma' = gammaHat, 'em_iters' = em_iter) )
}

##' estimates parameters from the general alternative model with non-count data via EM
##'
##' @param Zdst matrix of sums of binary response, whether prey species s was eaten or not during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
estEM1 <- function(Zdst, Ydst, J, I, em_maxiter){

    ## some numbers
    S <- ncol(Zdst)
    T <- nrow(Zdst)
    em_iter <- 1
    
    ## estimate gamma
    gammaHat <- Ydst/I

    ## initialize lambdaHat
    init <- est1(Zdst, Ydst, J, I)
    lambdaHat <- lambdaHat_old <- init$lambda

    ## iterate EM
    while ( TRUE ) {
        
        ## update lambdas that aren't estimated 0
        elambda <- exp(lambdaHat)
        EX <- lambdaHat*elambda/(elambda-1)
        lambdaHat <- Zdst*EX/J

        ## check convergence
        if ( converged(lambdaHat, lambdaHat_old) ) break
        lambdaHat_old <- lambdaHat
        em_iter <- em_iter+1
        
        ## limit iterations
        if ( em_iter > em_maxiter )
            stop(sprintf('max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
    }
    
    ## since estimated 0s don't follow from calculations above
    lambdaHat[which(Zdst == 0, arr.ind=T)] <- 0
    
    return( list('lambda' = lambdaHat, 'gamma' = gammaHat, 'em_iters' = em_iter) )
}

