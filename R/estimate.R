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
    list('gamma' = gammaHat, 'c' = cHat)
}

## unbalanced, completely observed data

##' estimates parameters from the null model when data unbalanced
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param indexC boolean; TRUE indexes c in null hypothesis by t
est0 <- function(Xdst, Ydst, J, I, indexC) {
    
    ## some numbers
    S <- ncol(Xdst)
    T <- nrow(Xdst)
    XYdst <- Xdst + Ydst
    stXdst <- sumST(Xdst)
    iter <- 1; maxiter <- 500

    ## not sure this is the right spot for these checks
    ## ensure J & I have dimension T or 1
    lJ <- length(J)
    lI <- length(I)
    if ( lJ != T ) stop("J indexed oddly says est0.")
    if ( lI != T ) stop("I indexed oddly says est0.")

    ## initialize some values
    if ( indexC ) {
        cHat <- cHat_old <- runif(length(J))
    } else {
        cHat <- cHat_old <- runif(1)
    }
    gammaHat <- gammaHat_old <- XYdst / (J*cHat + I)

    ## iteratively update; relies on concavity of log-lik
    while ( TRUE ) {

        ## update parameters
        if ( indexC ) {
            cHat <- sumSp(Xdst) / (J*sumSp(gammaHat))
        } else {
            cHat <- stXdst / sumT(J*sumSp(gammaHat))            
        }
        gammaHat <- XYdst / (J*cHat + I) # row-wise division

        ## check convergence
        if ( converged(gammaHat, gammaHat_old) &&
            converged(cHat, cHat_old) ) break

        ## if not converged, update estimates for next iteration
        gammaHat_old <- gammaHat
        cHat_old <- cHat
        iter <- iter+1
        
        if ( iter > maxiter )
            stop(sprintf('est0: %d not sufficient iterations for simultaneous equations.', maxiter))
    }
    list('gamma' = gammaHat, 'c' = cHat, 'iters' = iter)
}

##' estimates parameters from the general alternative model when data unbalanced
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est1 <- function(Xdst, Ydst, J, I) {

    ## some numbers
    ## not sure this is the right spot for these checks
    T <- nrow(Xdst)
    if ( length(J) != T ) stop("J indexed oddly says est1")
    if ( length(I) != T ) stop("I indexed oddly says est1")

    gammaHat <- Ydst / I                # row-wise division
    lambdaHat <- Xdst / J
    list('lambda' = lambdaHat, 'gamma' = gammaHat)
}

## non-count data

##' estimates parameters from the null model with non-count data via EM
##'
##' @param Zdst matrix of sums of binary response, whether prey species s was eaten or not during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param indexC boolean; TRUE indexes c in null hypothesis by t
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
estEM0 <- function(Zdst, Ydst, J, I, indexC, em_maxiter){

    ## initialize some values
    em_iter <- 1
    if ( indexC ) {
        cHat <- cHat_old <- runif(length(J))
    } else {
        cHat <- cHat_old <- runif(1)        
    }

    init <- est1(Zdst, Ydst, J, I)
    gammaHat <- gammaHat_old <- init$gamma 
    lambda <- elambda <- init$lambda

    ## iterate EM
    while ( TRUE ) {
        
        ## expected value of Xjst
        lambda <- cHat*gammaHat
        elambda <- exp(lambda)
        EX <- lambda*elambda / (elambda - 1)

        ## convenience
        ZEX <- Zdst*EX

        ## iterate only once simultaneous eqs
        gammaHat <- (ZEX + Ydst) / (cHat*J + I)
        if ( indexC ) {
            cHat <- sumSp(ZEX) / (J*sumSp(gammaHat))
        } else {
            cHat <- sumST(ZEX) / sumT(J*sumSp(gammaHat))
        }

        ## check convergence of EM
        if ( converged(cHat, cHat_old) &&
            converged(gammaHat, gammaHat_old) ) break
        
        ## if not converged, store updated estimates
        cHat_old <- cHat
        gammaHat_old <- gammaHat
        ## print(sprintf('EM iteration %d found values: c = %f', em_iter, cHat))
        em_iter <- em_iter+1
        
        ## limit iterations
        if ( em_iter > em_maxiter )
            stop(sprintf('H0: max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
    }
    list('c' = cHat, 'gamma' = gammaHat, 'em_iters' = em_iter)
}

##' estimates parameters from the general alternative model with non-count data via EM
##'
##' @param Zdst matrix of sums of binary response, whether prey species s was eaten or not during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param em_maxiter maximum number of iterations allowed for EM algorithm
estEM1 <- function(Zdst, Ydst, J, I, em_maxiter){

    ## estimate gamma
    gammaHat <- Ydst/I

    ## initialize some things
    em_iter <- 1
    lambdaHat <- lambdaHat_old <- gammaHat

    ## iterate EM
    while ( TRUE ) {
        ## update lambda
        elambda <- exp(lambdaHat)
        EX <- lambdaHat*elambda/(elambda-1)
        lambdaHat <- Zdst*EX/J

        ## check convergence
        if ( converged(lambdaHat, lambdaHat_old) ) break
        lambdaHat_old <- lambdaHat
        em_iter <- em_iter+1
        
        ## limit iterations
        if ( em_iter > em_maxiter )
            stop(sprintf('H1: max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
    }
    
    ## since estimated 0s don't follow from calculations above
    lambdaHat[which(Zdst == 0, arr.ind=T)] <- 0
    
    list('lambda' = lambdaHat, 'gamma' = gammaHat, 'em_iters' = em_iter)
}

## begin new functions ##


##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
est1t <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {
    
}

##' estimate parameters from hypothesis: lambda = c*gamma; S*T + 1 free parameters 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
estC <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    if (EM) {
        ## initialize some values
        em_iter <- 1
        cHat <- cHat_old <- runif(1)        

        init <- est1(Xdst, Ydst, J, I)
        gammaHat <- gammaHat_old <- init$gamma 
        lambda <- elambda <- init$lambda

        ## iterate EM
        while ( TRUE ) {
            
            ## expected value of Xjst
            lambda <- cHat*gammaHat
            elambda <- exp(lambda)
            EX <- lambda*elambda / (elambda - 1)

            ## convenience
            ZEX <- Xdst*EX

            ## iterate only once simultaneous eqs
            gammaHat <- (ZEX + Ydst) / (cHat*J + I)
            cHat <- sumST(ZEX) / sumT(J*sumSp(gammaHat))

            ## check convergence of EM
            if ( converged(cHat, cHat_old) &&
                converged(gammaHat, gammaHat_old) ) break
            
            ## if not converged, store updated estimates
            cHat_old <- cHat
            gammaHat_old <- gammaHat
            ## print(sprintf('EM iteration %d found values: c = %f', em_iter, cHat))
            em_iter <- em_iter+1
            
            ## limit iterations
            if ( em_iter > em_maxiter )
                stop(sprintf('H0: max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
        }
        
        ## calc standard error with est params
        SE <- seEM(NULL, gammaHat, cHat, Xdst, Ydst, J, I)

        ## calc log-lik with est params
        loglik <- llEM(Xdst, Ydst, NA, gammaHat, J, I, cHat)
        
        list('c' = cHat, 'gamma' = gammaHat, 'em_iters' = em_iter,
             'll' = loglik, 'se' = SE)
    } else {
        if (BALANCED) {
            i <- I[1]
            iXdY <- i*sumST(Xdst)/sumST(Ydst)
            gammaHat <- (Xdst + Ydst) / (iXdY + i)
            cHat <- iXdY/J[1]

            ## calc standard error with est params
            SE <- se(NULL, gammaHat, cHat, Xdst, Ydst, J, I)

            ## calc log-lik with est params
            loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cHat)
            
            list('gamma' = gammaHat, 'c' = cHat,
                 'll' = loglik, 'se' = SE)
        } else {
            ## some numbers
            S <- ncol(Xdst)
            T <- nrow(Xdst)
            XYdst <- Xdst + Ydst
            stXdst <- sumST(Xdst)
            iter <- 1; maxiter <- 500

            ## not sure this is the right spot for these checks
            ## ensure J & I have dimension T or 1
            if ( length(J) != T ) stop("J indexed oddly says est0.")
            if ( length(I) != T ) stop("I indexed oddly says est0.")

            ## initialize some values
            cHat <- cHat_old <- runif(1)
            gammaHat <- gammaHat_old <- XYdst / (J*cHat + I)

            ## iteratively update; relies on concavity of log-lik
            while ( TRUE ) {

                ## update parameters
                cHat <- stXdst / sumT(J*sumSp(gammaHat))            
                gammaHat <- XYdst / (J*cHat + I) # row-wise division

                ## check convergence
                if ( converged(gammaHat, gammaHat_old) &&
                    converged(cHat, cHat_old) ) break

                ## if not converged, update estimates for next iteration
                gammaHat_old <- gammaHat
                cHat_old <- cHat
                iter <- iter+1
                
                if ( iter > maxiter )
                    stop(sprintf('estC: %d not sufficient iterations for simultaneous equations.', maxiter))
            }
            
            ## calc standard error with est params
            SE <- se(NULL, gammaHat, cHat, Xdst, Ydst, J, I)

            ## calc log-lik with est params
            loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cHat)
            list('gamma' = gammaHat, 'c' = cHat, 'iters' = iter,
                 'll' = loglik, 'se' = SE)
        }
    }
}

##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
estCt <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    if (EM) {
        ## initialize some values
        em_iter <- 1
        cHat <- cHat_old <- runif(length(J))
        init <- est1(Xdst, Ydst, J, I)
        gammaHat <- gammaHat_old <- init$gamma 
        lambda <- elambda <- init$lambda

        ## iterate EM
        while ( TRUE ) {
            
            ## expected value of Xjst
            lambda <- cHat*gammaHat
            elambda <- exp(lambda)
            EX <- lambda*elambda / (elambda - 1)

            ## convenience
            ZEX <- Xdst*EX

            ## iterate only once simultaneous eqs
            gammaHat <- (ZEX + Ydst) / (cHat*J + I)
            cHat <- sumSp(ZEX) / (J*sumSp(gammaHat))
            
            ## check convergence of EM
            if ( converged(cHat, cHat_old) &&
                converged(gammaHat, gammaHat_old) ) break
            
            ## if not converged, store updated estimates
            cHat_old <- cHat
            gammaHat_old <- gammaHat
            ## print(sprintf('EM iteration %d found values: c = %f', em_iter, cHat))
            em_iter <- em_iter+1
            
            ## limit iterations
            if ( em_iter > em_maxiter )
                stop(sprintf('H0: max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
        }

        ## calc standard error with est params
        SE <- seEM(NULL, gammaHat, cHat, Xdst, Ydst, J, I)
        
        ## calc log-lik with est params
        loglik <- llEM(Xdst, Ydst, NA, gammaHat, J, I, cHat)
        
        list('c' = cHat, 'gamma' = gammaHat, 'em_iters' = em_iter,
             'll' = loglik, 'se' = SE)
    } else {
        ## some numbers
        S <- ncol(Xdst)
        T <- nrow(Xdst)
        XYdst <- Xdst + Ydst
        stXdst <- sumST(Xdst)
        iter <- 1; maxiter <- 500

        ## not sure this is the right spot for these checks
        ## ensure J & I have dimension T or 1
        if ( length(J) != T ) stop("J indexed oddly says est0.")
        if ( length(I) != T ) stop("I indexed oddly says est0.")

        ## initialize some values
        cHat <- cHat_old <- runif(length(J))
        gammaHat <- gammaHat_old <- XYdst / (J*cHat + I)

        ## iteratively update; relies on concavity of log-lik
        while ( TRUE ) {

            ## update parameters
            cHat <- sumSp(Xdst) / (J*sumSp(gammaHat))
            gammaHat <- XYdst / (J*cHat + I) # row-wise division

            ## check convergence
            if ( converged(gammaHat, gammaHat_old) &&
                converged(cHat, cHat_old) ) break

            ## if not converged, update estimates for next iteration
            gammaHat_old <- gammaHat
            cHat_old <- cHat
            iter <- iter+1
            
            if ( iter > maxiter )
                stop(sprintf('estCt: %d not sufficient iterations for simultaneous equations.', maxiter))
        }

        ## calc standard error with est params
        SE <- se(NULL, gammaHat, cHat, Xdst, Ydst, J, I)

        ## calc log-lik with est params
        loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cHat)
        list('gamma' = gammaHat, 'c' = cHat, 'iters' = iter,
             'll' = loglik, 'se' = SE)
    }
}

##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period    
estGen <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    if (EM) {
        
        ## estimate gamma
        gammaHat <- Ydst/I

        ## initialize some things
        em_iter <- 1
        lambdaHat <- lambdaHat_old <- gammaHat

        ## iterate EM
        while ( TRUE ) {
            ## update lambda
            elambda <- exp(lambdaHat)
            EX <- lambdaHat*elambda/(elambda-1)
            lambdaHat <- Xdst*EX/J

            ## check convergence
            if ( converged(lambdaHat, lambdaHat_old) ) break
            lambdaHat_old <- lambdaHat
            em_iter <- em_iter+1
            
            ## limit iterations
            if ( em_iter > em_maxiter )
                stop(sprintf('H1: max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
        }
        
        ## since estimated 0s don't follow from calculations above
        lambdaHat[which(Xdst == 0, arr.ind=T)] <- 0

        ## calc standard error with est params
        SE <- seEM(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)

        ## calc loglik wit est params
        loglik <- llEM(Xdst, Ydst, lambdaHat, gammaHat, J, I)                
        list('lambda' = lambdaHat, 'gamma' = gammaHat, 'em_iters' = em_iter,
             'll' = loglik, 'se' = SE)

    } else {
        ## some numbers
        ## not sure this is the right spot for these checks
        T <- nrow(Xdst)
        if ( length(J) != T ) stop("J indexed oddly says est1")
        if ( length(I) != T ) stop("I indexed oddly says est1")

        gammaHat <- Ydst / I                # row-wise division
        lambdaHat <- Xdst / J

        ## calc standard error with est params
        SE <- se(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)

        ## calc loglik with est params
        loglik <- ll(Xdst, Ydst, lambdaHat, gammaHat, J, I)
        
        list('lambda' = lambdaHat, 'gamma' = gammaHat,
             'll' = loglik, 'se' = SE)
    }

}
