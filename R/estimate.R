##' estimates parameters from hypothesis lambda = gamma; S*T free parameters
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
est1 <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    JI <- J + I
    if (EM) {
        em_iter <- 1
        init <- est1(Xdst, Ydst, J, I, FALSE, em_maxiter, BALANCED)
        gammaHat <- gammaHat_old <- init$gamma

        ## iterate EM
        while (EM) {
            
            ## some numbers
            lambda <- gammaHat
            elambda <- exp(lambda)
            EX <- lambda*elambda / (elambda - 1)
            ZEX <- Xdst*EX

            ## update parameters
            gammaHat <- (ZEX + Ydst) / JI

            ## check convergence of EM
            if ( converged(gammaHat, gammaHat_old) ) break
            
            ## if not converged, store updated estimates
            gammaHat_old <- gammaHat

            ## limit iterations
            em_iter <- em_iter+1
            if ( em_iter > em_maxiter )
                stop(sprintf('est1: max EM iterations, %d, reached. Please adjust accordingly.', em_maxiter))
        }

        ## calc standard error with est params
        egamma <- exp(-gammaHat)
        SE <- Ydst / gammaHat^2 + J*egamma/(egamma - 1)^2
        
        ## calc loglik with est params
        loglik <- llEM(Xdst, Ydst, gammaHat, gammaHat, J, I)
        
        list('lambda' = as.matrix(gammaHat), 'gamma' = as.matrix(gammaHat),
             'll' = loglik, 'se' = sqrt(1/SE))
        
    } else {
        
        XYdst <- Xdst + Ydst
        gammaHat <-  XYdst / JI

        ## calc standard error with est params
        SE <- gammaHat / sqrt(XYdst)

        ## calc loglik with est params
        loglik <- ll(Xdst, Ydst, gammaHat, gammaHat, J, I)

        list('lambda' = as.matrix(gammaHat), 'gamma' = as.matrix(gammaHat),
             'll' = loglik, 'se' = SE)        
    }
}

##' estimate parameters from hypothesis: lambda = c*gamma; S*T + 1 free parameters
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estC <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    if (EM) {
        ## initialize some values
        em_iter <- 1
        cHat <- cHat_old <- runif(1)        

        init <- est1(Xdst, Ydst, J, I, FALSE, em_maxiter, BALANCED)
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
        
        list('c' = cHat, 'gamma' = as.matrix(gammaHat), 'em_iters' = em_iter,
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
            
            list('gamma' = as.matrix(gammaHat), 'c' = cHat,
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
            
            list('gamma' = as.matrix(gammaHat), 'c' = cHat, 'iters' = iter,
                 'll' = loglik, 'se' = SE)
        }
    }
}

##' estimates parameters from hypothesis lambda = gamma c'_t; S*T + T free parameters
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estCt <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    if (EM) {
        ## initialize some values
        em_iter <- 1
        cHat <- cHat_old <- runif(length(J))
        init <- est1(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED)
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
        
        list('c' = cHat, 'gamma' = as.matrix(gammaHat), 'em_iters' = em_iter,
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
        
        list('gamma' = as.matrix(gammaHat), 'c' = cHat, 'iters' = iter,
             'll' = loglik, 'se' = SE)
    }
}

##' estimates parameters from hypothesis lambda = gamma c'_s; S*T + S free parameters
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estCs <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    ## some numbers
    S <- ncol(Xdst); s <- seq_len(S)
    T <- nrow(Xdst)

    if (EM) {
        ## initialize some values
        em_iter <- 1
        cHat <- cHat_old <- runif(S)
        init <- est1(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED)
        gammaHat <- gammaHat_old <- init$gamma 

        ## iterate EM
        while ( TRUE ) {
            
            ## expected value of Xjst
            lambda <- sapply(s, function(j) cHat[j]*gammaHat[,j]) # col-wise `*`
            elambda <- exp(lambda)
            EX <- lambda*elambda / (elambda - 1)

            ## convenience
            ZEX <- Xdst*EX

            ## iterate only once simultaneous eqs
            ZEXY <- ZEX + Ydst
            gammaHat <- sapply(s, function(j) ZEXY[,j] / (cHat[j]*J + I)) # col-wise `/`
            cHat <- sumT(ZEX) / sumT(J*gammaHat)
            
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
        
        list('c' = cHat, 'gamma' = as.matrix(gammaHat), 'em_iters' = em_iter,
             'll' = loglik, 'se' = SE)
    } else {

        XYdst <- Xdst + Ydst
        stXdst <- sumST(Xdst)
        iter <- 1; maxiter <- 500

        ## not sure this is the right spot for these checks
        ## ensure J & I have dimension T or 1
        if ( length(J) != T ) stop("J indexed oddly says est0.")
        if ( length(I) != T ) stop("I indexed oddly says est0.")

        ## initialize some values
        cHat <- cHat_old <- runif(S)
        gammaHat <- gammaHat_old <- sapply(s, function(j) XYdst[,j] / (J*cHat[j] + I))

        ## iteratively update; relies on concavity of log-lik
        while ( TRUE ) {

            ## update parameters
            cHat <- sumT(Xdst) / (sumT(J*gammaHat))
            gammaHat <- sapply(s, function(j) XYdst[,j] / (J*cHat[j] + I)) # col-wise `/`

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
            
        list('gamma' = as.matrix(gammaHat), 'c' = cHat, 'iters' = iter,
             'll' = loglik, 'se' = SE)
    }
}


##' estimates parameters from hypothesis lambda != gamma; 2*S*T free parameters
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period    
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estCst <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

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
        
        list('lambda' = as.matrix(lambdaHat), 'gamma' = as.matrix(gammaHat),
             'em_iters' = em_iter, 'll' = loglik, 'se' = SE)

    } else {
        ## some numbers
        ## not sure this is the right spot for these checks
        T <- nrow(Xdst)
        if ( length(J) != T ) stop("J indexed oddly says estGen")
        if ( length(I) != T ) stop("I indexed oddly says estGen")

        gammaHat <- Ydst / I                # row-wise division
        lambdaHat <- Xdst / J

        ## calc standard error with est params
        SE <- se(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)

        ## calc loglik with est params
        loglik <- ll(Xdst, Ydst, lambdaHat, gammaHat, J, I)
        
        list('lambda' = as.matrix(lambdaHat), 'gamma' = as.matrix(gammaHat),
             'll' = loglik, 'se' = SE)
    }

}
