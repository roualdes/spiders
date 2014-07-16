##' simulate data for predator preferences model
##'
##' @param S number of prey species
##' @param T number of time periods
##' @param J number of predators caught at each time
##' @param I effective number of traps at each time
##' @param lambda matrix of rates at which predator eats prey species; TxS
##' @param gamma matrix of rates at which prey species is seen in habitat; TxS
##' @export
simPref <- function(S, T, J, I, lambda, gamma) {
    
    ## some numbers
    ns <- seq_len(S)                    # index prey species
    nt <- seq_len(T)                    # index times

    ## initialize data frame
    eaten <- as.data.frame(matrix(NA, nrow=J*T, ncol=S+1))
    caught <- as.data.frame(matrix(NA, nrow=I*T, ncol=S+1))

    colnames(eaten) <- colnames(caught) <- c('time', paste('preySpecies', ns, sep=''))
    eaten[,1] <- rep(paste('time', nt, sep=''), each=J)
    caught[,1] <- rep(paste('time', nt, sep=''), each=I)

    ## fill data frame
    for ( i in ns ) {
        eaten[,i+1] <- rpois(J*T, lambda[,i])
        caught[,i+1] <- rpois(I*T, gamma[,i])            
    }
    list('eaten' = eaten,
         'caught' = caught)
}
