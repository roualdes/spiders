##' simulate data for predator preferences model
##'
##' @param S number of prey species
##' @param T number of times
##' @param J number of predators caught at each time; assumed constant for now
##' @param I number of traps at each time; assumed constant for now
##' @param lambda vector of rates at which predator eats prey species, indexed by s and t
##' @param gamma vector of rates at which prey species is seen in habitat, indexed by s and t
##' @export
simData <- function(S, T, J, I, lambda, gamma) {
    ## some quick numbers
    ST <- S*T                           # rows of data frame
    ns <- seq_len(S)                    # index prey species
    nt <- seq_len(T)                    # index times

    ## initialize data frame
    preyNames <- paste('preySpecies', ns, sep='')
    time <- paste('time', nt, sep='')
    df <- as.data.frame(matrix(NA, nrow=ST, ncol=4))
    colnames(df) <- c('prey', 'time', 'eaten', 'caught')

    ## fill data frame
    df[,1] <- rep(preyNames, each=T)
    df[,2] <- rep(time, length.out=ST)
    df[,3] <- rpois(ST, J*lambda)
    df[,4] <- rpois(ST, I*gamma)
    df
}
