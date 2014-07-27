##' format method for predPref objects
##'
##' @param x a predPref object, as output by function predPref
##' @param ... additional arguments affecting the formatting
##' @export
format.predPref <- function(x, ...) {

    ## some numbers
    S <- ncol(x$null$gamma)
    s <- seq_len(S)
    T <- nrow(x$null$gamma)
    t <- seq_len(T)
    ST <- S*T
    indices <- unlist(lapply(s,
                             function(y) sapply(t,
                                                function(x) paste(x, y, sep=''))))
    params <- c('lambda', 'gamma')
    
    ## initialize output list
    out <- vector('list', 2)
    names(out) <- c('null', 'alternative')
    out[[1]] <- as.data.frame(matrix(NA, nrow=ST+1, 3))
    out[[2]] <- as.data.frame(matrix(NA, nrow=ST*2, 3))
    colnames(out[[1]]) <- colnames(out[[2]]) <- c('parameter', 'estimate', 'standard.error')
    gamma_names <- paste('gamma', indices, sep='')
    out[[1]][,'parameter'] <- c('c', gamma_names)
    out[[2]][,'parameter'] <- c(paste('lambda', indices, sep=''), gamma_names)

    
    ## fill null hypothesis summary table
    out[[1]][,'estimate'] <- c(x$null$c, unlist(x$null$gamma))
    out[[1]][,'standard.error'] <- x$null$SE

    ## fill alternative hypothesis summary table
    out[[2]][,'estimate'] <- c(unlist(x$alt$lambda), unlist(x$alt$gamma))
    out[[2]][,'standard.error'] <- as.vector(x$alt$SE)

    lapply(out, function(z) format(z, ...))
}
