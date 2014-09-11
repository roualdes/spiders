##' summary method for predPref objects
##'
##' @param object predPref object as returned from predPref()
##' @param ... additional arguments
##' @param sig.level significance level used in hypothesis test
##' @export
summary.predPref <- function(object, ..., sig.level=0.05) {
    out <- vector('list', 7)
    names(out) <- c('loglikH0', 'loglikH1', 'p.value',
                    'estimates', 'df', 'Lambda', 'hypotheses')
    if ( object$p.value > sig.level ) {
        out[['estimates']] <- object$null
    } else {
        out[['estimates']] <- object$alt
    }
    out[['p.value']] <- object$p.value
    out[['df']] <- object$df
    out[['Lambda']] <- object$Lambda
    out[['loglikH0']] <- object$loglikH0
    out[['loglikH1']] <- object$loglikH1
    out[['hypotheses']] <- object$hypotheses
    class(out) <- 'summary.predPref'
    out
}

print.summary.predPref <- function(x) {
    cat("\nPredator Preferences Model:\n\n")
    cat("Log-likelihoods:\n\tH0 (", x$hypotheses[1], "):", x$loglikH0, "\n\tH1 (", x$hypotheses[2], "):", x$loglikH1, "\n")
    cat("Likelihood Ratio Test:\n\t-2*(llH0-llH1) =", x$Lambda, "on", x$df, "degrees of freedom\n")
    cat("\tp-value =", x$p.value, "\n")
    cat("Parameter Estimates:\n")

    ## some numbers
    S <- ncol(x$estimates$gamma); s <- seq_len(S)
    T <- nrow(x$estimates$gamma); t <- seq_len(T)
    indices <- unlist(lapply(s,
                             function(y) sapply(t,
                                                function(x) paste(x, '_', y, sep=''))))
    if ( !is.null(x$estimates$c) ) {
        est <- cbind(c(as.vector(x$estimates$c),
                       as.vector(x$estimates$gamma)),
                     as.vector(sqrt(diag(x$estimates$var))))
        rn <- c(paste('c_', seq_along(x$estimates$c), sep=''),
                paste('gamma_', indices, sep=''))
    } else {
        if ( converged(x$estimates$lambda, x$estimates$gamma) ) {
            est <- cbind(as.vector(x$estimates$gamma),
                     as.vector(sqrt(diag(x$estimates$var))))
            rn <- paste('gamma_', indices, sep='')
        } else {
            est <- cbind(c(as.vector(x$estimates$lambda),
                           as.vector(x$estimates$gamma)),
                         as.vector(sqrt(diag(x$estimates$var))))
            rn <- c(paste('lambda_', indices, sep=''),
                    paste('gamma_', indices, sep=''))            
        }
    }

    colnames(est) <- c('estimate', 'Std. Err.')
    rownames(est) <- rn
    printCoefmat(est)
}
