##' summary method for predPref objects
##'
##' @param object predPref object as returned from predPref()
##' @param ... additional arguments
##' @param sig.level significance level used in hypothesis test
##' @export
summary.predPref <- function(object, ..., sig.level=0.05) {
    out <- vector('list', 6)
    names(out) <- c('loglikH0', 'loglikH1', 'p.value', 'estimates', 'df', 'Lambda')
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
    class(out) <- 'summary.predPref'
    out
}

print.summary.predPref <- function(x) {
    cat("\nPredator Preferences Model:\n\n")
    cat("Log-likelihoods:\n\tH0:", x$loglikH0, "\n\tH1:", x$loglikH1, "\n")
    cat("Likelihood Ratio Test:\n\t-2*(llH0-llH1) =", x$Lambda, "on", x$df, "degrees of freedom\n")
    cat("\tp-value =", x$p.value, "\n")
    cat("Parameter Estimates:\n")
    print(x$estimates)
}
