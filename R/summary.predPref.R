##' summary method for predPref objects
##'
##' @param object predPref object as returned from predPref()
##' @param sig.level significance level used in hypothesis test
##' @param digits number of significant digits to be printed
##' @param ... additional arguments affecting the call to format within summary fn
##' @export
summary.predPref <- function(object, sig.level = 0.05, digits=NULL, ...) {
    out <- vector('list', 2)
    names(out) <- c('statistics', 'p.value')
    neat <- format(object, digits=digits, ...)
    if ( object$p.value > sig.level ) {
        out[['statistics']] <- neat$null
    } else {
        out[['statistics']] <- neat$alt
    }
    out[['p.value']] <- object$p.value
    out
}
