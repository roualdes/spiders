##' estimates linear contrasts of c, c_s, c_t, c_st
##'
##' @param x a prefPref object as fit by the eponymous function
##' @param b a vector to linearly transform c_st
##' @param mu a number to test the linear contrast against in the null
##' @param alternative string to specify alternative hypothesis
##' @param conf.level confidence level of the interval
##' @export
testC <- function(x, b, mu = 0, alternative = c("two.sided", "less", "greater"), conf.level = 0.95) {
    
    ## some check on input; stolen from t.test
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")

    ## get appropriate estimates
    alpha <- 1-conf.level
    if ( x$LRT$p.value < alpha ) {
        if ( !is.null(x$alt$c) ) {
            C <- matrix(x$alt$c)
            lenC <- nrow(C); lc <- seq_len(lenC)
            varC <- x$alt$var[lc,lc]
        } else {
            stop("Don't yet know how to handle c_{st} = lambda_{st}/gamma_{st}.")
        }
    } else {
        if ( !is.null(x$null$c) ) {
            C <- matrix(x$null$c)
            lenC <- nrow(C); lc <- seq_len(lenC)
            varC <- x$null$var[lc,lc]
        } else {
            stop("Don't yet know how to handle 1 = lambda_{st}/gamma_{st}.")
        }        
    }

    ## checks on contrast
    if ( length(b) != lenC ) {
        stop("length of contrast b and length of C don't match.")
    }
        

    ## calculations
    stat <- t(b)%*%C
    scale <- t(b)%*%varC%*%b
    stdev <- sqrt(scale)

    ## test and confidence interval
    a2 <- alpha/2
    alternative <- match.arg(alternative)
    if ( alternative == 'two.sided' ) {
        pval <- 2*pnorm(-abs(stat), mean=mu, sd=stdev)
        interval <- qnorm(c(a2, 1-a2), mean=stat, sd=stdev)
    } else if ( alternative == 'less' ) {
        pval <- pnorm(stat, mean=mu, sd=stdev)
        interval <- c(-Inf, qnorm(conf.level, mean=stat, sd=stdev))
    } else if ( alternative == 'greater' ) {
        pval <- pnorm(stat, mean=mu, sd=stdev, lower.tail=F)
        interval <- c(qnorm(alpha, mean=stat, sd=stdev), Inf)
    }

    ## info to fill in print method of class htest
    names(stat) <- 'b^t * C'
    names(mu) <- 'linear contrast'
    method <- paste('Linear Contrast: ', paste(t(b), collapse=' '), sep='')
    attr(interval, 'conf.level') <- conf.level
    est <- as.vector(C)
    names(est) <- paste('mean of c_', lc, sep='')

    out <- list(statistic = stat, p.value = pval,
                conf.int = interval, null.value = mu,
                alternative = alternative, method = method,
                estimate = est, data.name=x$data.name)
    class(out) <- 'htest'
    out
}
