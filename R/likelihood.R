##' log-likelihood of predator preferances model
##'
##' @param Xdst sum of number of eaten prey species s during occurrence t
##' @param Ydst sum of number of caught prey species s during occurrence t
##' @param lambda vector of parameters relative to eaten prey, indexed by s and t
##' @param gamma vector of parameters relative to caught prey, indexed by s and t
##' @param J number of predators caught at each time; assumed constant for now
##' @param I number of traps at each time; assumed constant for now
##' @param c scalar in our base null hypothesis
ll <- function(Xdst, Ydst, lambda, gamma, J, I, c=NULL) {
    ## sums over s and t
    sumGamma <- sum(gamma)
    lnGamma <- log(gamma)

    ## log likelihood
    if (missing(c)) {
        ll <- -J*sum(lambda) + sum(Xdst*log(lambda)) 
    } else {
        ll <- -J*c*sumGamma + sum(Xdst*(lnGamma+log(c)))
    }
    ll - I*sumGamma + sum(Ydst*lnGamma)
}
