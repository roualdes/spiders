##' estimates parameters of predator preferences model and calculates LRT
##'
##' @param dataframe a df of observed data
##' @param J number of predators caught at each time; assumed constant for now
##' @param I number of traps at each time; assumed constant for nwo
##' @param alpha LRT level of significance
##' @export
predpref <- function(dataframe, J, I, alpha=0.05) {
    Xdst <- dataframe[,'eaten']
    Ydst <- dataframe[,'caught']

    ## H1; 2ST parameters
    lambdaH1 <- Xdst/J
    gammaH1 <- Ydst/I
    llH1 <- ll(Xdst, Ydst, lambdaH1, gammaH1, J, I)

    ## H0; ST+1 parameters
    ratio <- sum(Xdst)/sum(Ydst)
    cH0 <- I*ratio/J
    gammaH0 <- (Xdst+Ydst)/(I*ratio + I)
    llH0 <- ll(Xdst, Ydst, NA, gammaH0, J, I, cH0)

    ## LRT
    df <- nrow(dataframe)-1
    LRT <- -2*(llH0 - llH1)
    
    list('lambdaH1' = lambdaH1, 'gammaH1' = gammaH1,
         'gammaH0' = gammaH0, 'cH0' = cH0,
         'loglikH1' = llH1, 'loglikH0' = llH0,
         '-2log(LRT)' = LRT, 'df' = df, 'p.value' = pchisq(LRT, df=df, lower.tail=F))
}
