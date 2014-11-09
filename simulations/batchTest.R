## libraries
library(BatchExperiments)
library(plyr)

## make / load registry
file.dir <- file.path('/Users/easy-e/spiders/simulations/job1')
src <- file.path('/Users/easy-e/spiders/simulations/utils.R')
reg <- makeExperimentRegistry("spidersRegistry", file.dir=file.dir, seed=123,
                              packages='spiders',
                              src.files=src)

### temporarily set parallel execution
setConfig(conf=list(cluster.functions=makeClusterFunctionsMulticore(ncpus=parallel::detectCores())))

## problems
genPrefs <- function(size, letter, EM) {
    ## get parameters
    J <- getSampleSize(size)            # eaten
    I <- getSampleSize(size)            # caught
    lambda <- getLambda(letter)         # eaten rate(s)
    gamma <- getGamma()                 # caught rate(s)

    ## some numbers
    S <- ncol(lambda)                   # number species
    T <- nrow(lambda)                   # number time points
    ## simulate data
    fdata <- simPref(S, T, J, I, lambda, gamma, EM=EM)

    ## return list of things we want later
    list(eaten = fdata$eaten, caught = fdata$caught, J=J, I=I)
}
addProblem(reg, 'genPrefs', dynamic=genPrefs, seed=125)
### remove problems: sapply(getProblemIds(reg), function(x) removeProblem(reg, x))

## algorithms
estPrefs <- function(dynamic, letter) {
    eat <- dynamic$eaten
    caught <- dynamic$caught
    hyp <- getHypothesis(letter)
    list(prefs=predPref(eat, caught, hypotheses=hyp, em_maxiter=15000), J=dynamic$J, I=dynamic$I)
}
addAlgorithm(reg, 'estPrefs', fun=estPrefs)
### remove algorithms: sapply(getAlgorithmIds(reg), function(x) removeAlgorithm(reg, x))

## design parameters
### problems
pr.params <- list(size=c('small', 'medium', 'large', 'huge'),
                  letter=c('t', 's', 'c'), EM=c(TRUE, FALSE))
pr.des <- makeDesign('genPrefs', exhaustive=pr.params)

## algorithms
alg.des <- makeDesign('estPrefs', exhaustive=list(letter=c('t', 's', 'c')))

## experiments
addExperiments(reg, pr.des, alg.des, repls=1)
### remove experiments: removeExperiments(reg, findExperiments(reg))

summarizeExperiments(reg)

## run
chunked <- chunk(getJobIds(reg), n.chunks=4, shuffle=TRUE)
submitJobs(reg, chunked)
waitForJobs(reg)
## showStatus(reg)

## summarize results
idx <- findExperiments(reg,
                       prob.pars=EM==TRUE && letter=='t',  # && size=='huge',
                       algo.pars=letter=='t')

summ <- function(job, res) {
    sumry <- summary(res$prefs)
    c(sumry$estimates$c)
}

reduceResultsMatrix(reg, ids=idx, fun=summ, rows=TRUE)

fdata <- generateProblemInstance(reg, 7)
pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('ct', 'cst'), em_maxiter=15000)
summary(pref)$estimates$c

