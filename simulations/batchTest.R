## libraries
library(BatchExperiments)
library(plyr)

## make / load registry
file.dir <- file.path('job1')
src <- file.path('utils.R')
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

## design parameters; experiments
hyps <- c('t', 's', 'c')
sz <- c('small', 'medium', 'large', 'huge')
em <- c(TRUE, FALSE)
for (h in seq_along(hyps)) {
    ## ensure only model hyps[h] is fit when it is the correct model to fit
    addExperiments(reg,
                   makeDesign('genPrefs',
                              exhaustive=list(letter=hyps[h], size=sz, EM=em)),
                   makeDesign('estPrefs',
                              exhaustive=list(letter=hyps[h])),
                   repls=500)
}
### remove experiments: removeExperiments(reg, findExperiments(reg))

summarizeExperiments(reg, show=c('prob', 'algo', 'letter', 'size', 'EM'))

## run
chunked <- chunk(getJobIds(reg), chunk.size=100, shuffle=TRUE)
submitJobs(reg, chunked)
## waitForJobs(reg)
## showStatus(reg)


