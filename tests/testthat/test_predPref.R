context('Output predPref.')

## set up some data
Predators <- 5*c(11,22,33,44,77)
Traps <- Predators
PreySpecies <- 2
Times <- 5
TxS <- c(Times, PreySpecies)
g <- matrix(sqrt(2), nrow=Times, ncol=PreySpecies)


## noEM
test_that('noEM, c vs c^t, null', {
    l <- matrix(2*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=F)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'index_c'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(length(pref$alt$c), Times)
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_more_than(pref$p.value, 0.1)     # correct conclusion
})

test_that('noEM, c vs general, null', {
    l <- matrix(2*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=F)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'general'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$lambda), TxS)
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_more_than(pref$p.value, 0.1)     # correct conclusion    
})

test_that('noEM, c vs c^t, alt', {
    l <- matrix((1:5)*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=F)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'index_c'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(length(pref$alt$c), Times)
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_less_than(pref$p.value, 0.1)     # correct conclusion    
})

test_that('noEM, c vs c^t, alt', {
    l <- matrix(exp(seq(0.5, 1.5))*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=F)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'general'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$lambda), TxS)
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_less_than(pref$p.value, 0.1)     # correct conclusion    
})



## EM
test_that('EM, c vs c^t, null', {
    l <- matrix(0.5*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=T)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'index_c'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(length(pref$alt$c), Times)
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_more_than(pref$p.value, 0.1)     # correct conclusion
})

test_that('EM, c vs general, null', {
    l <- matrix(0.5*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=T)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'general'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$lambda), TxS)
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_more_than(pref$p.value, 0.1)     # correct conclusion    
})

test_that('EM, c vs c^t, alt', {
    l <- matrix(seq(0.5,1.5,length.out=5)*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=T)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'index_c'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(length(pref$alt$c), Times)
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_less_than(pref$p.value, 0.1)     # correct conclusion    
})

test_that('EM, c vs c^t, alt', {
    l <- matrix(exp(seq(0.5,1,length.out=5))*sqrt(2), nrow=Times, ncol=PreySpecies)
    fdata <- simPref(PreySpecies, Times, Predators, Traps, l, g, EM=T)
    pref <- predPref(fdata$eaten, fdata$caught, hypotheses=c('c', 'general'), em_maxiter = 50000)
    expect_is(pref, 'predPref')             # inherits class
    expect_equal(length(pref$null$c), 1)    # hypotheses correctly chosen
    expect_equal(dim(pref$null$gamma), TxS) # dimensions of estimates
    expect_equal(dim(pref$alt$gamma), TxS)
    expect_less_than(pref$p.value, 0.1)     # correct conclusion    
})


