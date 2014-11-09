getSampleSize <- function(size) {
    Times <- 5
    sample_sizes <- list('small' = 20:50, 'medium' = 30:75, 'large' = 50:150, 'huge' = 100:200)
    sample(sample_sizes[[size]], Times)
}

getHypothesis <- function(letter) {
    hypotheses <- list('t' = c('Ct', 'Cst'), 's' = c('Cs', 'Cst'), 'c' = c('C', 'Ct'))
    hypotheses[[letter]]

}

getLambda <- function(letter) {
    PreySpecies <- 2
    Times <- 5
    lambda <- list('t' = matrix(1:5, nrow=Times, ncol=PreySpecies),
                   's' = matrix(c(rep(sqrt(2), Times), rep(pi, Times)),
                       nrow=Times, ncol=PreySpecies),
                   'c' = matrix(2*pi, nrow=Times, ncol=PreySpecies))
    lambda[[letter]]
}

getGamma <- function() {
    PreySpecies <- 2
    Times <- 5
    matrix(pi, nrow=Times, ncol=PreySpecies)
}






