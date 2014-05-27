##' sum specified columns by time
##'
##' @param data a data.frame
##' @param vars column variables in data to sum over
##' @param by extra variables to sum by
getMonthCounts <- function(data, vars, by) {
    if (missing(by)) by <- 'time' else by <- c('time', by)
    ddply(data, by, function(dfr, idx) colSums(as.matrix(dfr[,idx])), vars)
}

##' count number of spiders or traps in each unit of time
##'
##' @param data a data.table
getUnitCounts <- function(data) {
    ddply(data, .(time), summarize, total=length(time))$total
}


## sum over species to get a vector of values for each time period
sumSp <- function(mat) {
    matrix(rowSums(mat), nrow=nrow(mat))
}

## sum over times to get a vector of values for each species
sumT <- function(mat) {
    matrix(colSums(mat), ncol=ncol(mat))
}

sumST <- function(mat) {
    sum(mat)
}

## colors stolen from http://geography.uoregon.edu/datagraphics/color_scales.htm
cols <- c('orange1' = "#FFBF80", 'orange2' = "#FF8000",
          'yellow1' = "#FFFF99", 'yellow2' = "#FFFF33",
          'green1' = "#B2FF8C", 'green2' = "#33FF00",
          'blue1' = "#A6EDFF", 'blue2' = "#1AB2FF",
          'purple1' = "#CCBFFF", 'puple2' = "#664CFF",
          'red1' = "#FF99BF", 'red2' = "#E61A33")
