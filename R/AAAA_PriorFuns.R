###################################################################
### This file describes several functions used in whole package ###
###################################################################
### functions:
##      paraMerge, merge two lists.
##      p_function, calculate p value.
##      Cpp_p_funciton, calculate p value.
##      plotMat, calcute the plot design.
##      .onUnload, clean up after unloading.

### paraMerge
#' Merge two lists
#'
#' \code{paraMerge} merges two arguments list into a new one, in which the old parameters
#' are over-written by new parameters.
#'
#' @param oldP A list containing the old parameters.
#' @param newP A list containing the new parameters.
#' @return A list of merged parameters is returned.
#' @examples
#' paraMerge(list(t1 = 1, t2 = 2, t3 = 3, t4 = 4), list(t3 = 4, t4 = 5, t5 = 6, t6 = 7))
#' @export
paraMerge = function(oldP, newP)
{
    if(any(is.null(newP), length(newP) == 0)) return(oldP)
    if(any(is.na(newP))) return(NA)
    if(length(newP) == 1 && is.na(newP)) return(list(NA))
    oldP = as.list(oldP)
    newP = as.list(newP)
    oldP = as.environment(oldP)
    lapply(names(newP), function(nm) oldP[[nm]] = newP[[nm]])

    return(as.list(oldP))
}

### p_function
#' Calculate the p-value
#'
#' \code{p_function} calculate the p-value of data given the epsilon and the test type.
#' \code{Cpp_p_function}, an alternative for the same purpose, is much faster and recommended.
#'
#' @param epsilon A value used for equivalence tests. Only valid when ifF == FALSE.
#' @param ifF A bool value indicating test type. When ifF = TRUE, standard F test is used; otherwise, equivalence F test.
#' @param myData A list of vector of numerics.
#' @return P value.
#' @examples
#' a = rnorm(100)
#' b = rnorm(150)
#' d = rnorm(200)
#' e = runif(250)
#' p_function(0, FALSE, list(a, b, d, e))
#' Cpp_p_function(0, FALSE, list(a, b, d, e))
#' @export
p_function = function(epsilon = 0, ifF = FALSE, myData = list())
{
    mylist = myData
    k = length(mylist)
    n = sum(sapply(mylist, length))
    xbar = sum(sapply(mylist, sum)) / n
    num = sum(sapply(mylist, function(x) length(x) / n * k * (mean(x) - xbar)^2))
    denum = 1 / (n - k) * sum(sapply(mylist, function(x) sum((x - mean(x))^2)))
    myfvalue = num / denum
    if(round(num, 18) == 0) myfvalue = 0

    fvalue = myfvalue * n / k / (k - 1)
    adjep = n / k * epsilon^2
    pvalue = (1 - ifF) * pf(fvalue, k - 1, n - k, adjep) + ifF * (1 - pf(fvalue, k - 1, n - k))
    return(pvalue)
}

### Cpp_p_function
#' Calculate the p-value
#'
#' \code{p_function} calculate the p-value of data given the epsilon and the test type.
#' \code{Cpp_p_function}, an alternative for the same purpose, is much faster and recommended.
#'
#' @param epsilon A value used for equivalence tests. Only valid when ifF == FALSE.
#' @param ifF A bool value indicating test type. When ifF = TRUE, standard F test is used; otherwise, equivalence F test.
#' @param myData A list of vector of numerics.
#' @return P value.
#' @examples
#' a = rnorm(100)
#' b = rnorm(150)
#' d = rnorm(200)
#' e = runif(250)
#' p_function(0, FALSE, list(a, b, d, e))
#' Cpp_p_function(0, FALSE, list(a, b, d, e))
#' @export
Cpp_p_function = function(epsilon = 0, ifF = FALSE, myData = list()) {
    res = Cpp_fvalue(myData)
    return(ifelse(ifF, (1 - pf(res[1], res[3] - 1, res[2] - res[3])), pf(res[1], res[3] - 1, res[2] - res[3], epsilon * epsilon * res[2] / res[3])))
}

### plotMat
#' Calcute the layout coordinates
#'
#' \code{plotMat} calcutes the layout coordinates given numbers of figures and columns of figures.
#' The order should be from top to bottom and left to right.
#'
#' @param nfig An integer indicating the total number of figures.
#' @param figcol An integer indicating the columns of figures.
#' @param xshrink,yshrink The right or top shrink.
#' @return A list of positions is returned with each element containing a numeric vector
#' of length 4, standing for (xleft, xright, ybottom, ytop)
#' @examples
#' plotMat(10, 3)
#' @export
plotMat = function(nfig, figcol, xshrink = 1, yshrink = 1)
{
    return(Cpp_plotMat(nfig = nfig, figcol= figcol, xshrink = xshrink, yshrink = yshrink))
}

### .onUnload
# Clean up
.onUnload <- function (libpath) {
  library.dynam.unload("DEEE", libpath)
}
