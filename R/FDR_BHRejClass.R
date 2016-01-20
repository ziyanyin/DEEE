#####################################################
### This file defines BHRej class and its methods ###
#####################################################
### functions:
##      BHRej, defines the BHRej class.
##      print.BHRej, print method for BHRej class.
##      PIE, A generic function producing pie plots.
##      PIE.BHRej, pie plot method for BHRej class.

### BHRej
#' Produce a BHRej object
#'
#' \code{BHRej} implements the Bonjamini-Hocheberg procedure given a vector of p-values.
#' Heap sort algorithm is used to improve efficiency. See also \code{\link{BHRejVec}}.
#'
#' @param FDR A numeric value between 0 and 1 for FDR control.
#' @param pvector A vector of numerics containing pvalues.
#' @return A BHRej object is return, which is list containing:
#' \describe{
#'  \item{rejList}{A vector of integers containing the relative indexes of rejected genes.}
#'  \item{pVal}{A vector of numerics containing the tested p-values.}
#'  \item{FDR}{A numeric value between 0 and 1.}
#'  \item{nTest}{An integer indicating the number of p-values.}
#' }
#' @examples
#' BHRej(FDR = 0.1, runif(10000, 0, 0.1))
#' @export
BHRej = function(FDR = 0.1, pvector)
{
    rejList = Cpp_FDR(pvector, FDR)

    myRes = list(rejList = rejList, pVal = pvector, FDR = FDR, nTest = length(pvector))
    class(myRes) = "BHRej"
    return(myRes)
}

### print.BHRej
#' @export
print.BHRej = function(obj)
{
	rejl = length(obj$rejList)
    cat("This is a BHRej object.\n")
	tmp = data.frame(FDR = format(obj$FDR), Rej = rejl, Total = obj$nTest, RejRate = paste0(format(rejl / obj$nTest * 100, digits = 4), "%"))
    print(tmp, row.names = FALSE, prefix = "\t\t\t")
}

### PIE
#' Generic function of PIE plot
#'
#' \code{PIE} produces a pie plot. More details could be found in \code{\link{PIE.BHRej}}
#' and \code{\link{PIE.BHRejVec}}.
#'
#' @param x An object.
#' @param ... Others parameters.
#' @export
PIE = function(x, ...) UseMethod("PIE")

### PIE.BHRej
#' PIE plot
#'
#' \code{PIE.BHRej} produces a pie plot for objects of BHRej class. The first argument
#' contains imformation of DE and the second EE. They should be produced with the same data
#'  set, but different test types.
#'
#' @param DE,EE Objects of BHRej class.
#' @param graPar A list of parameters adjusting global graphic.
#' @param main A string specifying title.
#' @param piePar A list of parameters adjusting the pie plot.
#' @examples
#' DE = BHRej(FDR = 0.1, runif(1000, 0, 0.1))
#' EE = BHRej(FDR = 0.1, runif(1000, 0, 0.1))
#' PIE(DE = DE, EE = EE, piePar = list(radius = 0.9, main = "What"))
#' @export
PIE.BHRej = function(DE, EE, graPar = NULL, main =  NULL, piePar = list())
{
    myDE = DE
    myEE = EE

    if(!is(myDE, "BHRej") || !is(myEE, "BHRej")) stop("Invalid input object")

    parDefault = par(graPar)
    on.exit(par(parDefault))

    FDRde = myDE$FDR
    FDRee = myEE$FDR
    FDR = FDRee
    if(FDRde != FDRee) stop("FDR values of DE and EE are not the same.")

    nDE = myDE$nTest
    nEE = myEE$nTest
    if(nDE != nEE) stop("Different length of genes.")
    n = nDE
    n2 = length(myEE$rejList)
    n3 = length(myDE$rejList)
    n1 = n - n2 - n3
	if(n1 < 0) stop("Numbers of DE and EE genes is greater than the total number.")

    do.call("pie", paraMerge(list(x = c(n1, n2, n3), labels = paste0(format(round(c(n1, n2, n3) / n * 100, 1)), "%"), main = paste0("FDR: ", FDR)), piePar))
    invisible()
}
