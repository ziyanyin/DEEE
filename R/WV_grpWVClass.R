#####################################################
### This file defines grpWV class and its methods ###
#####################################################
#' @include AA_samVecClass.R
NULL

### functions:
##      grpWV, defines the grpWV class.
##      madWV, produces a grpWV object with certain restrictions.
##      print.grpWV, print method for grpWV class.
##      plot.grpWV, the plot method for grpWV class.

### grpWV
#' Produce a grpWV object
#'
#' \code{grpWV} defines and produces a grpWV object. The difference between grpWV and
#' \code{\link{weave}} is their way to apply FUN. In weave, FUN is applied to every
#' single column; while in grpWV, FUN takes the whole columns in a group as argument.
#' See also \code{\link{madWV}}, in which a grpWV object also is produced, but with
#' FUN = ``mad'', interval = c(0, 1) and ifeq = ``>''. See \code{\link{plot.grpWV}}
#' for plotting instructions.
#'
#' @param SVobj A samVec object containing data and other information.
#' @param degree A integer specifying the degree of the polynomail logistic regression that
#' applied to fit the data.
#' @param interval A value between 0 and 1 to restrict the rankit used.
#' @param FUN A string specifying the function used in grpWV.
#' @param ifeq An bool value indicating whether >= or > is used in comparison. This argument
#' is only valid when it comes tied values.
#' @return A grpWV object is returned, which is a list containing
#' \describe{
#'  \item{fitY}{A list of lists of numerics containing the fitted response values.}
#'  \item{rankit}{A list of numerics containing the rankit of the means of the data matrix.}
#'  \item{order}{A list of integers containing the order the rankit.}
#'  \item{nWV}{An integers indicating the numbers of grpWV.}
#'  \item{interval}{An interval indicating the range of rankit that is used.}
#'  \item{Function}{A string indicating the function used in grpWV.}
#'  \item{ifeq}{An bool value indicating whether >= or > is used in comparison. This
#'   argument is only valid when it comes tied values.}
#'  \item{selCol}{A list of vector of integers containing absolute indexes of choosen
#'  columns.}
#'  \item{colInd}{A list of vector of integers containing relative indexes of columns.}
#'  \item{labels}{A vector of characters, integers or short string marking the groups.}
#'  \item{dataType}{A long string containing useful information of the data.}
#' }
#' @examples
#' data(GCwPADataA)
#' t1 = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"),
#'  dataType = "Example")
#' madRes = grpWV(SVobj = t1, degree = 6, interval = c(0, 1), FUN = "mad", ifeq = TRUE)
#' plot(madRes, legPar = list(cex = 1, ncol = 1), mainPar = list(main = "Example"))
#' @export
grpWV = function(SVobj, degree = 6, interval = c(0, 1), FUN = "mad", ifeq = FALSE)
{
    if(!is(SVobj, "samVec")) stop("Invalid SVobj input")
    mydata = SVobj$data
	selCol = SVobj$selCol
	colInd = SVobj$colInd
    nWV = SVobj$nGroup
    nrow = SVobj$nRow

    newX = MatRankit(mydata)
    myOd = Cpp_order(newX) + 1

    selInd = which((newX <= max(interval)) & (newX >= min(interval)))
    newX = newX[selInd]
    polyx = poly(newX, degree = degree, raw = TRUE)

    fitY = lapply(1:nWV, function(idcol) {
        y = sapply(1:nrow, function(i) as.numeric(get(ifelse(ifeq, ">=", ">")) (match.fun(FUN)(mydata[i, colInd[[idcol]]]), match.fun(FUN)(mydata[i, -colInd[[idcol]]]))))
        newY = y[selInd]
        predict(glm(newY ~ polyx, family = "binomial"), type = "response")
        })
    names(fitY) = 1:nWV
    res = list(fitY = fitY, rankit = newX, order = myOd, nWV = nWV, interval = interval, Function = FUN, ifeq = ifelse(ifeq, ">=", ">"), selCol = SVobj$selCol, colInd = SVobj$colInd, labels = SVobj$labels, dataType = SVobj$dataType)
    class(res) = "grpWV"
    return(res)
}

### madWV
#' Produce a grpWV object
#'
#' \code{madWV} defines and produces a grpWV object with FUN = ``mad'', interval = c(0, 1)
#' and ifeq = ``>''. See also \code{\link{grpWV}}, which produces a grpWV object with no
#' restriction but is much slower. See \code{\link{plot.grpWV}}
#' for plotting instructions.
#'
#' @param SVobj A samVec object containing data and other information.
#' @param degree A integer specifying the degree of the polynomail logistic regression that
#' applied to fit the data.
#' @return A grpWV object is returned, which is a list containing
#' \describe{
#'  \item{fitY}{A list of lists of numerics containing the fitted response values.}
#'  \item{rankit}{A list of numerics containing the rankit of the means of the data matrix.}
#'  \item{order}{A list of integers containing the order the rankit.}
#'  \item{nWV}{An integers indicating the numbers of grpWV.}
#'  \item{interval}{c(0, 1).}
#'  \item{Function}{"mad".}
#'  \item{ifeq}{">".}
#'  \item{selCol}{A list of vector of integers containing absolute indexes of choosen
#'  columns.}
#'  \item{colInd}{A list of vector of integers containing relative indexes of columns.}
#'  \item{labels}{A vector of characters, integers or short string marking the groups.}
#'  \item{dataType}{A long string containing useful information of the data.}
#' }
#' @examples
#' data(GCwPADataA)
#' testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"),
#'  dataType = "Example")
#' gw1 = madWV(SVobj = testset, degree = 6)
#' plot(gw1, legPar = list(cex = 1, ncol = 1), mainPar = list(main = "Example"))
#' @export
madWV = function(SVobj, degree = 4)
{
    if(!is(SVobj, "samVec")) stop("Invalid SVobj input")
    mydata = SVobj$data
	selCol = SVobj$selCol
	colInd = SVobj$colInd
    nWV = SVobj$nGroup
    nrow = SVobj$nRow

    newX = MatRankit(mydata)
    myOd = Cpp_order(newX) + 1

	res = list()
	polyx = poly(newX, degree = degree, raw = TRUE)
	fitY = lapply(1:nWV, function(ref) {
		cat(" ", ref, "\t out of ", nWV, "\n")
		y = MADWVY(mydata, colInd[[ref]] - 1)
		predict(glm(y ~ polyx, family = "binomial"), type = "response")
	})
	names(fitY) = 1:nWV

	res = list(fitY = fitY, rankit = newX, order = myOd, nWV = nWV, interval = c(0, 1), Function = "mad", ifeq = ">", selCol = SVobj$selCol, colInd = SVobj$colInd, labels = SVobj$labels, dataType = SVobj$dataType)
    class(res) = "grpWV"
    return(res)
}

### print.grpWV
#' @export
print.grpWV = function(x, ...)
{
    obj = x
    cat("This is a grpWV object.\n")
    cat("Data type: ", obj$dataType, ". Function:", obj$Function, ",", obj$ifeq, ".")
	cat("Rankit interval: (", paste0(obj$interval, collapse = ", "), "). \n")
    cat(paste0("Total Weaves: ", obj$nWV, ".\n"))
    cat(paste0("Contains: ", paste0(obj$labels, collapse = ", "), ".\n"))
}

### plot.grpWV
#' Plot method for grpWV class
#'
#' \code{grpWV} is the plot method for grpWV class.
#'
#' @param x A grpWV object.
#' @param indSet A vector of integers indicating which weaves are plotted.
#' @param legPar A list of legend parameters.
#' @param graPar A list of global graphics parameters.
#' @param mainPar A list of background (including main title) parameters.
#' @param colSet,ltySet Vectors of numerics or characters specifying col types or line
#' @param ... ignored
#' @examples
#' data(GCwPADataA)
#' testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"))
#' madw1 = madWV(SVobj = testset, degree = 6)
#' plot(madw1, legPar = list(cex = 1, ncol = 1), mainPar = list(main = "Example"))
#' @export
plot.grpWV = function(x, indSet = NULL, legPar = list(cex = 0.6, ncol = 3), graPar = list(), mainPar = list(), colSet = NULL, ltySet = NULL, ...)
{
    obj = x
	if(is.null(indSet)) indSet = 1:obj$nWV
    refSet = 1:length(obj$colInd)
	if(is.null(colSet)) colSet = refSet
	if(is.null(ltySet)) ltySet = refSet
    parDe = par(graPar)
    on.exit(par(parDe))

    do.call("plot", paraMerge(list(x = obj$interval, y = c(0.5, 0.5), xlab = "Rankit", ylab = paste0("Pr(", obj$Function, "(samples)", obj$ifeq, obj$Function, "(other samples))"), xlim = obj$interval, ylim = c(0, 1), type = "l", col = "gray", lty = 10), mainPar))
    if(is.null(legPar$ifleg) || legPar$ifleg)
    {
        legLoc = NULL
        if(is.null(legPar$x)) legLoc = c(mean(obj$interval), 1.04)
        legPar$ifleg = NULL
        do.call("legend", paraMerge(list(x = legLoc[1], y = legLoc[2], legend = obj$labels, lty = ltySet, col = colSet, xpd = NA), legPar))
    }

	fitY = obj$fitY
    x = obj$rankit
	myOd = obj$order

	newX = x[myOd]
	sapply(indSet, function(i) {
		lines(x = newX, y = fitY[[i]][myOd], col = colSet[i], lty = ltySet[i])
		})
    grid()
    invisible()
}

