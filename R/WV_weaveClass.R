#####################################################
### This file defines weave class and its methods ###
#####################################################
#' @include AA_samVecClass.R
NULL

### functions:
##      weave, defines the weave class.
##      medWV, produces a weave object with certain restrictions.
##      print.weave, print method for weave class.
##      plot.weave, the plot method for weave class.

### weave
#' Produce a weave object.
#'
#' \code{weave} defines and produces a weave object. The difference between weave and
#' \code{\link{grpWV}} is their way to apply FUN. In weave, FUN is applied to every
#' single column; while in grpWV, FUN takes the whole columns in a group as argument.
#' See also \code{\link{medWV}}, in which a weave object also is produced, but with
#' FUN = ``median'', interval = c(0, 1) and ifeq = ``>''. See \code{\link{plot.weave}}
#' for plotting instructions.
#'
#' @param dataset A samVec object containing data and other information.
#' @param degree A integer specifying the degree of the polynomail logistic regression that
#' applied to fit the data.
#' @param FUN A string specifying the function used in weave.
#' @param interval A value between 0 and 1 to restrict the rankit used.
#' @param ifeq An bool value indicating whether >= or > is used in comparison. This argument
#' is only valid when it comes tied values.
#' @return A weave object is returned, which is a list containing
#' \describe{
#'  \item{fitY}{A list of lists of numerics containing the fitted response values.}
#'  \item{rankit}{A list of numerics containing the rankit of the means of the data matrix.}
#'  \item{order}{A list of integers containing the order the rankit.}
#'  \item{nWV}{An integers indicating the numbers of weave.}
#'  \item{interval}{An interval indicating the range of rankit that is used.}
#'  \item{Function}{A string indicating the function used in weave.}
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
#' testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"),
#'  dataType = "Example")
#' medw1 = weave(dataset = testset, degree = 6, FUN = "median", interval = c(0, 1))
#' plot(medw1, legPar = list(cex = 1, ncol = 1), plotPar = list(main = "Example"))
#' @export
weave = function(dataset, degree = 4, FUN = "median", interval = c(0, 1), ifeq = FALSE)
{
    if(!is(dataset, "samVec")) stop("Invalid dataset input")
    mydata = dataset$data

    nWV = dataset$nCol
    nrow = dataset$nRow

    newX = MatRankit(mydata)
    myOd = Cpp_order(newX) + 1

    selInd = which((newX <= max(interval)) & (newX >= min(interval)))
    newX = newX[selInd]
    polyx = poly(newX, degree = degree, raw = TRUE)

    fitY = lapply(1:nWV, function(idcol) {
        cat(" ", idcol, "\t out of ", nWV, "\n")
        y = sapply(1:nrow, function(i) as.numeric(get(ifelse(ifeq, ">=", ">"))(mydata[i, idcol], match.fun(FUN)(mydata[i, -idcol]))))
        newY = y[selInd]
        predY = predict(glm(newY ~ polyx, family = "binomial"), type = "response")
        predY
        })
    names(fitY) = 1:nWV

    res = list(fitY = fitY, rankit = newX, order = myOd, nWV = nWV, interval = interval, Function = FUN, ifeq = ifelse(ifeq, ">=", ">"), selCol = dataset$selCol, colInd = dataset$colInd, labels = dataset$labels, dataType = dataset$dataType)
    class(res) = "weave"
    return(res)
}

### medWV
#' Produce a weave object with restrictions.
#'
#' \code{medWV} defines and produces a weave object with FUN = ``median'',
#' interval = c(0, 1) and ifeq = ``>''. See also \code{\link{weave}}, which produces a weave
#' object with no restriction but is much slower. See \code{\link{plot.weave}}
#' for plotting instructions.
#'
#' @param dataset A samVec object containing data and other information.
#' @param degree A integer specifying the degree of the polynomail logistic regression that
#' applied to fit the data.
#' @return A weave object is returned, which is a list containing
#' \describe{
#'  \item{fitY}{A list of lists of numerics containing the fitted response values.}
#'  \item{rankit}{A list of numerics containing the rankit of the means of the data matrix.}
#'  \item{order}{A list of integers containing the order the rankit.}
#'  \item{nWV}{An integers indicating the numbers of weave.}
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
#' m1 = medWV(dataset = testset, degree = 6)
#' plot(m1, legPar = list(cex = 1, ncol = 1), plotPar = list(main = "Example"))
#' @export
medWV = function(dataset, degree = 4)
{
    if(!is(dataset, "samVec")) stop("Invalid dataset input")
    mydata = dataset$data

    nWV = dataset$nCol
    nrow = dataset$nRow
    newX = MatRankit(mydata)
    myOd = Cpp_order(newX) + 1

	res = list()
	polyx = poly(newX, degree = degree, raw = TRUE)
	fitY = lapply(1:nWV, function(ref) {
		cat(" ", ref, "\t out of ", nWV, "\n")
		y = MedWVY(mydata, ref - 1)
		predict(glm(y ~ polyx, family = "binomial"), type = "response")
	})

	names(fitY) = 1:nWV

	res = list(fitY = fitY, rankit = newX, order = myOd, nWV = nWV, interval = c(0, 1), Function = "median", ifeq = ">", selCol = dataset$selCol, colInd = dataset$colInd, labels = dataset$labels, dataType = dataset$dataType)
    class(res) = "weave"
    return(res)
}

### print.weave
#' @export
print.weave = function(x, ...)
{
    obj = x
    cat("This is a weave object.\n")
    cat("Data type: ", obj$dataType, ". Function:", obj$Function, ",", obj$ifeq, ". ")
	cat("Rankit interval: (", paste0(obj$interval, collapse = ", "), "). \n")
    cat(paste0("Total Weaves: ", obj$nWV, " in ", length(obj$labels), " Groups."), "\n")
    cat(paste0("Contains: ", paste0(obj$labels, collapse = ", "), "."))
}

### plot.weave
#' Plot method for weave class
#'
#' \code{plot.weave} is the plot method for weave class.
#'
#' @param x A weave object.
#' @param indSet A vector of integers indicating which weaves are plotted.
#' @param legPar A list of legend parameters.
#' @param graPar A list of global graphics parameters.
#' @param plotPar A list of background (including main title) parameters.
#' @param ... ignored
#' @examples
#' data(GCwPADataA)
#' testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"))
#' m1 = medWV(dataset = testset, degree = 6)
#' plot(m1, legPar = list(cex = 1, ncol = 1), plotPar = list(main = "Example"))
#' @export
plot.weave = function(x, indSet = 1:obj$nWV, legPar = list(cex = 0.6, ncol = 3), graPar = list(), plotPar = list(), ...)
{
    obj = x
    parDe = par(graPar)
    on.exit(par(parDe))

    do.call("plot", paraMerge(list(x = obj$interval, y = c(0.5, 0.5), xlab = "Rankit", ylab = paste0("Pr(sample", obj$ifeq, obj$Function, "(other samples))"), xlim = obj$interval, ylim = c(0, 1), type = "l", col = "gray", lty = 10), plotPar))

	fitY = obj$fitY
    x = obj$rankit
    int = obj$interval
	myOd = obj$order

	newX = x[myOd]
	refset = rep(1:length(obj$colInd), sapply(obj$colInd, length))
	sapply(indSet, function(i) {
		lines(x = newX, y = fitY[[i]][myOd], col = refset[i])
		})
    grid()

    if(is.null(legPar$ifleg) || legPar$ifleg)
    {
        legLoc = NULL
        if(is.null(legPar$x)) legLoc = c(mean(int), 1.04)
        legPar$ifleg = NULL
        do.call("legend", paraMerge(list(x = legLoc[1], y = legLoc[2], legend = obj$labels, lty = 1, col = 1:length(obj$labels)), legPar))
    }

    invisible()
}
