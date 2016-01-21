####################################################
### This file describe PCA class and its methods ###
####################################################
#' @include AA_samVecClass.R
NULL

### functions:
##      PCA, defines the PCA class.
##      print.PCA, print method for PCA class.
##      plot.PCA, the plot method for PCA class.

### PCA
#' Produce a PCA object.
#'
#' \code{PCA} defines and produces a PCA object. See \code{\link{plot.PCA}}
#' for plotting instructions.
#'
#' @param SVobj A samVec object containing data and other information.
#' @param scale A bool value indicating whether scaling is used or not?
#' @return A PCA object is returned, which is list containing:
#' \describe{
#'  \item{pcx1}{A matrix containing the first principle conponent.}
#'  \item{pcx2}{A matrix containing the second principle conponent.}
#'  \item{colInd}{A list of vector of integers containing relative indexes of columns.}
#'  \item{labels}{A vector of characters, integers or short string marking the groups.}
#'  \item{dataType}{A long string containing useful information of the data.}
#'  \item{nGroup}{An interger indicating the number of groups specified by selCol.}
#'  \item{nCol}{An integer indicating the number of columns of data.}
#'  \item{nRow}{A integers indicating the number of rows of data.}
#' }
#' @examples
#' data(GCwPADataA)
#' testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"))
#' pca1 = PCA(testset)
#' plot(pca1, mainPar = list(labels = "Example", y = 1.55), pcaPar = list(cex = 1, font = 2))
#' @export
PCA = function(SVobj, scale = F) {
    if(!is(SVobj, "samVec")) stop("Invalid SVobj input")
	X = SVobj$data
	fit = prcomp(t(X), scale = scale)
	pcx1 = t(X)%*%fit$rotation[,1]
	pcx2 = t(X)%*%fit$rotation[,2]
	res = list(pcx1 = pcx1, pcx2 = pcx2, colInd = SVobj$colInd, labels = SVobj$labels, dataType = SVobj$dataType, nGroup = SVobj$nGroup, nCol = SVobj$nCol, nRow = SVobj$nRow)
	class(res) = "PCA"
	return(res)
}

### Print methods of class PCA
#' @export
print.PCA = function(x, ...)
{
    obj = x
    cat("This is a PCA object\n")
    cat("Sample:", obj$dataType, ". Group Labels:", paste0(obj$labels, collapse = ", "), ".\n")
    cat(paste0("Group size: ", obj$nGroup, ". Total Columns: ", obj$nCol, ". Number of Rows: ", obj$nRow, ".\n"))
}

### plot.PCA
#' The plot method for a PCA object.
#'
#' \code{plot.PCA} plots all the points given their first and second components. Different
#' points in different groups are given different colors and labels.
#'
#' @param x A PCA object.
#' @param mainPar A list of parameters specifying the main title.
#' @param graPar A list of global graphics parameters.
#' @param pcaPar A list of parameters adjusting the PCA points.
#' @param ... ignored
#' @examples
#' data(GCwPADataA)
#' testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"))
#' pca1 = PCA(testset)
#' plot(pca1, mainPar = list(labels = "Example", y = 1.55), pcaPar = list(cex = 1, font = 2))
#' @export
plot.PCA = function(x, mainPar = list(), graPar = list(), pcaPar = list(), ...) {
    obj = x
    pcx1 = obj$pcx1
    pcx2 = obj$pcx2
	plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
	do.call("text", paraMerge(list(x = 1, y = 1.5, labels = paste0("PCA ", obj$dataType), font = 2, cex = 2.5, xpd = NA), mainPar))

	parDe = paraMerge(list(new = par("new")), par(graPar))
	on.exit(par(parDe))

    par(new = TRUE)
	plot(pcx1, pcx2, pch = " ", main = "", xlab = "", ylab = "")
	do.call("text", paraMerge(list(x = pcx1, y = pcx2, labels = rep(obj$labels, sapply(obj$colInd, length)), col = rep(1:obj$nGroup, sapply(obj$colInd, length)), cex = 2), pcaPar))
	invisible()
}
