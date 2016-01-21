###############################################################
### This file defines Class of samVec and its methods ###
###############################################################
### functions:
##      samVec, define the samVec class.
##      print.samVec, print method for samVec class.
##      plot.samVec, plot method for samVec class.

### samVec
#' Produce a samVec object
#'
#' \code{samVec} defines the samVec object given necessary data and information.
#' The objects of samVec class are the base for further analysis. 
#' See \code{\link{plot.samVec}} for plotting instructions.
#'
#' @param datMat A matrix.
#' @param selCol A list of vector of integers containing absolute indexes of choosen columns.
#' @param labels A vector of characters or short strings marking the groups. The length of
#' labels should be the same with the length of selCol. If not specifying, a vector of
#' integer is given.
#' @param dataType A string contains data type or other useful imformation.
#' @return \code{samVec} returns a samVec object. A samVec is a list contains:
#' \describe{
#'  \item{data}{A matrix containing the data imported by datMat.}
#'  \item{nCol}{An integer indicating the number of columns of data.}
#'  \item{nGroup}{An interger indicating the number of groups specified by selCol.}
#'  \item{colInd}{A list of vector of integers containing relative indexes of columns.}
#'  \item{selCol}{A list of vector of integers containing absolute indexes of choosen columns.}
#'  \item{nRow}{A integers indicating the number of rows of data.}
#'  \item{labels}{A vector of characters, integers or short string marking the groups.}
#'  \item{dataType}{A long string containing useful information of the data.}
#' }
#' @examples
#' datamatrix = matrix(rnorm(100000), 1000, 100)
#' testset = samVec(datamatrix, selCol = list(1:5, 11:15, 21:25), dataType = "Example")
#' plot(testset, main = testset$dataType)
#' @export
samVec = function(datMat, selCol = NULL, labels = NULL, dataType = "MAQC")
{
    if(!is.matrix(datMat)) stop("datMat has to be a matrix")
    datMat = datMat[, unlist(selCol)]
    nCol = ncol(datMat)
    nGroup = length(selCol)
    if(!is.list(selCol)) stop("selCol must be a list.")
    if(is.null(labels))
    {
        labels = 1:nGroup
    }

	colInd = split(1:nCol, factor(rep(1:length(selCol), sapply(selCol, length))))
    if(length(labels) != nGroup)
    {
        warning("Unequal numbers of labels and selCol.")
        labels = rep_len(labels, nGroup)
    }
    testset = list(data = datMat, nCol = nCol, nGroup = nGroup)
  	testset$colInd = colInd
    testset$nRow = nrow(datMat)
    testset$labels = labels
    testset$dataType = dataType
    testset$selCol = selCol
    class(testset) = "samVec"

    return(testset)
}

#' @export
print.samVec = function(x, ...)
{
    obj = x
    cat("This is a samVec object\n")
    cat(paste0("Data Type: ", obj$dataType, ". Selected Columns: ", paste0(obj$selCol, collapse = ", ")), ".\n")
    tmpstr = data.frame(paste(obj$labels, "   "), I(paste0(obj$nRow, " * ", sapply(obj$selCol, length))))
    names(tmpstr) = c("Group Labels", "Group Dim")
    cat(paste0("Group size: ", length(obj$selCol), ". Total Columns: ", obj$nCol, ". Number of Rows: ", obj$nRow, ".\n"))
    print(tmpstr)
}

#' Plot methods of class samVec
#'
#' \code{plot.samVec} divides the data by their groups and calcutes their densities
#' separately. Their densities are plotted together in the same figure.
#'
#' @param x A samVec object.
#' @param graPar A list of parameters adjusting the global graphics.
#' @param colSet A vector of strings or integers specifying colors.
#' @param ltySet A vector of integers specifying line types.
#' @param main A string specifying main title.
#' @param legPar A list of parameters specifying the legend.
#' @param ... ignored
#' @examples
#' datamatrix = matrix(rnorm(100000), 1000, 100)
#' testset = samVec(datamatrix, selCol = list(1:5, 11:15, 21:25), dataType = "Example")
#' plot(testset, main = testset$dataType)
#' @export
plot.samVec = function(x, graPar = NULL, colSet = 1:x$nGroup, ltySet = rep(1, x$nGroup), main = "KDE", legPar = list(), ...)
{
    obj = x
    parDe = paraMerge(list(xpd = par("xpd")), par(graPar))
    on.exit(par(parDe))
    myData = obj$data
    tmpRes = lapply(obj$colInd, function(x) density(as.vector(myData[, x])))

    plot(0, 0, xlim = range(sapply(tmpRes, function(tmpx) tmpx$x)), ylim = range(sapply(tmpRes, function(tmpx) tmpx$y)), type = "n", xlab = "", ylab = "Densities", main = main)
    if(is.null(legPar$ifleg) || legPar$ifleg)
    {
        legPar$ifleg = NULL
        do.call("legend", c(list("topright", legend = paste0("site ", obj$labels), col = colSet, lty = ltySet), legPar))
    }

    sapply(seq_along(tmpRes), function(i) {
        lines(tmpRes[[i]], col = colSet[i], lty = ltySet[i])
    })

    par(xpd = FALSE)
    grid()

    invisible()
}
