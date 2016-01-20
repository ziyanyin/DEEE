#############################################################
### This file defines the necessary function for heat map ###
#############################################################
### function:
##      DEEEHeat, the basic heat map function.
##      DEEEAxis, the function plotting axes.
##      DEEELegend, the function plotting axes.
##      HeatMap, the function combining DEEEHeat, DEEEAxis and DEEELegend together.

### DEEEHeat
#' Core function for heatmap.
#'
#' \code{DEEEHeat} is the basic function to construct a heatmap. This function plot the
#' core heat map and
#' return the coordinates containing the border of the plot.
#'
#' @param mydata A matrix containing data.
#' @param mycol A vector of strings indicating color range.
#' @param graPar A list of global graphics parameters.
#' @param textNum A list specifying texts plotted in each cell of the heatmap.
#' NA for no number.
#' @param textX A list containing parameters specifing X-axis. NA for no plotting.
#' @param textY A list containing parameters specifing Y-axis. NA for no plotting.
#' @param Seg A bool value indicating whether plot separating lines or not.
#' @return A vector of numerics indicating the border of the plot, (xleft, xright, ybottom,
#' ytop) is returned.
#' @examples
#' testMatrix = matrix(sample(1:100), 10, 10)
#' colnames(testMatrix) = 1:10
#' rownames(testMatrix) = 1:10
#' DEEEHeat(testMatrix)
#' @export
DEEEHeat = function(mydata, mycol = NULL, graPar = list(), textNum = NULL, textX = NULL, textY = NULL, Seg = TRUE)
{
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

    nx = ncol(mydata)
    ny = nrow(mydata)

    vdata = as.vector(mydata)

    if(is.null(mycol)) mycol = c("red", "yellow", "green")

    colFun = colorRampPalette(mycol)
    totalCol = colFun(round(max(vdata)))

    total = nx * ny
    nxseq = rep(seq(1, nx), each = ny)
    nyseq = rep(seq(ny, 1), nx)

    parDe = par(c(graPar, list(new = TRUE)))
    on.exit(par(parDe))

    idata = t(apply(t(mydata), 1, rev))
    image(x = 0:nx, y = 0:ny, z = idata, col = totalCol, xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(-2, nx + 2), ylim = c(-2, ny + 2), bty = "n")

    if(Seg) segments(c(0:nx, rep(0, ny + 1)), c(rep(0, nx + 1), 0:ny), c(0:nx, rep(nx, ny + 1)), c(rep(ny, nx + 1), 0:ny), col = "white")
    if(!(is.null(textNum) || is.na(textNum))) do.call("text", paraMerge(list(x = nxseq - 0.5, y = nyseq - 0.5, labels = round(vdata, 1)), textNum))
    if(!(is.null(textX) || is.na(textX))) do.call("text", paraMerge(list(x = 1:nx - 0.5 , y = -1, labels = rownames(mydata), xpd = NA), textX))
    if(!(is.null(textY) || is.na(textY))) do.call("text", paraMerge(list(x = -2, y = seq(nrow(mydata), 1) - 0.5, labels = rownames(mydata), xpd = NA), textY))

    return(c(0, ny, 0, nx))
}

### DEEEAxis
#' Add axis
#'
#' \code{DEEEAxis} adds axis to existing heapmaps.
#'
#' @param lm A numeric indicating left margin.
#' @param rm A numeric indicating right margin.
#' @param tm A numeric indicating top margin.
#' @param bm A numeric indicating bottom margin.
#' @param sep A vector of numeric used as delimiter of axis columns.
#' @param labels A vector of strings added to each separate part.
#' @param col A vector of strings specifying the colors added to each separate part.
#' @param bdPar A list of parameters specifying the borders of axis.
#' @param txtPar A list of parameters specifying the labels in the axis.
#' @examples
#' testMatrix = matrix(sample(1:100), 10, 10)
#' colnames(testMatrix) = 1:10
#' rownames(testMatrix) = 1:10
#' DEEEHeat(testMatrix)
#' DEEEAxis(lm = 0, rm = 10, bm = 10, tm = 11, sep = c(1, 3, 5), col = 2:5)
#' @export
DEEEAxis = function(lm = 0, rm, tm, bm, sep = numeric(), labels = letters[1:(length(sep) + 1)], col = NULL, bdPar = list(ifbd = TRUE, col = "white"), txtPar = NULL)
{
    nsep = length(sep)
    if((nsep > 0) && (max(sep) >= rm)) warning("Sep right margin outside!")
    rect(c(lm, sep), rep(bm, nsep + 1), c(sep, rm), rep(tm, nsep + 1), col = col, border = NA, xpd = NA)
    if(is.null(bdPar$ifbd) || bdPar$ifbd) {
        bdPar$ifbd = NULL
        do.call("segments", paraMerge(list(x0 = lm, y0 = bm, x1 = rm, y1 = bm, xpd = NA, col = "white"), bdPar))
        do.call("segments", paraMerge(list(x0 = lm, y0 = bm, x1 = lm, y1 = tm, xpd = NA, col = "white"), bdPar))
        do.call("segments", paraMerge(list(x0 = lm, y0 = tm, x1 = rm, y1 = tm, xpd = NA, col = "white"), bdPar))
        do.call("segments", paraMerge(list(x0 = rm, y0 = bm, x1 = rm, y1 = tm, xpd = NA, col = "white"), bdPar))
        }

    if(nsep > 0)
    {
        if(is.null(bdPar$ifbd) || bdPar$ifbd) {
            bdPar$ifbd = NULL
            do.call("segments", paraMerge(list(x0 = sep, y0 = rep(bm, nsep), x1 = sep, y1 = rep(tm, nsep), xpd = NA, col = "white"), bdPar))
        }

        tmp = c(lm, rep(sep, each = 2), rm)
        mPos = vector("numeric", nsep + 1)
        for(i in 1:(nsep + 1)) mPos[i] = (tmp[i * 2 - 1] + tmp[i * 2]) / 2

        do.call("text", paraMerge(list(x = mPos, y = rep((bm + tm) / 2, nsep), labels = labels, xpd = NA), txtPar))
    }
    else
    {
        do.call("text", paraMerge(list(x = (lm + rm) / 2, y = (bm + tm) / 2, labels = labels, xpd = NA), txtPar))
    }

    invisible()
}

### DEEELegend
#' Add legends
#'
#' \code{DEEELegend} adds legends to existing heatmap.
#'
#' @param pos A vector of length four specifying (xleft, ybottom, xright, ytop),
#' indicating the position of legend.
#' @param mycol A vector of strings giving the color range of legend. Usually not need to
#' be specified.
#' @param cex A numeric indicating the cex of text in legend.
#' @param labels A vector of strings added to the legend.
#' @param textPos A list of two vectors of numerics indicating the position of the label
#' positions.
#' @examples
#' testMatrix = matrix(sample(1:100), 10, 10)
#' colnames(testMatrix) = 1:10
#' rownames(testMatrix) = 1:10
#' DEEEHeat(testMatrix)
#' DEEEAxis(lm = 0, rm = 10, bm = 10, tm = 11, sep = c(1, 3, 5), col = 2:5)
#' DEEELegend(pos = c(10.5, 5, 11, 10), labels = c(0, 50, 100), textPos = list(x = c(11.5,
#'  11.5, 11.5), y = c(5, 7.5, 10)))
#' @export
DEEELegend = function(pos, textPos = NULL, cex = 1, labels, mycol = NULL)
{
    if(is.null(mycol)) mycol = c("red", "yellow", "green")
    legCol = colorRampPalette(mycol)(100)
    yp = seq(pos[2], pos[4], length.out = 101)

    image(x = pos[c(1, 3)], y = yp, z = matrix(1:100, 1, 100), col = legCol, xpd = NA, bty = "n", xaxt = "n", yaxt = "n", add = T, xlab  ="", ylab = "")

    #rect(pos[1], pos[2], pos[3], pos[4], xpd = NA)
    if(is.null(textPos)) textPos = list(x = c(pos[3] / 2 + pos[1] / 2, pos[3] / 2 + pos[1] / 2), y = c(pos[4] + 0.5, pos[2] - 0.5))
    text(textPos, labels = labels, xpd = NA, cex = cex)

    invisible()
}

### HeatMap
#' Plot heat map
#'
#' \code{HeatMap} constructs a heapmat by combining DEEEHeat, DEEEAxis and DEEELegend
#' together based on specified arguments. It is needed to pay attention that all paramters
#'  should be lists.
#'
#' @param mydata A matrix containing data.
#' @param mainPar A list of parameters of the title.
#' @param heatPar A list of parameters specifying the main heatmap.
#' @param legPar A list of parameters of specifying DEEELegend.
#' @param axisPar A list of lists of parameters of axes. Every single list in it specifying
#' one axis.
#' @param globalPar A list sepecifying for the global graphical parameters.
#' @return A vector of numerics indicating the border of the core plot, (xleft, xright,
#' ybottom, ytop) is returned.
#' @examples
#' HeatMap(mydata = matrix(sample(1:100), 10, 10),
#'                      mainPar = list(y = 1.65),
#'                      heatPar = list(mycol = c("blue", "yellow", "green")),
#'                      legPar = list(pos = c(10.5, 0, 11, 10), labels = c(0, 25, 50, 100),
#'                                      textPos = list(x = c(rep(11.5, 4)), y = c(0, 2.5, 5.0, 10.0))),
#'                      axisPar = list(a1 = list(sep = c(1, 3, 5), col = 2:4),
#'                                     a2 = list(bm = 10, tm = 11, sep = c(1, 5, 8), col = 4:7)),
#'                      globalPar = list(mar = c(0, 2, 2, 2)))
#' @export
HeatMap = function(mydata, mainPar = list(), heatPar = list(), legPar = list(), axisPar = list(), globalPar = list())
{
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
    do.call("text", paraMerge(list(x = 1, y = 1.6, labels = "MainTitle", cex = 2, xpd = NA), mainPar))
    nx = ncol(mydata)
    ny = nrow(mydata)

    parDe = par(paraMerge(globalPar, list(new = TRUE)))
    on.exit(par(parDe))
    if(is.null(heatPar$mycol)) heatPar$mycol = c("blue", "red", "green")
    mycol = heatPar$mycol
    do.call("DEEEHeat", paraMerge(list(mydata = mydata, mycol = mycol), heatPar))

    legPar = paraMerge(list(pos = c(nx + 0.5, ny / 2, nx + 1, ny), labels = c(round(max(mydata), -1), 0)), legPar)
    legPar$mycol = mycol
    do.call("DEEELegend", legPar)

    if(length(axisPar) != 0)
    {
        for(i in seq_along(axisPar))
        {
            tmpAxisPar = axisPar[[i]]
            tmpAxisPar$lm = ifelse(is.null(tmpAxisPar$lm), 0, tmpAxisPar$lm)
            tmpAxisPar$rm = ifelse(is.null(tmpAxisPar$rm), nx, tmpAxisPar$rm)
            tmpAxisPar$bm = ifelse(is.null(tmpAxisPar$bm), ny + i, tmpAxisPar$bm)
            tmpAxisPar$tm = ifelse(is.null(tmpAxisPar$tm), ny + i + 1, tmpAxisPar$tm)
            do.call("DEEEAxis", tmpAxisPar)
        }
    }
    return(c(0, ny, 0, nx))
}




