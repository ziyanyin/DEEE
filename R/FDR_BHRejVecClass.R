###########################################################
### This file define the BHRejVec class and its methods ###
###########################################################
#' @include FDR_BHRejClass.R
#' @include AA_samVecClass.R
NULL

### functions:
##      BHRejVec, defines the BHRejVec class.
##      print.BHRej, print method for BHRejVec class.
##      PIE.BHRejVec, pie plot method for BHRejVec class.

### BHRejVec
#' Produce a BHRejVec object
#'
#' \code{BHRejVec} defines BHRejVec class and produce a BHRejVec object. A BHRejVec object
#' is a combination of BHRej objects of different FDR values. Also, a samVec object
#' is used to construct a BHRejVec object.
#'
#' @param FDR A vector of numeric values between 0 and 1.
#' @param testset A samVec object containing data and necessary information.
#' @param ifF Test type. When ifF == TRUE, standard F tests are used and; otherwise,
#' equivalence F tests.
#' @param eps A numeric value specifying the equivalence test tolerance. Only valid when
#' ifF == FALSE.
#' @return A BHRejVec object is returned. A BHRejVec is a list contains
#' \describe{
#'  \item{BHProc}{A list of BHRej objects with different FDR values}
#'  \item{labels}{The labels indicating groups specified in the samVec object.}
#'  \item{dataType}{A string the same with the dataType in the samVec object}
#'  \item{testType}{A string of "DE" or "EE" based on ifF or not.}
#'  \item{FDR}{A numeric vector containing the FDR values used.}
#' }
#' @examples
#' data(GCwPADataA)
#' t1 = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"),
#'  dataType = "Example")
#' DE = BHRejVec(FDR = c(0.1, 0.01), testset = t1)
#' EE = BHRejVec(FDR = c(0.1, 0.01), testset = t1, ifF = TRUE)
#' @export
BHRejVec = function(FDR = 0.1, testset, ifF = FALSE, eps = 0.5)
{
    if(!is(testset, "samVec")) stop("Invalid input; must be samVec object")
    datalist = lapply(testset$colInd, function(x) testset$data[, x])
	tmpres = Cpp_fvalues(datalist)
	k = tmpres[1]
	n = tmpres[2]
	nrow = tmpres[3]
	tmpres = tmpres[-(1:3)]
	if(ifF)
	{
		pvalues = sapply(tmpres, function(x) 1 - pf(x, k - 1, n - k))
	}
	else
	{
		pvalues = sapply(tmpres, function(x) pf(x, k - 1, n - k, eps * eps * n / k))
	}
	rejlist = lapply(FDR, function(tmpq) BHRej(FDR = tmpq, pvalues))

    res = list()
    res$BHProc = rejlist

    res$testType = ifelse(ifF, "DE", "EE")
    res$labels = testset$labels
    res$dataType = testset$dataType
	res$FDR = FDR
    class(res) = "BHRejVec"

    return(res)
}

### pring.BHRejVec
#' @export
print.BHRejVec = function(myobj)
{
    DEEE = myobj$testType
    obj = myobj$BHProc
    cat("This is a BHRejVec object for", myobj$dataType, DEEE, "analysis. \n")
    FDR = format(myobj$FDR)
    Rejections = sapply(obj, function(x) length(x$rejList))
    RejRate = paste0(format(Rejections / obj[[1]]$nTest * 100, digits = 4), "%")
    tmp = data.frame(FDR = FDR, Rej = Rejections, Total = sapply(obj, function(sgl) sgl$nTest), RejRate = RejRate)
    print(tmp, row.names = FALSE, prefix = "\t\t\t")
}

#plot.BHRejVec = function(Obj, graPar = NULL, plotPar = list(main = paste0(attr(Obj, "dataType"), " ", attr(Obj, "pool"), ", BH procedure for ", attr(Obj, "equi"))), leg = TRUE, legPar = list(cex = 1))
#{
#    parDefault = par(graPar)
#    on.exit(par(parDefault))

#    myObj = Obj$data
#    obj = myObj[[1]]
#    do.call("plot", c(list(x = 1:length(obj$pvalue[, 2]) * 100 / length(obj$pvalue[, 2]), y = obj$pvalue[, 2], type = "l", xlim = range(1, 100), ylim = c(0, 1), ylab = "Sorted Pvalue", xlab = "%"), plotPar))
#    for(i in 1:length(myObj))
#    {
#        tmpn = attr(obj, "length")
#        obj = myObj[[i]]
#        lines(1:tmpn * 100 / tmpn, obj$qvalue, col = i, lwd = 0.5 * par("lwd"))
#    }
#    if(leg == TRUE) do.call("legend", c(list("topleft", legend = format(sapply(myObj, function(x) attr(x, "FDR"))), lty = 1, col = 1:length(myObj), title = "FDR"), legPar))
#}

### PIE
#' PIE plot BHRejVec class.
#'
#' \code{PIE.BHRejVec} produces pie plots given two BHRejVec objects. THey must be produced
#' with the same samVec objec and have the same FDR values but different test types.
#'
#' @param myDE,myEE Objects of BHRejVec class to be ploted.
#' @param figcol An integer indicating the number of columns of figures plotted.
#' @param graPar A list of parameters adjusting global graphics setting.
#' @param mainPar A list of parameter adjusting the main title.
#' @param innerPar A list of parameters adjusting the inner pie plots setting.
#' @param innerPIE A list of parameters adjusting the pie plots setting.
#' @param legPar A list of parameters adjusting legend.
#' @param xshrink,yshrink Numeric values between 0 and 1 adjusting the height of figures and
#' width.
#' @return NULL.
#' @examples
#' data(GCwPADataA)
#' t1 = samVec(GCwPADataA, selCol = list(1:5, 6:10, 11:15, 16:20, 21:25, 26:30), labels = LETTERS[1:6])
#' EE = BHRejVec(FDR = c(0.1, 0.05, 0.01, 0.005, 0.001), testset = t1)
#' DE = BHRejVec(FDR = c(0.1, 0.05, 0.01, 0.005, 0.001), testset = t1, ifF = TRUE)
#' pdf("test.pdf", width = 12, height = 8)
#' PIE(DE, EE, figcol = 3, mainPar = list(y = 1.5, main = "TITLE", cex = 2.5),
#'  innerPIE = list(radius = 1, col = c("blue", "red", "green")),
#'  innerPar = list(mar = c(1.1, 2.1, 6.1, 0)), legPar = list(x = "bottomright", cex = 1.5))
#' dev.off()
#' @export
PIE.BHRejVec = function(myDE, myEE, figcol = 1, graPar = list(), mainPar = list(), innerPar = list(), innerPIE = list(), legPar = list(), xshrink = 1, yshrink = 0.95)
{
    tmppar = paraMerge(list(fig = par("fig"), new = par("new")), par(graPar))
    on.exit(par(tmppar))
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    do.call("text", paraMerge(list(x = 1, y = 1.5, labels = mainPar$main, font = 2, cex = 2, xpd = NA), mainPar))

    if(is.null(innerPIE$col)) tmpcol = c("lightblue", "yellow", "firebrick")
    else tmpcol = innerPIE$col

    if(is.null(legPar$ifleg) || legPar$ifleg)
    {
        legPar$ifleg = NULL
        do.call("legend", paraMerge(list(x = "topright", fill = tmpcol, legend = c("Other", "EE", "DE"), xpd = NA, bty = "n"), legPar))
    }

    if(!is(myEE, "BHRejVec")) stop("Invalid input object")
    DE = myDE$BHProc
    EE = myEE$BHProc
    if(length(DE) != length(EE)) stop("Lenght of DE and EE data must be the same")

    FDR = myDE$FDR

    pltdes = plotMat(length(FDR), figcol = figcol, xshrink = xshrink, yshrink = yshrink)
    sapply(seq_along(FDR), function(tmpi)
    {
        par(fig = pltdes[[tmpi]], new = TRUE)
        tmpEE = EE[[tmpi]]
        tmpDE = DE[[tmpi]]
        do.call("PIE.BHRej", paraMerge(list(DE = tmpDE, EE = tmpEE, graPar = innerPar), list(piePar = paraMerge(list(col = tmpcol), innerPIE))))
        invisible()
    })
    invisible()
}

