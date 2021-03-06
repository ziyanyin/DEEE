% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WV_grpWVClass.R
\name{madWV}
\alias{madWV}
\title{Produce a grpWV object}
\usage{
madWV(SVobj, degree = 4)
}
\arguments{
\item{SVobj}{A samVec object containing data and other information.}

\item{degree}{A integer specifying the degree of the polynomail logistic regression that
applied to fit the data.}
}
\value{
A grpWV object is returned, which is a list containing
\describe{
 \item{fitY}{A list of lists of numerics containing the fitted response values.}
 \item{rankit}{A list of numerics containing the rankit of the means of the data matrix.}
 \item{order}{A list of integers containing the order the rankit.}
 \item{nWV}{An integers indicating the numbers of grpWV.}
 \item{interval}{c(0, 1).}
 \item{Function}{"mad".}
 \item{ifeq}{">".}
 \item{selCol}{A list of vector of integers containing absolute indexes of choosen
 columns.}
 \item{colInd}{A list of vector of integers containing relative indexes of columns.}
 \item{labels}{A vector of characters, integers or short string marking the groups.}
 \item{dataType}{A long string containing useful information of the data.}
}
}
\description{
\code{madWV} defines and produces a grpWV object with FUN = ``mad'', interval = c(0, 1)
and ifeq = ``>''. See also \code{\link{grpWV}}, which produces a grpWV object with no
restriction but is much slower. See \code{\link{plot.grpWV}}
for plotting instructions.
}
\examples{
data(GCwPADataA)
testset = samVec(GCwPADataA, selCol = list(1:5, 11:15, 21:25), labels = c("E", "R", "T"),
 dataType = "Example")
gw1 = madWV(SVobj = testset, degree = 6)
plot(gw1, legPar = list(cex = 1, ncol = 1), mainPar = list(main = "Example"))
}

