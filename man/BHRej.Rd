% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDR_BHRejClass.R
\name{BHRej}
\alias{BHRej}
\title{Produce a BHRej object}
\usage{
BHRej(FDR = 0.1, pvector)
}
\arguments{
\item{FDR}{A numeric value between 0 and 1 for FDR control.}

\item{pvector}{A vector of numerics containing pvalues.}
}
\value{
A BHRej object is return, which is list containing:
\describe{
 \item{rejList}{A vector of integers containing the relative indexes of rejected genes.}
 \item{pVal}{A vector of numerics containing the tested p-values.}
 \item{FDR}{A numeric value between 0 and 1.}
 \item{nTest}{An integer indicating the number of p-values.}
}
}
\description{
\code{BHRej} implements the Bonjamini-Hocheberg procedure given a vector of p-values.
Heap sort algorithm is used to improve efficiency. See also \code{\link{BHRejVec}}.
See \code{\link{PIE.BHRej}} for plotting instructions.
}
\examples{
DE = BHRej(FDR = 0.1, runif(1000, 0, 0.1))
EE = BHRej(FDR = 0.1, runif(1000, 0, 0.1))
PIE(DE, EE, piePar = list(radius = 0.9, main = "What"))
}

