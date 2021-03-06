% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RandMat.R
\name{RandMat}
\alias{RandMat}
\title{RandMat}
\usage{
RandMat(dim1, dim2)
}
\arguments{
\item{dim1}{number of rows}

\item{dim2}{number of columns}
}
\value{
The random matrix (dim1xdim2) with only one nonzero element per row.
}
\description{
Generates a random binary and row stochastic matrix.
}
\examples{
RandMat(3, 5)
}
