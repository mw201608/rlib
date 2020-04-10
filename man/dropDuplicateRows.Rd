\name{dropDuplicateRows}
\docType{package}
\alias{dropDuplicateRows}
\title{Drop Duplication}
\description{
Drop duplication by selecting one with the largest mean or standard deviation
}
\usage{
dropDuplicateRows(x, Rows, keep.row.names=TRUE, By=c('mean','sd'))
}
\arguments{
 \item{x}{numeric matrix with variables in the rows and samples in the columns.}
 \item{Rows}{a character vector of group labels for each row variable of \code{x}.}
 \item{keep.row.names}{logical: whether to keep original row names?}
 \item{By}{which statistics to use as criterion.}
}
\details{
This function treats duplication as a single group and select one with the largest mean or standard deviation from each group.
}
\value{
A matrix.
}
\examples{
##not run###
}
