\docType{methods}
\name{dcor}
\alias{dcor}
\alias{ddcor}
\alias{ddcor,matrix,matrix-method}
\alias{ddcor,numeric,numeric-method}
\title{Differential Correlation}
\description{
Compare two correlation coefficients using Fisher's Z transformation
}
\usage{
dcor(r1, r2, n1, n2, two.tailed=FALSE)
ddcor(cor1, cor2, n1, n2, two.tailed=FALSE, adjust='BH', FDR=0.05, ...)
\S4method{ddcor}{matrix,matrix}(cor1, cor2, n1, n2, two.tailed=FALSE, 
adjust='BH', FDR=0.05, selection=c('whole','upper.tri','lower.tri'), diag=FALSE, ...)
\S4method{ddcor}{numeric,numeric}(cor1, cor2, n1, n2, two.tailed=FALSE, 
adjust='BH', FDR=0.05, ...)
}
\arguments{
 \item{r1}{numeric, a vector of correlation coefficients.}
 \item{r2}{numeric, a vector of correlation coefficients.}
 \item{n1}{sample size in computing \code{r1} or \code{cor1}.}
 \item{n2}{sample size in computing \code{r2} or \code{cor2}.}
 \item{two.tailed}{logical: whether to give two.tailed test p value.}
 \item{cor1}{numeric, a vector or a matrix of correlation coefficients.}
 \item{cor2}{numeric, a vector or a matrix of correlation coefficients.}
 \item{adjust}{character, method to compute FDR. See \code{\link{p.adjust}}}
 \item{FDR}{numeric, false discovery rate threshold.}
 \item{selection}{character, indicating whether to use the whole, upper triangle or lower triangle of the input matrices \code{cor1} and \code{cor2}.}
 \item{diag}{logical, whether to include diagonal when \code{selection} is upper.tri or lower.tri.}
 \item{...}{preserved for additional parameters.}
}
\details{
Function \code{dcor} performs differential correlation test using Fisher's Z transformation. See reference http://www.sciencedirect.com/science/article/pii/S0378111912014497.

Function \code{ddcor} implements DGCA (reference) with two matrices of correlation coefficients.
}
\value{
\code{dcor} returns a vector or a list with two entries: z (z score) and p (p value).
\code{ddcor} returns a data.frame.
}
\examples{
##not run###
}
