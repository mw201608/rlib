\name{labels2colors}
\docType{package}
\alias{labels2colors}
\title{Convert Labels to Colors}
\description{
Drop duplication by selecting one with the largest variation
}
\usage{
labels2colors(labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey",
	commonColorCode = TRUE)
}
\arguments{
 \item{labels}{Vector or matrix of non-negative integer or other (such as character) labels.}
 \item{zeroIsGrey}{If TRUE, labels 0 will be assigned color grey. Otherwise, labels below 1 will trigger an error.}
 \item{colorSeq}{Color sequence corresponding to labels. If not given, a standard sequence will be used.}
 \item{naColor}{Color encodes missing values.}
 \item{commonColorCode}{logical: if labels is a matrix, should each column have its own colors?}
}
\details{
This function is derived from WGCNA:::labels2colors.
}
\value{
A vector or array of character strings of the same length or dimensions as \code{labels}.
}
\examples{
##not run###
}
