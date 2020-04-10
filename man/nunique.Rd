\name{nunique}
\alias{nunique}
\alias{deEmptyNA}
\title{Count Unique Elements or Remove Missing Values}
\description{
Handle unique elements and missing values
}
\usage{
nunique(x, na.rm=FALSE)
deEmptyNA(x, drop.empty.string=TRUE)
}
\arguments{
 \item{x}{vector.}
 \item{na.rm}{whether missing values should be omitted.}
 \item{drop.empty.string}{whether to drop empty character string.}
}
\details{
Function \code{nunique} counts the unique elements and \code{deEmptyNA} ommits missing values and/or empty character.
}
\value{
Function \code{nunique} returns a number while \code{deEmptyNA} returns a vector.
}
\examples{
##not run###
}
