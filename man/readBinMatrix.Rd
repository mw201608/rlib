\name{writeBinMatrix}
\alias{writeBinMatrix}
\alias{readBinMatrix}
\title{Matrix Binary IO}
\description{
Read and write matrix data from/to a binary file
}
\usage{
writeBinMatrix(x,filename)
readBinMatrix(filename)
}
\arguments{
 \item{x}{a numeric matrix.}
 \item{filename}{character string, file name}
}
\details{
Please be warned that row names and column names will be ignored when writting matrix into a binary file.
}
\seealso{\code{\link{writeBin}},\code{\link{readBin}}}
\examples{
##not run###
}
