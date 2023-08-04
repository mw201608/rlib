\name{writeBinMatrix}
\alias{writeBinMatrix}
\alias{readBinMatrix}
\alias{readSingleColumnFromBinMatrix}
\alias{readSingleColumnFromRemoteBinMatrix}
\title{Matrix Binary IO}
\description{
Read and write 2-D matrix data from/to a binary file
}
\usage{
writeBinMatrix(x, filename, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = "little", useBytes = TRUE )
readBinMatrix(filename, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = "little")
readSingleColumnFromBinMatrix(j, id, filename, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = "little", 
    verbose = FALSE)
readSingleColumnFromRemoteBinMatrix(j, id, url, dimBytes = 8, charBytes = 20, dataBytes = 8, endian = 'little',
    dims = NULL, verbose = FALSE)
}
\arguments{
 \item{x}{a numeric matrix.}
 \item{filename}{character string, file name.}
 \item{dimBytes}{integer. The number of bytes for dimension size.}
 \item{charBytes}{integer. The max number of bytes per element for row and column names.}
 \item{dataBytes}{integer. The number of bytes per element for the data.}
 \item{endian}{\code{big} or \code{little}. See \code{\link{writeBin}}.}
 \item{useBytes}{See \code{\link{writeBin}}.}
 \item{j}{integer, a column index.}
 \item{id}{character string, a column id. Ignored when \code{j} is present.}
 \item{verbose}{Whether to print reading information.}
 \item{url}{character string, URL of a remote file location.}
 \item{dims}{vector, dimensions of the data matrix.}
}
\details{
\code{writeBinMatrix} writes a 2-D numeric matrix into a binary file, which can be read back to R by \code{readBinMatrix}, \code{readSingleColumnFromBinMatrix}, or \code{readSingleColumnFromRemoteBinMatrix} with matching arguments.

To write matrix into a binary file, \code{writeBinMatrix} first writes two integers referring to the number of rows and the number of columns of the matrix respectively.
Next, if \code{charBytes} is a positive value, \code{writeBinMatrix} writes the column names and then row names. Lastly, \code{writeBinMatrix} writes matrix data in a column-major mode.

\code{readSingleColumnFromBinMatrix} is a wraper of \code{readBinMatrix}. When both \code{j} and \code{id} are missing, it calls \code{readBinMatrix}. When either \code{j} or \code{id} is not missing, it reads the \code{j}-th column or the column whose name is \code{id}.
}
\value{
No return value for \code{writeBinMatrix}.

\code{readBinMatrix} returns a matrix by default. However, when \code{dataBytes} is set to 0, it returns either the dimension of the data matrix in the case of \code{charBytes} = 0, or otherwise dimnames of the data matrix stored in the input file. 

\code{readSingleColumnFromBinMatrix} returns a vector if \code{j} or \code{id} is present, otherwise it returns the value from \code{readBinMatrix}.

\code{readSingleColumnFromRemoteBinMatrix} returns a vector.
}
\seealso{\code{\link{writeBin}}, \code{\link{readBin}}}
\examples{
##not run###
}
