\name{write.tsv}
\alias{write.tsv}
\alias{read.tsv}
\alias{writeCsv}
\alias{fread2}
\alias{read.zip}
\alias{read.gzip}
\alias{read.rds}
\alias{save.rds}
\title{Input and Output Utilities}
\description{
Wrapper functions to read or write from/to a file or connection.
}
\usage{
write.tsv(x,file="", row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE, verbose=TRUE)
read.tsv(file, as.is=TRUE, header = TRUE, sep = "\t", quote = "\"",
	dec = ".", fill = TRUE, comment.char = "", verbose=TRUE, ...)
writeCsv(x, file="", row.names=FALSE, quote=TRUE, verbose=TRUE)
fread2(file, verbose=TRUE, showProgress=FALSE, ...)
read.zip(file, FUN = read.table, ...)
read.gzip(file, FUN = read.delim, ..., verbose = TRUE)
read.rds(file, refhook = NULL)
save.rds(object, file = "", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
}
\arguments{
 \item{x}{the object to be written, preferably a matrix or data frame.}
 \item{file}{the name of the file which the data are to be read from or written to.}
 \item{row.names}{either a logical value indicating whether the row names of x are to be written along with x, or a character vector of row names to be written.}
 \item{col.names}{either a logical value indicating whether the column names of x are to be written along with x, or a character vector of column names to be written.}
 \item{as.is}{whether to suppress the conversion of character variables (which are not converted to logical, numeric or complex) to factors.}
 \item{header}{a logical value indicating whether the file contains the names of the variables as its first line.}
 \item{sep}{the field separator character.}
 \item{quote}{the set of quoting characters.}
 \item{dec}{the character used in the file for decimal points.}
 \item{fill}{logical. If TRUE then in case the rows have unequal length, blank fields are implicitly added.}
 \item{comment.char}{a character vector of length one containing a single character or an empty string.}
 \item{verbose}{whether to print the name of the file to be read from or write to.}
 \item{showProgress}{whether to displays progress on the console.}
 \item{FUN}{function to process the file handler.}
 \item{object}{R object to serialize.}
 \item{ascii}{See \code{\link{saveRDS}}.}
 \item{version}{See \code{\link{saveRDS}}.}
 \item{compress}{See \code{\link{saveRDS}}.}
 \item{refhook}{a hook function for handling reference objects.}
 \item{...}{additional arguments, maybed used by \code{FUN}.}
}
\details{
These functions are wrapper of existing read/write functions. \code{read.zip} will read a single file from a zip archive. Currently \code{read.zip} will not work when there are multiple files in the zip archive. \code{readGzOnTheFly} will read from a remote or local gzip archive with a single file.
}
\value{
For \code{read.tsv}, a data frame (data.frame) containing a representation of the data in the file. For \code{fread2}, a data.table by default. For \code{read.zip} and \code{read.gzip}, the return value depends on \code{FUN}.
}
\seealso{\code{\link{read.table}}, \code{\link{write.table}}, \code{\link{write.csv}}, \code{\link[data.table]{fread}}, \code{\link{unzip}}, \code{\link{readRDS}}, \code{\link{saveRDS}}}
\examples{
##not run###
}
