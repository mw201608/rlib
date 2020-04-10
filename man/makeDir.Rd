\name{makeDir}
\alias{makeDir}
\title{File system utilities}
\description{
Make directory recursively.
}
\usage{
makeDir(path, showWarnings = TRUE, recursive = TRUE, mode = "0777")
}
\arguments{
 \item{path}{a character vector containing a single path name.}
 \item{showWarnings}{logical; should the warnings on failure be shown?}
 \item{recursive}{logical. Should elements of the path other than the last be created? If true, like the Unix command mkdir -p.}
 \item{mode}{the mode to be used on Unix-alikes: it will be coerced by as.octmode. For Sys.chmod it is recycled along paths.}
}
\details{
These functions are wrapper of existing read/write functions.
}
\value{
\code{makeDir} returns invisibly a logical vector indicating if the operation succeeded for each of the files attempted. If the directory already exists, it does nothing.
}
\seealso{\code{\link{dir.create}}}
\examples{
##not run###
}
