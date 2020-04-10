\name{rbind2m}
\alias{rbind2m}
\alias{cbind2m}
\title{Combine two matrices or data.frame}
\description{
Take two matrices or data-frame arguments which may have different rownames or colnames and combine by columns or rows, respectively. 
}
\usage{
rbind2m(x,y)
cbind2m(x,y)
}
\arguments{
 \item{x, y}{matrix or data.frame.}
}
\details{
\code{x} and \code{y} must be the same data type.
The two functions will be re-implemented by the generic rbind and cbind functions in the future.
}
\value{
A matrix or data.frame.
}
\seealso{\code{\link{rbind}}, \code{\link{cbind}}}
\examples{
##not run###
}
