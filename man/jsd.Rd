\name{JSD}
\alias{JSD}
\title{Jensen-Shannon divergence}
\description{
Divergence measures based on the Shannon entropy
}
\usage{
JSD(mat, pseudocount=0.000001)
}
\arguments{
 \item{mat}{matrix, a contingency table.}
 \item{pseudocount}{a small number added to zero entry.}
}
\value{
Returns an object of class "\code{\link{dist}}". 
}
\seealso{
\code{\link{dist}}
}
\examples{
##not run
}
