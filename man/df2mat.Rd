\name{df2mat}
\alias{df2mat}
\title{Data Frame to Matrix}
\description{
The function df2mat convert data.frame to matrix.
}
\usage{
df2mat(x,val,columnFactors,rowFactors,default=NA)
}
\arguments{
 \item{x}{a data.frame.}
 \item{val}{an integer or character string specifying the column in \code{x} holding the actual data values.}
 \item{columnFactors}{a vector of integers or character strings specifying the variables in code{x} from which the combinations form the column names of the output matrix.}
 \item{rowFactors}{a vector of integers or character strings specifying the variables in \code{x} from which the combinations form the row names of the output matrix.}
 \item{default}{the default value of the output matrix when it is missing in \code{x}.}
}
\details{
The column variables specified by \code{columnFactors} and \code{columnFactors} must be factors.
}
\value{
A matrix.
}
\seealso{\code{\link{data.frame}}}
\examples{
##not run###
}
