\name{r2p}
\alias{r2p}
\title{Compute correlation P value}
\description{
Compute P value significance of correlation coefficient
}
\usage{
r2p(r, n, alternative=c("two.sided", "less", "greater"))
}
\arguments{
 \item{r}{a vector of correlation coefficients.}
 \item{n}{a vector of sample sizes for the corresponding correlation coefficients.}
 \item{alternative}{indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}" or "\code{less}". See \code{\link{cor.test}}.}
}
\details{
This function computes p value significance of correlation coefficient using t distribution approximation.
}
\value{
A vector of p values.
}
\seealso{
\code{\link{cor.test}}
}
\examples{
##not run###
}
