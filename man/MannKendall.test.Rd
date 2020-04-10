\name{MannKendall.test}
\alias{MannKendall.test}
\title{Mann Kendall Test}
\description{
Compute the Mann-Kendall Trend Test
}
\usage{
MannKendall.test(x)
}
\arguments{
 \item{x}{a vector of time ordered data.}
}
\details{
The Mann-Kendall test is a nonparametric test for monotonic trends, such as concentrations that are either consistently increasing or decreasing over time. Note that P value is obtained from normal approximation.
}
\value{
A vector with two values: tau, p.value
}
\seealso{
\code{\link{Jonckheere.test}}
}
\examples{
##not run###
}
