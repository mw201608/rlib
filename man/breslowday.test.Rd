\name{breslowday.test}
\alias{breslowday.test}
\title{Breslow and Day Test}
\description{
Function to perform the Breslow and Day (1980) test with Tarone correction
}
\usage{
breslowday.test(x)
}
\arguments{
 \item{x}{matrix, a 2 x 2 x K contingency table.}
}
\details{
Uses the equations in Lachin (2000) p. 124-125. Note that the results of the Tarone corrected test do not correspond to the numbers in the Lachin book. Adapted from R programm by Michael Hoehle <http://www-m4.ma.tum.de/pers/hoehle>.
}
\value{
A vector with three values:
\item{X2.HBD}{the Breslow and Day test statistic}
\item{X2.HBDT}{the Breslow and Day test statistic corrected by Tarone}
\item{p}{p value based on the Tarone test statistic using a chi-square distribution with K-1 degree of freedom.}
}
\examples{
##not run
}
