\name{ratio.test}
\alias{ratio.test}
\title{Ratio Test}
\description{
Compute ratio and its standard error
}
\usage{
ratio.test(x, y, covxy=0)
}
\arguments{
 \item{x}{a vector of numeric values in the first sample.}
 \item{y}{a vector of numeric values in the second sample.}
 \item{covxy}{numerc, covariance between variable x and variable y.}
}
\details{
This function computes estimates of ratio and standard error between the means of two given samples. Confidence interval can be used to test if ratio is significantly different from 1, eg a 95% CI will be (ratio-1.96*se, ratio+1.96*se).
}
\references{
A biologist`s guide to statistical thinking and analysis, by David S. Fay and Ken Gerow (http://www.wormbook.org/chapters/www_statisticalanalysis/statisticalanalysis.html).
}
\value{
A vector with two values: ratio, and se (standard error) of ratio.
}
\examples{
##not run###
(fit=ratio.test(c(100,200,300),c(150,190,250)))
#95% confidence interval
c(lower=fit['ratio']-1.96*fit['se'], upper=fit['ratio']+fit['se'])
}
