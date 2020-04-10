\name{Jonckheere.test}
\alias{Jonckheere.test}
\title{Jonckheere Trend Analaysis}
\description{
Compute the Jonckheere Trend Test
}
\usage{
Jonckheere.test(y,x)
}
\arguments{
 \item{y}{a vector of numeric data values.}
 \item{x}{a vector or a factor specifying the ordered group of each values in y.}
}
\value{
A vector with two values: J (test statistics), Z (normalized J value), and p.value (two-sided p value calculated from normal approximation).
}
\references{
Reference: Hollander and Wolfe. 1999. Nonparametric statistical methods (Second edition). John Wiley & Sons, Inc., New York.
}
\seealso{
\code{\link{MannKendall.test}}
}
\examples{
\dontrun{
set.seed(1234)
#generate data in 7 groups with increasing trend
x=rep(1:7,each=3)
y1=rnorm(1:length(x),mean=x)
#generate data without an obvious trend
y2=rnorm(1:length(x),mean=2)
#plot data
plot(x,y1,xlab='X',ylab='y')
points(x,y2,col=2)
abline(lm(y1~x))
abline(lm(y2~x),col=2)
#test for trend in group ordered 1 to 7
(fit1=Jonckheere.test(y1,x))
(fit2=Jonckheere.test(y2,x))
#test if two trends are significantly different
(z <- as.vector(fit1['Z'] - fit2['Z']))
(p=2*pnorm(z,lower.tail= z < 0))
}
}
