\name{amvce}
\alias{amvce}
\title{Animal model variance component estimation}
\description{
Linear mixded model (or quantitative genetic animal model) variance component estimation using REML
}
\usage{
amvce(y, X, A, tolerance=1.0e-3, maxIter=100, verbose=FALSE)
}
\arguments{
 \item{y}{a vector of numeric values of dependent variable.}
 \item{X}{a matrix of fixed effects.}
 \item{A}{a square matrix of random effect correlation structure.}
 \item{tolerance}{optional, a relative tolerance to test converging.}
 \item{maxIter}{optional, maximum number of iteration.}
 \item{verbose}{optional, whether to print progress.}
}
\details{
This function implements mixed effect model fitting by REstricted Maximum Likelihood (REML) using the Expectation Maximization (EM) algorithm. For technical details of the analysis model, please refer to the reference of Yang et al 2014.
}
\value{
A list with elements:
\item{coefficients}{estimates of fixed effect coefficients}
\item{apVb}{estimate of variance of fixed effect coefficients}
\item{VC}{a vector of variance component estimates}
\item{loglik}{log-likelihood}
}
\references{
Yang S, Liu Y, Jiang N, Chen J, Leach L, Luo Z and Wang M. Genome-wide eQTLs and heritability for gene expression traits in unrelated individuals. BMC Genomics 2014, 15:13
}
\author{
Minghui Wang <m.h.wang@live.com>
}
\examples{
\dontrun{
#Create example data
set.seed(1234567)
n=200
M1=matrix(rnorm(n*300), nrow=n)
K1=cov(t(M1))
A=K1/mean(diag(K1))
covY=2*A+1*diag(n)
y=10+crossprod(chol(covY),rnorm(n))
X=matrix(rep(1, n),ncol=1)
#Fit linear mixed model
(fit2=amvce(y,X,A))
}
}
