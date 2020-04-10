\name{recoverGeneSymbolFromExcelDate}
\alias{recoverGeneSymbolFromExcelDate}
\title{Recover Gene Symbol}
\description{
Recover gene symbols from excel date format
}
\usage{
recoverGeneSymbolFromExcelDate(genes, species='human')
}
\arguments{
 \item{genes}{a vector of character strings.}
 \item{human}{Default 'human'. If not 'human', recovered symbols will have capital case only for the first letter.}
}
\details{
Microsoft excel will automatically converts any date-like string to date format. Unfortunately, this behavior will mess up some gene symbols eg, DEC1, MARCH1, SEPT1. This function tries to recover the date format into original gene symbols.
}
\value{
A vector.
}
\examples{
##not run###
}
