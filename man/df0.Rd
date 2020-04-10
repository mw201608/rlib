\name{df0}
\alias{df0}
\title{Data Frame}
\description{
The function df0 creates data frames with/without converting character vectors to factors. 
}
\usage{
df0(..., row.names = NULL, check.rows = FALSE,
	check.names = TRUE, fix.empty.names = TRUE,
	stringsAsFactors = FALSE)
}
\arguments{
 \item{...}{these arguments are of either the form value or tag = value. Component names are created based on the tag (if present) or the deparsed argument itself.}
 \item{row.names}{NULL or a single integer or character string specifying a column to be used as row names, or a character or integer vector giving the row names for the data frame.}
 \item{check.rows}{if TRUE then the rows are checked for consistency of length and names.}
 \item{check.names}{logical. If TRUE then the names of the variables in the data frame are checked to ensure that they are syntactically valid variable names and are not duplicated. If necessary they are adjusted (by make.names) so that they are. This option might not be available for old R version.}
 \item{fix.empty.names}{logical indicating if arguments which are "unnamed" (in the sense of not being formally called as someName = arg) get an automatically constructed name or rather name "". Needs to be set to FALSE even when check.names is false if "" names should be kept.}
 \item{stringsAsFactors}{logical: should character vectors be converted to factors?}
}
\details{
See \code{Details} of function \code{\link{data.frame}}.
}
\value{
A data frame, a matrix-like structure whose columns may be of differing types (numeric, logical, factor and character and so on).
}
\seealso{\code{\link{data.frame}}}
\examples{
##not run###
}
