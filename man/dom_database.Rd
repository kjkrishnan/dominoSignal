% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dom_database}
\alias{dom_database}
\title{Access database}
\usage{
dom_database(dom, name_only = TRUE)
}
\arguments{
\item{dom}{A domino object that has been created}

\item{name_only}{A boolean for whether to return only the name of the database used
or the entire database that is stored. Default TRUE.}
}
\value{
A vector of unique databases used in building the domino object OR
a data frame that includes the database information used in the domino object creation
}
\description{
A function to pull database information from a domino object
}
\examples{
database_name <- dom_database(dominoSignal:::pbmc_dom_built_tiny)
full_database <- dom_database(dominoSignal:::pbmc_dom_built_tiny, name_only = FALSE)

}
