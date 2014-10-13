\name{getMatrices}
\alias{getMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Obtain a list of PWMs.
}
\description{
Returns a list of PWMs, by default all the PWMs in the database. Alternativelly, filtered by the ids provided by filter.
}
\usage{
getMatrices(filter, dbname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{
list of PWMs to filter results.}
  \item{dbname}{
%%     ~~Describe \code{dbname} here~~
}
}
\author{
Diego Diez
}
\examples{
pwms = getMatrices()
}
