\name{getTFterms}
\alias{getTFterms}
\title{
Get terms associated with a specified TFclass subset.
}
\description{
Returns a vector of names (not ids) with the members of a particular subset in the TFclass database. By default it returns the Class subset.
}
\usage{
getTFterms(subset = "Class", dbname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subset}{a subset in TFclass (default Class).}
  \item{dbname}{SQLite file to use as database.}
}
\author{
Diego Diez
}