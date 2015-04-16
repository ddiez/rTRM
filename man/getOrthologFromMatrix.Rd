\name{getOrthologFromMatrix}
\alias{getOrthologFromMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain gene identifiers for a target organism associated with a list of PWMs.
}
\description{
Obtain gene identifiers for a target organism associated with a list of PWMs.
}
\usage{
getOrthologFromMatrix(filter, organism = "human", dbname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{
vector of matrices to filter results.
}
  \item{organism}{
target organism.
}
  \item{dbname}{
database- usually not need to specify.
}
}

\author{
Diego Diez
}