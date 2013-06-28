\name{processBiogrid}
\alias{processBiogrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process a data.frame with BioGRID data into a network for a target organism
}
\description{
Process a data.frame with BioGRID data into a network for a target organism.
}
\usage{
processBiogrid(dblist, org = "human", simplify = TRUE, type = "physical", mimic.old = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dblist}{data.frame containing the BioGRID data.}
  \item{org}{
target organism (default: "human")
}
  \item{simplify}{
whether to eliminate redundant edges (default TRUE)
}
  \item{type}{
type of interaction (physical or genetic) to include (default: "physical")
}
  \item{mimic.old}{mimic old behavior of processBiogrid() when interactions for multiple species could be retrieved. Used only for testing.
}
}
\value{
An igraph object.
}
\author{
Diego Diez
}
