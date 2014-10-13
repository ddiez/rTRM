\name{writeTRMreport}
\alias{writeTRMreport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Export a table with TRM nodes and associated information.}
\description{
This function generates a data.frame with the nodes in the provided graph and associated annotations.
}
\usage{
writeTRMreport(graph, file, organism, target, query, sort.by = "symbol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
a graph object.
}
  \item{file}{
file name.
}
  \item{organism}{
organisms for the annotations.
}
  \item{target}{
target transcription factor.
}
  \item{query}{
query transcription factors.
}
  \item{sort.by}{
order the columns of the data.frame by (default: "symbol").
}
}
\author{
Diego Diez
}