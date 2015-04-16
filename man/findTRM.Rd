\name{findTRM}
\alias{findTRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Indentifies a TRM associated with a target node and one or more query nodes.
}
\description{
This the main function used to identify TRMs. It takes a graph object and use it to search in the neighborhood of a target node for query nodes that are separated a maximum distance (controlled by max.bridge parameter).
}
\usage{
findTRM(g, target, query, method = "nsa", max.bridge = 1, extended = FALSE, strict = FALSE, type = "igraph")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
the network used to identify TRMs (tipically a PPI network)
}
  \item{target}{
character variable with the name of a target node.
}
  \item{query}{
character vector with the list of query nodes.
}
  \item{method}{
method to use.
}
  \item{max.bridge}{
maximum number of nodes allowed between the target and query nodes.
}
\item{extended}{
whether to allow distance restrictions to include both target and query nodes.
}
\item{strict}{
whether to return a single component (using decompose.graph())
}
\item{type}{
type of graph object to return, either an "igraph" (the default) or a "graphNEL"
}
}
\details{
Currently only "first" and "nsa" methods are available. First is used for tests and returns the first neighborhood of the target node. Method "nsa" implements the TRM finding algorithm.
}
\value{
A network in igraph format or NULL.
}

\author{
Diego Diez
}
\examples{
# load example network.
load(system.file(package = "rTRM", "extra/example.rda"))

# define target and query nodes.
target = "N6"
query = c("N7", "N12", "N28")

# find TRM:
s = findTRM(g, target = target, query = query, method = "nsa", max.bridge = 1)
}
