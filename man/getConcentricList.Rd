\name{getConcentricList}
\alias{getConcentricList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns a list with nodes membership to be used in a graph with a concentric layout
}
\description{
Specify target and enriched motifs and returns a list with circle membership. This information is used by layout.concentric to position the nodes in plots.
}
\usage{
getConcentricList(g, t, e, max.size = 60, order.by = "label")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
graph to layout (extract the nodes).
}
  \item{t}{
list of target nodes (will go in the center).
}
  \item{e}{
list of enriched nodes (will go in the periphery).
}
  \item{max.size}{
maximum number of nodes per layer.
}
  \item{order.by}{ordering attribute for list before split.
}

}
\author{
Diego Diez
}
