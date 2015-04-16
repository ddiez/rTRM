\name{layout.arc}
\alias{layout.arc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Layouts a graph using arcs.
}
\description{
Generates a layout for graphs that places in the center the target transcription factors, in the sides the enriched transcription factors and in between of them the bridge proteins.
}
\usage{
layout.arc(g, target, query)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
the graph object to layout.
}
  \item{target}{
list of target nodes (typically target transcription factors.)
}
  \item{query}{
list of query nodes (typically enriched transcription factors.)
}
}
\value{
A matrix with the x and y locations of each node in the target graph.
}
\author{
Diego Diez
}
