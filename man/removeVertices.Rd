\name{removeVertices}
\alias{removeVertices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove nodes from a graph and returns the largest component
}
\description{
Remove nodes from a graph and returns the largest component
}
\usage{
removeVertices(g, filter, keep.hanging = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
graph to remove nodes.
}
  \item{filter}{
%%     ~~Describe \code{filter} here~~
}
  \item{keep.hanging}{
(logical) whether to return the largest component or not.
}
}

\author{
Diego Diez
}