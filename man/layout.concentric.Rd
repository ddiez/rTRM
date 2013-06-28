\name{layout.concentric}
\alias{layout.concentric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a concentric layout for graphs}
\description{
Generates a matrix with x,y coordinates for each node in a target graph, which layouts the nodes using concentric circles.}
\usage{
layout.concentric(g, concentric = NULL, radius = NULL, order.by)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{graph (igraph) to layout.
}
  \item{concentric}{
list with the components of each layer.
}
  \item{radius}{
radious of each layer.
}
  \item{order.by}{
graph attributes to order nodes by.
}
}

\author{
Diego Diez
}