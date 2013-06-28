\name{plotGraph}
\alias{plotGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot an graph in igraph format.
}
\description{
This function plots graphs of the class igraph.
}
\usage{
plotGraph(g, layout = layout.fruchterman.reingold, mar = .5, vertex.pch = 21, vertex.cex, vertex.col, vertex.bg, vertex.lwd, edge.col, edge.lwd, edge.lty, label = TRUE, label.col, label.cex, label.pos = NULL, label.offset = 1.5, adjust.label.col = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{a network object.}
  \item{layout}{graph layout, either a function or the output of a layout function.}
  \item{mar}{plot margin.}
  
  \item{vertex.pch}{node size.}
  \item{vertex.cex}{node size.}
  \item{vertex.col}{node line color.}
  \item{vertex.bg}{node background color.}
  \item{vertex.lwd}{node line width.}
  
  \item{edge.col}{edge color.}
  \item{edge.lwd}{edge line width.}
  \item{edge.lty}{edge line type.}
  
  \item{label}{logical; whether to plot labels.}
  \item{label.col}{label color.}
  \item{label.cex}{label expansion.}
  \item{label.pos}{label position.}
  \item{label.offset}{label offset.}
  
  \item{adjust.label.col}{
whether to automatically adjust label color depending on the luminance of the node's color/s.

}
}

\author{
Diego Diez
}