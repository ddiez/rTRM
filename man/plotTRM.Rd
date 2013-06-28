\name{plotTRM}
\alias{plotTRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot an annotated TRM.
}
\description{
This function plots the output findTRM() after it has been annotated with cluster information with annotateTRM() function. Cluster membership is plotted using a pie plot.
}
\usage{
plotTRM(g, layout = layout.fruchterman.reingold, mar = .5, vertex.col, vertex.cex, vertex.lwd, edge.col, edge.lwd, edge.lty, label = TRUE, label.cex, label.col, label.pos = NULL, label.offset = 1.5, adjust.label.col = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{a network object with cluster information (attribute piecolor).}
  \item{layout}{graph layout, either a function or the output of a layout function.}
  \item{mar}{plot margin.}
  
  \item{vertex.col}{node color.}
  \item{vertex.cex}{node size.}
  \item{vertex.lwd}{node border line width.}
  
  \item{edge.col}{edge color.}
  \item{edge.lwd}{edge line width.}
  \item{edge.lty}{edge line type.}
  
  \item{label}{logical; whether to plot labels.}
  \item{label.cex}{label expansion.}
  \item{label.col}{label color.}
  \item{label.pos}{label position.}
  \item{label.offset}{label offset.}
  \item{adjust.label.col}{whether to automatically adjust label color depending on the luminance of the node's color.}
}

\author{
Diego Diez
}
