\name{plotTRMlegend}
\alias{plotTRMlegend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the legend of a TRM with information about the cluster families.
}
\description{
This function just plots a legend with the cluster membership of the provided list of genes. The legend includes de most prominent families of each cluster and there is some name polishing as well.
}
\usage{
plotTRMlegend(x, title = NULL, cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{list of family names or igraph object.}
  \item{title}{title for the legend.}
  \item{cex}{numeric value controlling the size of the legend's text.}
}

\author{
Diego Diez
}
