\name{annotateFreq}
\alias{annotateFreq}
\title{
Annotate a graph with frequency of nodes/edges in other graphs.
}
\description{
Returns an annotated graph with node size and edge width proportional at the number of occurrences of nodes/edges in a suplied list of graphs.
}
\usage{
annotateFreq(g, graph_list)
}
\arguments{
  \item{g}{target graph to annotate.}
  \item{graph_list}{list of graph to extract information from.}
}
\details{
Commonly graph_list refers to a list of predicted TRMs (with findTRM) and g is the combined TRM. This function annotates the nodes/edges in g to known their frequency in the original list of graphs.
}
\author{
Diego Diez
}