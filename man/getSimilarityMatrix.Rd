\name{getSimilarityMatrix}
\alias{getSimilarityMatrix}

\title{
Compute similarity matrix of list of graphs.
}
\description{
This function computes pair-wise similarity based on common nodes (default) or edges between the graphs passed as a list.
}
\usage{
getSimilarityMatrix(g_list, type = "nodes")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g_list}{list of graph objects.}
  \item{type}{type of similarity, either node (default) or edge.}
}
\author{
Diego Diez
}