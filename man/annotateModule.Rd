\name{annotateModule}
\alias{annotateModule}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate a network module with information
}
\description{
Uses information about expression, enrichment and parent PPI network to annotate a subgraph.
}
\usage{
annotateModule(g, enrich, trm, targets, ppi, exprs, tfs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
graph to annotate in igraph format.
}
  \item{enrich}{
list of enriched transcription factors (or motifs).
}
  \item{trm}{
TRM to compare with (to identify bridges).
}
  \item{targets}{
list of target transcription factors (typically those with ChIP-seq data).
}
  \item{ppi}{
parent PPI network (to check membership of nodes).
}
  \item{exprs}{
list of entrezgene ids representing expressed genes.
}
  \item{tfs}{
}
}

\author{
Diego Diez
}
