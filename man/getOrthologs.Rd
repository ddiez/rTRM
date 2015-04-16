\name{getOrthologs}
\alias{getOrthologs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Obtain the mapping to Entrez Gene identifiers in the given organism.
}
\description{
Obtain the mapping to Entrez Gene identifiers in the given organism.
}
\usage{
getOrthologs(filter, organism, dbname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{entrezgene identifiers for the original mapping (PWM to gene). These can belong to diverse species and correspond to the "entrezgene" column obtained with getMaps() function.}
  \item{organism}{target organisms, currently supported "human" and "mouse"}
  \item{dbname}{
%%     ~~Describe \code{dbname} here~~
}
}
\details{
If organism is not specified the entire table of orthologs (with all supported species) is returned.
}
\value{
A data.frame object with ortholog information.
}
\author{
Diego Diez
}
\examples{
getOrthologs(organism = "human")
}