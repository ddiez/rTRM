\name{getTFclassFromEntrezgene}
\alias{getTFclassFromEntrezgene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Applies getTFclass sequentially to a vector of entrezgene identifiers.
}
\description{
Applies getTFclass sequentially to a vector of entrezgene identifiers.
}
\usage{
getTFclassFromEntrezgene(x, subset = "Class", tfclass, dbname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of entrezgene identifiers.}
  \item{subset}{
level in the ontology (subset in TFclass terminology. By default "Class")
}
  \item{tfclass}{
data.frame with tfclass data to pass to the recursive function.
}
  \item{dbname}{
  SQLite file to use as database.
  }
}

\author{
Diego Diez
}