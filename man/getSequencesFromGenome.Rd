\name{getSequencesFromGenome}
\alias{getSequencesFromGenome}

\title{
Retrieves a set of sequences from a BSgenome object and optionally appends a label to each sequence id.
}
\description{
This is just a wrapper to getSeq() in package Biostrings that facilitates adding a label to each sequence.}
\usage{
getSequencesFromGenome(BED, genome, append.id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BED}{
file with peak locations in BED format.}
  \item{genome}{
a BSgenome object (e.g. Mmusculus)
}
  \item{append.id}{
optional label to append to each sequence id.}
}
\author{
Diego Diez
}