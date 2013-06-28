\name{getBiogridData}
\alias{getBiogridData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Downloads network data from BioGRID in TAB2 format.}
\description{
This function is used to generate igraph network objects from BioGRID data. It downloads the database into a data.frame object that can be used later with processBiogrid()
}
\usage{
getBiogridData(release)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{release}{
release of BioGRID to download.
}
}
\details{
The release to download must be specified as currently there is no way to download automatically the latests release.
}
\value{
An data.frame object.
}
\author{
Diego Diez
}
