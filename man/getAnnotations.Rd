\name{getAnnotations}
\alias{getAnnotations}
%%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Obtain the 'pwm' table fromt the database, containing PWM's annotations.
}
\description{
Obtain the 'pwm' table fromt the database, containing PWM's annotations.
}
\usage{
getAnnotations(filter, dbname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{
	one or more PWM ids.
}
  \item{dbname}{
	the location of the database (to load custom datbases).
}
}

\author{
Diego Diez
}
\examples{
ann = getAnnotations()
}