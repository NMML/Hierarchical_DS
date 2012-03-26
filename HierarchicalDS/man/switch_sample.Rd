\name{switch_sample}
\alias{switch_sample}
\title{function to sample from a specified probability density function}
\usage{
  switch_sample(n, pdf, cur.par, RE)
}
\arguments{
  \item{n}{number of samples desired}

  \item{pdf}{probability density function (pois1, poisson,
  normal, unif.disc, unif.cont)}

  \item{cur.par}{a vector giving parameters for the
  specified distribution; only the first is used for single
  parameter distributions}

  \item{RE}{random effects, if present}
}
\value{
  a vector of length n samples from the desired
  distribution
}
\description{
  function to sample from a specified probability density
  function
}
\author{
  Paul B. Conn
}
\keyword{density}
\keyword{probability}

