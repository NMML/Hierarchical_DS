\name{switch_sample_prior}
\alias{switch_sample_prior}
\title{function to sample from hyperpriors of a specified probability density function; note that
initial values for sigma of lognormal random effects are fixed to a small value (0.05) to
prevent numerical errors}
\usage{
  switch_sample_prior(pdf, cur.par)
}
\arguments{
  \item{pdf}{probability density function (pois1, poisson,
  normal, unif.disc, unif.cont)}

  \item{cur.par}{a vector giving parameters for the
  specified distribution; only the first is used for single
  parameter distributions}
}
\value{
  a vector of length n samples from the desired
  distribution
}
\description{
  function to sample from hyperpriors of a specified
  probability density function; note that initial values
  for sigma of lognormal random effects are fixed to a
  small value (0.05) to prevent numerical errors
}
\author{
  Paul B. Conn
}
\keyword{density}
\keyword{probability}

