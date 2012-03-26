\name{switch_pdf}
\alias{switch_pdf}
\title{function to calculate the joint pdf for a sample of values from one of a number of pdfs}
\usage{
  switch_pdf(x, pdf, cur.par, RE)
}
\arguments{
  \item{x}{values to be evaluated}

  \item{pdf}{probability density function (pois1, poisson,
  pois1_ln, poisson_ln, normal, multinom)}

  \item{cur.par}{a vector giving parameters for the
  specified distribution; only the first is used for single
  parameter distributions}

  \item{RE}{random effects, if present}
}
\value{
  total log likelihood of points
}
\description{
  function to calculate the joint pdf for a sample of
  values from one of a number of pdfs
}
\author{
  Paul B. Conn
}
\keyword{density}
\keyword{probability}

