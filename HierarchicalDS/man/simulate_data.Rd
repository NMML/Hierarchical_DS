\name{simulate_data}
\alias{simulate_data}
\title{function to simulate double observer spatial distance sampling data subject to possible zero inflation and species misidentification}
\usage{
  simulate_data(S, Observers, ZIP = TRUE, misID = TRUE,
    tau.pois = 15, tau.bern = 20)
}
\arguments{
  \item{S}{number of spatial strata (a single transect is
  placed in each strata and assumed to cover the whole
  strata)}

  \item{Observers}{A (2 x S) matrix giving the observer
  identity for each transect}

  \item{ZIP}{If TRUE, simulate abundance using a
  zero-inflated Poisson model}

  \item{misID}{If TRUE, assume species misidentification}

  \item{tau.pois}{Precision of the ICAR process for Poisson
  abundance}

  \item{tau.bern}{Precision of the ICAR process for zero
  inflation}
}
\value{
  a distance sampling dataset
}
\description{
  function to simulate double observer spatial distance
  sampling data subject to possible zero inflation and
  species misidentification
}
\author{
  Paul B. Conn
}
\keyword{distance}
\keyword{sampling,}
\keyword{simulation}

