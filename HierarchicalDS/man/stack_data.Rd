\name{stack_data}
\alias{stack_data}
\title{function to stack data (going from three dimensional array to a two dimensional array including only "existing" animals}
\usage{
  stack_data(Data, Obs.transect, n.transects,
    stacked.names, factor.ind)
}
\arguments{
  \item{Data}{three-d dataset}

  \item{Obs.transect}{current number of observations of
  animals in each transect (vector)}

  \item{n.transects}{number of transects}

  \item{stacked.names}{column names for new stacked
  dataset}

  \item{factor.ind}{a vector of indicator variables (1 =
  factor/categorical variable, 0 = continuous variable)}
}
\value{
  a stacked dataset
}
\description{
  function to stack data (going from three dimensional
  array to a two dimensional array including only
  "existing" animals
}
\author{
  Paul B. Conn
}
\keyword{data}
\keyword{stack}

