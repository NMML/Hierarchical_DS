\name{get_mod_matrix}
\alias{get_mod_matrix}
\title{function to produce a design matrix given a dataset and user-specified formula object}
\usage{
  get_mod_matrix(Cur.dat, stacked.names, factor.ind,
    Det.formula, Levels)
}
\arguments{
  \item{Cur.dat}{current dataset}

  \item{stacked.names}{column names for current dataset}

  \item{factor.ind}{a list of indicator variables (1 =
  factor/categorical variable, 0 = continuous variable)}

  \item{Det.formula}{a formula object}

  \item{Levels}{A list object giving the number of levels
  for factor variables}
}
\value{
  a design matrix
}
\description{
  function to produce a design matrix given a dataset and
  user-specified formula object
}
\author{
  Paul B. Conn
}
\keyword{matrix}
\keyword{model}

