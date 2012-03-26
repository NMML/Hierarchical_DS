\name{generate_inits}
\alias{generate_inits}
\title{generate initial values for MCMC chain if not already specified by user}
\usage{
  generate_inits(DM.hab, DM.det, G.transect, Area.trans,
    Area.hab, Mapping, point.ind, spat.ind, grp.mean)
}
\arguments{
  \item{DM.hab}{design matrix for habitat model}

  \item{DM.det}{design matrix for detection model}

  \item{G.transect}{a vector of the number of groups of
  animals in area covered by each transect}

  \item{Area.trans}{a vector giving the proportion of a
  strata covered by each transect}

  \item{Area.hab}{a vector of the relative areas of each
  strata}

  \item{Mapping}{a vector mapping each transect to it's
  associated strata}

  \item{point.ind}{is point independence assumed
  (TRUE/FALSE)}

  \item{spat.ind}{is spatial independence assumed?
  (TRUE/FALSE)}

  \item{grp.mean}{pois1 parameter for group size}
}
\value{
  a list of initial parameter values
}
\description{
  generate initial values for MCMC chain if not already
  specified by user
}
\author{
  Paul B. Conn
}
\keyword{initial}
\keyword{mcmc}
\keyword{values,}

