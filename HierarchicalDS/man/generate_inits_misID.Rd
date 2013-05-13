\name{generate_inits_misID}
\alias{generate_inits_misID}
\title{generate initial values for misID model if not already specified by user}
\usage{
  generate_inits_misID(DM.hab.pois, DM.hab.bern, DM.det,
    N.hab.pois.par, N.hab.bern.par, G.transect, Area.trans,
    Area.hab, Mapping, point.ind, spat.ind, grp.mean,
    misID, misID.mat, N.par.misID)
}
\arguments{
  \item{DM.hab.pois}{a list of design matrices for the
  Poisson habitat model (elements are named sp1,sp2, etc.)}

  \item{DM.hab.bern}{If a hurdle model, a list of design
  matrices for the Bernoulli habitat model (elements are
  named sp1,sp2, etc.) (NULL if not hurdle)}

  \item{DM.det}{design matrix for detection model}

  \item{N.hab.pois.par}{vector giving number of parameters
  in the Poisson habitat model for each species}

  \item{N.hab.bern.par}{vector giving number of parameters
  in the Bernoulli habitat model for each species (NULL if
  not hurdle)}

  \item{G.transect}{a matrix of the number of groups of
  animals in area covered by each transect; each row gives
  a separate species}

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

  \item{grp.mean}{a vector giving the pois1 parameter for
  group size (one entry for each species)}

  \item{misID}{if TRUE, indicates that misidentification is
  incorporated into modeling}

  \item{misID.mat}{a matrix specifying which elements of
  the misID matrix are linked to model equations}

  \item{N.par.misID}{a vector giving the number of
  parameters for each misID model (in multinomial logit
  space)}
}
\value{
  a list of initial parameter values
}
\description{
  generate initial values for misID model if not already
  specified by user
}
\author{
  Paul B. Conn
}
\keyword{initial}
\keyword{mcmc}
\keyword{values,}

