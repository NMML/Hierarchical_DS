\name{simulate_data}
\alias{simulate_data}
\title{function to simulate distance sampling data from simple model with increasing abundance
intensity, no assumed spatial structure, and point independence.  If no parameters are given, uses internally defined values.}
\usage{
  simulate_data(S, Observers, misID = TRUE, X.site = NULL,
    n.species = 2, Beta.hab = NULL, Beta.det = NULL,
    detect.model = ~Observer + Distance + Group + Species,
    dist.cont = FALSE, n.bins = 5, cor.par = 0.5,
    Grp.par = c(3, 1), Conf.model = NULL, Conf.par = NULL)
}
\arguments{
  \item{S}{number of spatial strata (a single transect is
  placed in each strata and assumed to cover the whole
  strata)}

  \item{Observers}{A (2 x S) matrix giving the observer
  identity for each transect}

  \item{misID}{If TRUE (default), misidentification is
  assumed}

  \item{X.site}{model.matrix for habitat covariates
  (defaults to a linear, quadratic effects of transect # on
  log scale)}

  \item{n.species}{Number of species to simulate (current
  max is 2) (default is 2)}

  \item{Beta.hab}{A (# of species X # parameters) matrix
  giving parameters for habitat-abundance relationship
  (default is linear increase for species 1, quadratic for
  species 2)}

  \item{detect.model}{A formula for the detection model.
  Default is ~Observer+Distance+Group+Species (formula
  should consist of these key words)}

  \item{Beta.det}{A vector giving parameters for the
  detection function; # of parameters must match
  model.matrix()!}

  \item{dist.cont}{If TRUE, uses continuous distances on
  (0,1).  If FALSE, uses discrete distance classes}

  \item{n.bins}{If dist.cont=FALSE, how many bins to use
  for distances.  Default is 5.}

  \item{cor.par}{Correlation at maximum distance.  Default
  is 0.5}

  \item{Grp.par}{A vector with an entry for each species
  giving parameters for group size (assumed zero-truncated
  Poisson). Default is 3 and 1, corresponding to mean group
  sizes of 4 and 2 for each species}

  \item{Conf.model}{If misID=TRUE, this gives a list vector
  with formulas for the confusion matrix. Formulas are
  provided for "getting it right" and "getting an unknown".
  Default is c(~Observer+Group+Distance+Species,
  ~Observer). The same models are assumed for each
  species.}

  \item{Conf.par}{If misID=TRUE, ths is a list vector with
  parameters for each confusion model.  Default is
  list(c(2,.2,.4,.3,-.1,-.2,-.4,-.8,.5),c(-1,-.2,.2))}
}
\value{
  a distance sampling dataset
}
\description{
  function to simulate distance sampling data from simple
  model with increasing abundance intensity, no assumed
  spatial structure, and point independence.  If no
  parameters are given, uses internally defined values.
}
\author{
  Paul B. Conn
}
\keyword{distance}
\keyword{sampling,}
\keyword{simulation}

