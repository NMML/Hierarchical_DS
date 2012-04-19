\name{mcmc_ds}
\alias{mcmc_ds}
\title{Function for MCMC analysis}
\usage{
  mcmc_ds(Par, Data, cur.iter, adapt, Control, DM.hab,
    DM.det, Q, Prior.pars, Meta)
}
\arguments{
  \item{Par}{A list comprised of the following parameters:
  "det": a vector giving the current iteration's linear
  model parameters for the detection model; "hab": a vector
  giving the current iteration's linear model parameters
  for abundance intensity; "cor": a correlation parameter
  for detections that's an increasing function of distance
  (correlation at the maximum distance); "Nu": a vector
  giving the log of the abundance intensity for each
  strata; "G": a vector giving the number of groups of
  animals in each strata; "N": a vector giving the number
  of animals in each strata "Cov.par": an (n.species X n X
  n.ind.cov) array holding parameters of individual
  covariate distributions.}

  \item{Data}{A four dimensional array; the first dimension
  gives species, the second gives the transect, the third
  indexes a (possible) observation, and the fourth
  dimension gives observations and covariates associated
  with a given animal.  These final columns are: Observer
  ID,Y(observation=0/1),Observed species,Obs
  covariates,Distance,Ind covariates}

  \item{cur.iter}{Number of iterations to run}

  \item{adapt}{If adapt==TRUE, run MCMC in adapt mode,
  optimizing MCMC proposal distributions prior to primary
  MCMC}

  \item{Control}{A list object including the following
  slots: "iter": number of MCMC iterations; "burnin":
  number of MCMC burnin iterations; "thin": if specified,
  how many iterations to skip between recorded posterior
  samples; "adapt": if adapt==TRUE, this gives the number
  of additional MCMC iterations should be performed to
  adapt MCMC proposals to optimal ranges prior to final
  MCMC run; "MH.cor": Metropolis-hastings tuning parameter
  for updating the correlation parameter (if
  Meta$point.ind==TRUE); "MH.nu": MH tuning parameter for
  Nu parameters (Langevin-Hastings multivariate update);
  "MH.beta": A matrix of tuning parameters for betas of the
  abundance process (nrows=number of species, ncol = max
  number of columns of habitat DM); "RJ.N": A vector giving
  the maximum number of additions and deletions proposed in
  an iteration of the RJMCMC algorithm for each transect}

  \item{DM.hab}{A design matrix for the log of abundance
  intensity}

  \item{DM.det}{A design matrix for the probit of detection
  probability}

  \item{Q}{An inverse precision matrix for the spatial ICAR
  process}

  \item{Prior.pars}{A list object giving parameters of
  prior distribution.  Includes the following slots
  "a.eta": alpha parameter for prior precision of spatial
  process (assumed Gamma(a.eta,b.eta)) "b.eta": beta
  parameter for prior precision of spatial process (assumed
  Gamma(a.eta,b.eta)) "a.nu": alpha parameter for prior
  precision of overdispersion process (assumed
  Gamma(a.nu,b.nu)) "b.nu": beta parameter for prior
  precision of overdispersion process (assumed
  Gamma(a.nu,b.nu)) "beta.sd": standard deviation for
  regression coefficients (assumed Normal(0,beta.sd^2)
  "misID.mu": a list vector, each entry gives normal prior
  means for misID regression coefficients for the
  corresponding model in Meta$misID.mat (can be set to null
  if no misID) "misID.sd": a list vector, each entry gives
  normal prior sd for misID regression coefficients for the
  corresponding model in Meta$misID.mat (can be set to null
  if no misID)}

  \item{Meta}{A list object giving a number of other
  features of the dataset, including: "n.transects" Number
  of transects "n.species" Number of species "S" Number of
  strata cells "spat.ind" Indicator for spatial dependence
  "Area.hab" Vector giving relative area covered by each
  strata "Area.trans" Vector giving fraction of area of
  relevant strata covered by each transect "Adj" Adjacency
  matrix giving connectivity of spatial grid cells
  "Mapping" Vector mapping each transect into a parent
  strata "Covered.area" Vector giving the fraction of each
  strata covered by transects "n.Observers" Vector giving
  the number of observers that operated on each transect
  "M" Vector giving maximum number of groups in each
  transect "stacked.names" Character vector giving column
  names for the dataset "factor.ind" Indicator vector
  specifying whether data columns are factors (1) or
  continuous (0) "Det.formula" a formula object specifying
  the model for the detection process "Levels" a list
  object, where slot names are comprised of detection model
  names; each slot gives total # of levels in the combined
  dataset "i.binned" indicator for whether distances are
  recorded in bins (1) or are continuous (0) "dist.pl"
  gives the column in Data where distances are located
  "G.transect" vector holding current number of groups of
  animals present in area covered by each transect
  "N.transect" vector holding current number of animals
  present in covered area by each transect "grps" indicator
  for whether observations are for groups rather than
  individuals "n.bins" number of distance bins (provided
  i.binned=1) "Bin.length" vector giving relative size of
  distance bins "n.ind.cov" Number of individual covariates
  (distance is not included in this total, but group size
  is) "Cov.prior.pdf" character vector giving the
  probability density function associated with each
  individual covariate (type ? hierarchical_DS for more
  info) "Cov.prior.parms" An (n.species X n X n.ind.cov)
  array providing "pseudo-prior" parameters for individual
  covarate distributions (only the first row used if a
  signle parameter distribution) "Cov.prior.fixed"
  indicator vector for whether parameters of each covariate
  distribution should be fixed within estimation routine
  "Cov.prior.n" (#species X #covariates) Matrix giving
  number of parameters in each covariate pdf "point.ind"
  Indicator for whether point independence assumed (if no,
  then no correlation modeled b/w multiple observers as
  function of distance) "fix.tau.nu" Indicator for whether
  tau.nu should be fixed (1) or estimated(0) "srr"
  Indicator for whether a spatially restricted regression
  model should be employed (1) or not (0) "srr.tol"
  Threshold eigenvalue level for SRR; only eigenvectors
  with higher eigenvalues than srr.tol are included in SRR
  formulation "misID" If TRUE, misidentification of species
  is modeled "misID.mat" With true state on rows and
  assigned state on column, each positive entry provides an
  index to misID.models (i.e. what model to assume on
  multinomial logit space); a 0 indicates an impossible
  assigment; a negative number designates which column is
  to be obtained via subtraction "misID.models" A formula
  vector providing linar model-type formulas for each
  positive value of misID.mat.  If the same model is used
  in multiple columns it is assumed that all fixed effects
  (except the intercept) are shared "N.par.misID" A vector
  specifying the number of parameters needed for each misID
  model "N.hab.par" A vector specifying the number of
  parameters needed for each species' habitat model}
}
\value{
  returns a list with the following slots: "MCMC": An
  'mcmc' object (see 'coda' R package) containing posterior
  samples; "Accept": A list object indicating the number of
  proposals that were accepted for parameters updated via
  Metropolis- or Langevin-Hastings algorithms; "Control": A
  list object giving MCMC tuning parameters (which are
  updated if the 'adapt' alorithm is used) "Obs.N": Records
  latent abundance in each transect; dimension is
  (n.species X # samples X # transects) "Pred.N": Posterior
  predictive distribution for abundance in each transect;
  obtained by sampling a Poisson distribution given current
  parameter values "Post": Holds posterior samples for
  strata specific group sizes ("Post$G") and abundance
  ("Post$N")
}
\description{
  Function for MCMC analysis
}
\author{
  Paul B. Conn
}
\keyword{areal,}
\keyword{augmentation,}
\keyword{data}
\keyword{distance}
\keyword{jump}
\keyword{mcmc,}
\keyword{reversible}
\keyword{sampling,}

