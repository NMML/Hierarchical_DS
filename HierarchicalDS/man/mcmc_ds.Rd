\name{mcmc_ds}
\alias{mcmc_ds}
\title{Function for MCMC analysis}
\usage{
  mcmc_ds(Par, Data, cur.iter, adapt, Control, DM.hab.pois,
    DM.hab.bern = NULL, DM.det, Q, Prior.pars, Meta)
}
\arguments{
  \item{Par}{A list comprised of the following parameters:
  "det": a vector giving the current iteration's linear
  model parameters for the detection model; "hab.pois": a
  vector giving the current iteration's linear model
  parameters for Poisson abundance intensity; each row
  gives parameters for a particular species "hab.bern": a
  vector giving the current iteration's linear model
  parameters for Bernoulli part of ZIP model for abundance
  (if Meta$ZIP=TRUE) "cor": a correlation parameter for
  detections that's an increasing function of distance
  (correlation at the maximum distance); "Nu": a vector
  giving the log of the abundance intensity for each
  strata; "Eta.pois": If Meta$spat.ind==FALSE, spatial
  random effects for Poisson abundance model; one for each
  cell and for each species "Eta.bern": If
  Meta$spat.ind==FALSE & Meta$ZIP=TRUE, spatial random
  effects for Bernoulli abundance model; one for each cell
  and for each species "tau.eta.pois": If
  Meta$spat.ind==FALSE, precision for spatial ICAR model(s)
  for the Poisson component "tau.eta.bern": If
  Meta$spat.ind==FALSE & Meta$ZIP=TRUE, precision for
  spatial ICAR model(s) for the Bernoulli component
  "tau.nu": Precision for Nu (overdispersion relative to
  the Poisson distribution) "G": a vector giving the number
  of groups of animals in each strata; "N": a vector giving
  the number of animals in each strata "MisID": a list,
  each entry i of which is a vector holding the parameters
  for the ith misID model "Cov.par": an (n.species X n X
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
  objects: "iter": number of MCMC iterations; "burnin":
  number of MCMC burnin iterations; "thin": if specified,
  how many iterations to skip between recorded posterior
  samples; "adapt": if adapt==TRUE, this gives the number
  of additional MCMC iterations should be performed to
  adapt MCMC proposals to optimal ranges prior to final
  MCMC run; "MH.cor": Metropolis-hastings tuning parameter
  for updating the correlation parameter (if
  Meta$point.ind==TRUE); "MH.nu": MH tuning parameter for
  Nu parameters (Langevin-Hastings multivariate update);
  "RJ.N": A vector giving the maximum number of additions
  and deletions proposed in an iteration of the RJMCMC
  algorithm for each transect "iter.fix.N" Number of
  iterations to skip RJMCMC step}

  \item{DM.hab.pois}{A design matrix for the Poisson model
  for abundance intensity (log scale)}

  \item{DM.hab.bern}{If Meta$ZIP=TRUE, a design matrix for
  the Bernoulli zero model (probit scale)}

  \item{DM.det}{A design matrix for the probit of detection
  probability}

  \item{Q}{An inverse precision matrix for the spatial ICAR
  process}

  \item{Prior.pars}{A list object giving parameters of
  prior distribution.  Includes the following objects
  "a.eta": alpha parameter for prior precision of spatial
  process (assumed Gamma(a.eta,b.eta)) "b.eta": beta
  parameter for prior precision of spatial process (assumed
  Gamma(a.eta,b.eta)) "a.nu": alpha parameter for prior
  precision of overdispersion process (assumed
  Gamma(a.nu,b.nu)) "b.nu": beta parameter for prior
  precision of overdispersion process (assumed
  Gamma(a.nu,b.nu)) "beta.tau": Prior precision for
  regression coefficients "misID.mu": a list vector, each
  entry gives normal prior means for misID regression
  coefficients for the corresponding model in
  Meta$misID.mat (can be set to null if no misID)
  "misID.sd": a list vector, each entry gives normal prior
  sd for misID regression coefficients for the
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
  "M" Matrix with species-specific rows giving maximum
  possible value for number of groups present in each
  transect (in practice just set high enough that values at
  M and above are never sampled during MCMC) and can be
  fine tuned as needed#' "stacked.names" Character vector
  giving column names for the dataset "factor.ind"
  Indicator vector specifying whether data columns are
  factors (1) or continuous (0) "detect" If TRUE, detection
  parameters are estimated; if FALSE assumes a census
  "Det.formula" a formula object specifying the model for
  the detection process "Levels" a list object, whose
  elements are comprised of detection model names; each
  element gives total # of levels in the combined dataset
  "i.binned" indicator for whether distances are recorded
  in bins (1) or are continuous (0) "dist.pl" gives the
  column in Data where distances are located "G.transect"
  vector holding current number of groups of animals
  present in area covered by each transect "N.transect"
  vector holding current number of animals present in
  covered area by each transect "grps" indicator for
  whether observations are for groups rather than
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
  number of parameters in each covariate pdf "ZIP" If TRUE,
  fit a ZIP model to abundance "point.ind" Indicator for
  whether point independence assumed (if no, then no
  correlation modeled b/w multiple observers as function of
  distance) "last.ind" If TRUE (and point.ind=TRUE), point
  independence operates by assuming 0 dependence at the
  farthest bin "cor.const" If TRUE, forces estimates of
  correlation associated with point independence to be
  positive if last.ind==FALSE or negative if last.ind==TRUE
  (default is FALSE) "fix.tau.nu" Indicator for whether
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
  positive value of misID.mat.  "misID.symm" If TRUE,
  classification probabilities assumed to be symmetric
  (e.g. pi^{2|1}=pi^{1|2}) "N.par.misID" A vector
  specifying the number of parameters needed for each misID
  model "N.hab.pois.par" A vector specifying the number of
  parameters needed for each species' Poisson abundance
  model "N.hab.bern.par" If fitting a ZIP model, this
  vector specifying the number of parameters needed for
  each species' Bernoulli zero model "post.loss" If TRUE,
  observed and predicted detections are compiled for
  posterior predictive loss}
}
\value{
  returns a list with the following objects: "MCMC": An
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
  parameter values (with possible zero inflation) "Post":
  Holds posterior samples for strata specific group sizes
  ("Post$G") and abundance ("Post$N") "Obs.det": if
  Meta$post.loss=TRUE, an array holding observed detection
  types for posterior predictive loss calculations dim =
  c(n.transects,n.obs.types,n.obs.types) "Pred.det": if
  Meta$post.loss=TRUE, an array holding predicted detection
  types for posterior predictive loss calculations dim =
  c(n.mcmc.iter,n.transects,n.obs.types,n.obs.types)
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

