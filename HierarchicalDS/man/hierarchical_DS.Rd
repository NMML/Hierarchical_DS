\name{hierarchical_DS}
\alias{hierarchical_DS}
\title{Primary function for hierarchical, areal analysis of distance sampling data.  This function
pre-processes data and calls other functions to perform the analysis, and is the only function
the user needs to call themselves.}
\usage{
  hierarchical_DS(Dat, Adj, Area.hab = 1, Mapping,
    Area.trans, Observers, Bin.length, Hab.cov, Obs.cov,
    Hab.pois.formula, Hab.bern.formula = NULL, Det.formula,
    detect = TRUE, Cov.prior.pdf, Cov.prior.parms,
    Cov.prior.fixed, Cov.prior.n, n.obs.cov = 0,
    pol.eff = c(1:2), ZIP = FALSE, point.ind = TRUE,
    spat.ind = FALSE, last.ind = FALSE, cor.const = FALSE,
    fix.tau.nu = FALSE, srr = TRUE, srr.tol = 0.5,
    misID = FALSE, misID.models = NULL, misID.mat = NULL,
    misID.symm = TRUE, Inits = NULL, grps = FALSE, M,
    Control, adapt = TRUE, Prior.pars, post.loss = TRUE)
}
\arguments{
  \item{Dat}{A data frame with the following columns:
  (1)transect ID; (2)match number currently, a maximum of 2
  observers on each transect; (3)(Observer ID);
  (4)(Observation (0/1)); (5) Observed species (integer -
  the max integer being 'unknown' if applicable) [NOTE:
  modeled as factor, but need to be input as integers to
  account for unknown species observations] (6-x)(Observer
  covariates); (things like survey conditions or observer
  skill; things that don't change during a transect.  Note
  these also need to be provided in Obs.cov)
  (x+1)(Distance; if all integers, assumed to be discrete
  bins; if continuous, assumed standardized to (0,1)
  interval); (x+2-??)(Group size and other individual
  covariates thought to influence detection; if group size
  is one of them, it's assumed to be column x+2); Note that
  column names can be used to tag covariates, and that
  object types (e.g. numeric, factor) will be preserved in
  analysis}

  \item{Adj}{Adjacency matrix for habitat cells (diagonal
  matrix implies spatial independence)}

  \item{Area.hab}{A vector giving the area of each
  geographical strata (default is equal area)}

  \item{Mapping}{A vector giving the habitat cell id number
  for each transect}

  \item{Area.trans}{A vector giving the effective area
  covered by each transect as fraction of total area in the
  strata it is located}

  \item{Observers}{A (2 x number of transects) matrix
  giving the observers IDs that were present for each
  transect (the 2nd row is to contain NAs if only 1
  observer was present)}

  \item{Bin.length}{If distances are binned, this vector
  gives the relative length of each distance bin (vector
  must sum to one)}

  \item{n.obs.cov}{Number of observer covariates (e.g.,
  seat position, visibility, etc.)}

  \item{Hab.cov}{A data.frame object giving covariates
  thought to influence abundance intensity at strata level;
  column names index individual covariates}

  \item{Obs.cov}{A (max # of observers X # of transects X #
  of observer covariates) size array giving observer
  covariate values for each transect flown}

  \item{Hab.pois.formula}{A list of formulas (one for each
  species) giving the specific model for Poisson abundance
  intensity at the strata level (e.g.,
  ~Vegetation+Latitude) for each species}

  \item{Hab.bern.formula}{If ZIP=TRUE, a list of formulas
  (one for each species) giving the specific model for the
  zero component for abundance intensity at the strata
  level (e.g., ~Vegetation+Latitude) for each species}

  \item{detect}{If TRUE (the default), detectability is
  estimated; if FALSE, assumes detection probability is 1.0
  (i.e. a strip transect with perfect detection).}

  \item{Det.formula}{A formula giving the model for
  detection probability (e.g.
  ~Distance+Group+Visibility+Observer). Note that there are
  several "reserved" variable names.  "Distance",
  "Observer", "Species", and "Group" are reserved variable
  names.}

  \item{Cov.prior.pdf}{If individual covariates are
  provided, this character matrix gives the form of the
  prior pdfs for each covariate current possibilities are
  "poisson",
  "pois1","poisson_ln","pois1_ln",uniform.disc","multinom","uniform.cont",
  or "normal".  "pois1" is 1+x where x~poisson;
  "poisson_ln" and "pois1_ln" are lognormal poisson models
  that incorporate overdispersion.  Note the dimension of
  this matrix are (# species X # of covariates)}

  \item{Cov.prior.parms}{A (s X k X n) array where s is the
  number of species, n is the number of individual
  covariates (other than distance), and k is the maximum
  number of parameters considered for a single covariate
  (NAs can be used to fill this matrix out for covariate
  priors that have <k parameters).  If Cov.prior.fixed=1
  for a given entry, the prior parameters supplied in each
  column apply to the prior pdf itself, and are treated as
  fixed.  If Cov.prior.fixed=0, the model will attempt to
  estimate the posterior distribution of model parameters,
  given hyperpriors.  In this case, it is actually the
  hyperpriors that are being specified.  For "poisson", and
  "pois1", it is assumed that lambda~gamma(alpha,beta), so
  alpha and beta must be supplied.  For "poisson_ln", and
  "pois1_ln", the model is lambda_i=exp(-sigma*Z_i+theta),
  so it is priors for theta and sigma that are specified
  (in that order).  Theta is assumed to have a
  normal(mu,s^2) distribution, and sigma is assumed to have
  a uniform(0,a) distribution; thus, priors are specified
  for these models as (mu,s, and a).  For the multinomial
  pdf, prior parameters of the dirichlet distribution must
  be specified if Cov.prior.fixed=1.}

  \item{Cov.prior.fixed}{An indicator matrix specifying
  which (if any) individual covariate distributions should
  be fixed during estimation}

  \item{Cov.prior.n}{An (# species X # indiv. covariates)
  matrix giving the number of parameters in each covariate
  pdf}

  \item{pol.eff}{For continuous distance, which polynomial
  degrees to model (default is c(1:2); an intercept is
  always estimated when "Distance" is listed in
  "Det.formula")}

  \item{ZIP}{If TRUE, estimate ZIP model for abundance that
  includes a Bernoulli model for zeros and a Poisson + 1
  model for positive values (default is FALSE)}

  \item{point.ind}{Estimate a correlation parameter for
  detection probability that's an increasing function of
  distance?}

  \item{spat.ind}{If TRUE, assumes spatial independence (no
  spatial random effects on abundance intensity) default is
  FALSE}

  \item{last.ind}{If point independence is modeled
  (point.ind=TRUE), last.ind=TRUE will set observer
  covariance to zero for the greatest distance and maximal
  correlation in first bin (default is FALSE)}

  \item{cor.const}{If TRUE, forces estimates of correlation
  associated with point independence to be positive if
  last.ind==FALSE or negative if last.ind==TRUE (default is
  FALSE)}

  \item{fix.tau.nu}{If TRUE, fixes tau.nu during estimation
  (the value to fix it to can be provided in "Inits")}

  \item{srr}{If TRUE, uses spatially retricted regression,
  where smoothing occurs on residuals and all spatial
  effects are orthogonal to the linear predictors (by
  default, analysis is limited to the highest 50
  eigenvalues of the decomposition of the residual
  projection matrix to reduce computing time)}

  \item{srr.tol}{Threshold eigenvalue level for SRR; only
  eigenvectors with higher eigenvalues than srr.tol are
  included in SRR formulation (default is 0.5)}

  \item{misID}{If TRUE, updates species for observed
  animals and estimates misID parameters}

  \item{misID.mat}{With true state on rows and assigned
  state on column, each positive entry provides an index to
  misID.models (i.e. what model to assume on multinomial
  logit space); a 0 indicates an impossible assigment; a
  negative number designates which column is to be obtained
  via subtraction}

  \item{misID.models}{A formula vector providing linear
  model-type formulas for each positive value of
  misID.mat.}

  \item{misID.symm}{If TRUE, the constraint
  pi^{i|j}=pi^{j|i} is implemented; in this case, entries
  for pi^{j|i} are all assumed to be given a '-1' in
  misID.mat}

  \item{grps}{If FALSE, detections are assumed to all be of
  individual animals}

  \item{M}{Matrix with species-specific rows giving maximum
  possible value for number of groups present in each
  transect (in practice just set high enough that values at
  M and above are never sampled during MCMC) and can be
  fine tuned as needed}

  \item{Control}{A list object including the following
  objects: "iter": number of MCMC iterations; "burnin":
  number of MCMC burnin iterations; "thin": if specified,
  how many iterations to skip between recorded posterior
  samples; "adapt": if adapt==TRUE, this gives the number
  of additional MCMC iterations should be performed to
  adapt MCMC proposals to optimal ranges prior to final
  MCMC run; "MH.cor": Metropolis-hastings tuning parameter
  for updating the correlation parameter (if
  point.ind==TRUE); "MH.nu": MH tuning parameters for Nu
  parameters (dimension = # species X # of unique strata
  sampled) "RJ.N": A matrix giving the maximum number of
  additions and deletions proposed in an iteration of the
  RJMCMC algorithm for each species (row) and each transect
  (column) "iter.fix.N": Number of iterations to skip
  RJMCMC step at beginning of estimation (useful for cases
  when estimation is unstable)}

  \item{Inits}{An (optional) list object providing initial
  values for model parameters, with the following objects:
  "hab.pois": Initial values for habitat linear predictor
  parameters for poisson model; "hab.bern": If ZIP=TRUE,
  initial values for habitat linear predictor parameters
  for bernoulli zero model; "det": Initial values for
  detection model (includes distance, observer, env.
  variables, and individual covariates); "cor.par": If
  point.ind==TRUE, this is an initial value for the
  correlation parameter (which must be in (0,1)); "Nu":
  Gives log(lambda) for each spatial strata; "Eta.pois": If
  spat.ind==FALSE, spatial random effects for Poisson
  abundance model; one for each cell and for each species
  "Eta.bern": If spat.ind==FALSE & ZIP=TRUE, spatial random
  effects for Bernoulli abundance model; one for each cell
  and for each species "tau.eta.pois": If spat.ind==FALSE,
  precision for spatial ICAR model(s) for the Poisson
  component "tau.eta.bern": If spat.ind==FALSE & ZIP=TRUE,
  precision for spatial ICAR model(s) for the Bernoulli
  component "tau.nu": Precision for Nu (overdispersion
  relative to the Poisson distribution) One need not
  specify an initial value for all parameter types (if less
  are specified, the others are generated randomly)}

  \item{adapt}{If adapt==TRUE, run an additional
  Control$adapt number of MCMC iterations to optimize MCMC
  proposal distributions prior to primary MCMC}

  \item{Prior.pars}{A list object giving parameters of
  prior distribution.  Includes the following objects
  "a.eta": alpha parameter for prior precision of spatial
  process (assumed Gamma(a.eta,b.eta)) "b.eta": beta
  parameter for prior precision of spatial process (assumed
  Gamma(a.eta,b.eta)) "a.nu": alpha parameter for prior
  precision of overdispersion process (assumed
  Gamma(a.nu,b.nu)) "b.nu": beta parameter for prior
  precision of overdispersion process (assumed
  Gamma(a.nu,b.nu)) "beta.tau": prior precision for
  regression coefficients (assumed
  Normal(0,(beta.tau*X'X)^(-1))}

  \item{post.loss}{If TRUE, calculates observed values and
  posterior predictions for detection data to use with
  posterior predictive loss functions}
}
\value{
  returns a list with the following objecs: MCMC: A list
  object containing posterior samples; Accept: A list
  object indicating the number of proposals that were
  accepted for parameters updated via Metropolis-Hastings;
  Control: A list object giving MCMC tuning parameters
  (which are updated if the 'adapt' alorithm is used)
}
\description{
  Primary function for hierarchical, areal analysis of
  distance sampling data.  This function pre-processes data
  and calls other functions to perform the analysis, and is
  the only function the user needs to call themselves.
}
\examples{
print("example analysis included in the script example_analysis.R")
}
\author{
  Paul B. Conn \email{paul.conn@noaa.gov}
}
\keyword{areal}
\keyword{augmentation,}
\keyword{data}
\keyword{distance}
\keyword{jump}
\keyword{mcmc,}
\keyword{model,}
\keyword{reversible}
\keyword{sampling,}

