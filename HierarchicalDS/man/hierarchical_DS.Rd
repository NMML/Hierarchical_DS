\name{hierarchical_DS}
\alias{hierarchical_DS}
\title{Primary function for hierarchical, areal analysis of distance sampling data.  This function
pre-processes data and calls other functions to perform the analysis, and is the only function
the user needs to call themselves.}
\usage{
  hierarchical_DS(Dat, Adj, Area.hab = 1, Mapping,
    Area.trans, Observers, Bin.length, Hab.cov, Obs.cov,
    Hab.formula, Det.formula, Cov.prior.pdf,
    Cov.prior.parms, Cov.prior.fixed, Cov.prior.n,
    n.obs.cov = 0, pol.eff = c(1:2), point.ind = TRUE,
    spat.ind = FALSE, fix.tau.nu = FALSE, srr = TRUE,
    srr.tol = 0.5, misID = FALSE, misID.models = NULL,
    misID.mat = NULL, Inits = NULL, Levels = NA,
    grps = FALSE, M, Control, adapt = TRUE, Prior.pars)
}
\arguments{
  \item{Dat}{A matrix or data frame with the following
  columns: (1)transect ID; (2)match number #currently, a
  maximum of 2 observers on each transect; (3)(Observer
  ID); (4)(Observation (0/1)); (5) Observed species
  (integer - the max integer being 'unknown' if applicable)
  (5-x)(Observer covariates); (x+1)(Distance; if all
  integers, assumed to be discrete bins; if continuous,
  assumed standardized to (0,1) interval); (x+2-??)(Group
  size and other individual covariates thought to influence
  detection; if group size is one of them, it's assumed to
  be column x+2); Note that column names can be used to tag
  covariates}

  \item{Adj}{Adjacency matrix for habitat cells (diagonal
  matrix implies spatial independence)}

  \item{Area.hab}{A vector giving the area of each
  geographical strata (default is equal area)}

  \item{Mapping}{A vector giving the habitat cell id # for
  each transect}

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

  \item{Hab.formula}{A formula vector giving the specific
  model for abundance intensity at the strata level (e.g.,
  ~Vegetation+Latitude) for each species}

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

  \item{point.ind}{Estimate a correlation parameter for
  detection probability that's an increasing function of
  distance?}

  \item{spat.ind}{If TRUE, assumes spatial independence (no
  spatial random effects on abundance intensity) default is
  FALSE}

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

  \item{misID.models}{A formula vector providing linar
  model-type formulas for each positive value of misID.mat.
  If the same model is used in multiple columns it is
  assumed that all fixed effects (except the intercept) are
  shared}

  \item{grps}{If FALSE, detections are assumed to all be of
  individual animals}

  \item{M}{Vector giving maximum possible value for number
  of groups present in each transect (in practice just set
  high enough that values at M and above are never sampled
  during MCMC) and can be fine tuned as needed}

  \item{Levels}{An optional list object with slots
  corresponding to factor variable names - giving the
  possible levels for factors (if not included, the
  function attempts to ascertain from data)}

  \item{Control}{A list object including the following
  slots: "iter": number of MCMC iterations; "burnin":
  number of MCMC burnin iterations; "thin": if specified,
  how many iterations to skip between recorded posterior
  samples; "adapt": if adapt==TRUE, this gives the number
  of additional MCMC iterations should be performed to
  adapt MCMC proposals to optimal ranges prior to final
  MCMC run; "MH.cor": Metropolis-hastings tuning parameter
  for updating the correlation parameter (if
  point.ind==TRUE); "MH.nu": MH tuning parameter for Nu
  parameters (Langevin-Hastings multivariate update);
  "MH.beta": A matrix of tuning parameters for betas of the
  abundance process (nrows=number of species, ncol = max
  number of columns of habitat DM); "RJ.N"}{A vector giving
  the maximum number of additions and deletions proposed in
  an iteration of the RJMCMC algorithm for each transect}

  \item{Inits}{An (optional) list object providing initial
  values for model parameters, with the following slots:
  "Beta.hab": Initial values for habitat linear predictor
  parameters; "Beta.det": Initial values for detection
  model (includes distance, observer, env. variables, and
  individual covariates); "cor.par": If point.ind==TRUE,
  this is an initial value for the correlation parameter
  (which must be in (0,1)); "Nu": Gives log(lambda) for
  each spatial strata; "Eta": If spat.ind==FALSE, spatial
  random effects; one for each strata; "tau.eta": If
  spat.ind==FALSE, precision for spatial ICAR model;
  "tau.nu": Precision for Nu (overdispersion relative to
  the Poisson distribution) One need not specify an initial
  value for all parameter types (if less are specified, the
  others are generated randomly)}

  \item{adapt}{If adapt==TRUE, run an additional
  Control$adapt number of MCMC iterations to optimize MCMC
  proposal distributions prior to primary MCMC}

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
  regression coefficients (assumed Normal(0,beta.sd^2)}
}
\value{
  returns a list with the following slots: MCMC: A list
  object containing posterior samples; Accept: A list
  object indicating the number of proposals that were
  accepted for parameters updated via Metropolis- or
  Langevin-Hastings algorithms; Control: A list object
  giving MCMC tuning parameters (which are updated if the
  'adapt' alorithm is used)
}
\description{
  Primary function for hierarchical, areal analysis of
  distance sampling data.  This function pre-processes data
  and calls other functions to perform the analysis, and is
  the only function the user needs to call themselves.
}
\examples{
# This example simulates and analyzes data from 10 transects and one study site; however
# the time to run is still about ~1 day, so it is recommended NOT to actually run it but do use "data(sim.data)" instead
S=10
set.seed(12345)
n.real.transects=S
n.transects=S #each transect spans two cells
#generate Observers matrix
Observers=matrix(NA,2,n.transects)
Obs.cov=array(0,dim=c(2,n.transects,1))
Obs.cov[1,,]=1
n.obs.cov=0
for(i in 1:n.real.transects)Observers[,i]=sample(c(1,2),size=2,replace=FALSE)
# provide levels of for each covariate of type "Factor"
Levels=list(Observer=c("1","2"),Distance=c("1","2","3","4","5"),Species=c("1","2"))
n.species=1
Out=simulate_data(S=S,Observers=Observers,misID=FALSE,n.species=1,Beta.hab=matrix(c(log(20),1),1,2),X.site=cbind(rep(1,S),log(c(1:S)/S)),detect.model=~Observer+Distance+Group,Beta.det=c(1.2,-.4,-.8,-1.1,-1.3,-1.5,.1))  #misID currently not implemented/debugged fully
Dat=Out$Dat
Mapping=c(1:S)
Area.trans=rep(1,S)
n.bins=length(unique(Dat[,"Distance"]))
Area.hab=rep(1,S)
Bin.length=rep(1,n.bins)
Adj=linear_adj(S)

Hab.cov=data.frame(matrix(log(c(1:S)/S),ncol=1)) #covariate on abundance intensity
colnames(Hab.cov)=c("Cov1")
Hab.formula=list(~Cov1)
Det.formula=~Observer+Distance+Group
misID.mat=NULL
misID.models=NULL

#set priors for individual covariates
Cov.prior.parms=array(0,dim=c(n.species,2,1))
Cov.prior.parms[1,,1]=c(1.1,1)
Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
Cov.prior.pdf=Cov.prior.fixed
Cov.prior.pdf[1]="pois1"
Cov.prior.n=matrix(2,1,1)

pol.eff=NULL #not currently used since using distance bins
point.ind=TRUE
spat.ind=TRUE #dont' include spatial dependence unless there really is spatial structure!
fix.tau.nu=FALSE
srr=FALSE
srr.tol=0.2
misID=FALSE
grps=TRUE
#set initial values for M (max number of groups allowed in a given transect)
M=t(Out$G.true*3)
M[which(M<20)]=20
Control=list(iter=30100,burnin=100,thin=10,MH.cor=0.2,MH.nu=c(.01),MH.misID=NULL,RJ.N=matrix(rep(5,S*n.species),n.species,S),adapt=1000)
#provide initial values for habitat parameters to improve stability
hab=matrix(0,n.species,2) #covariates are intercept, index,
hab[1,]=c(log(20),0)
Inits=list(hab=hab,tau.nu=c(100))
#additional prior parameters
Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.sd=100) #(1,.01) prior makes it closer to a uniform distribution near the origin
adapt=TRUE

set.seed(8327329)    #chain1
#set.seed(8327330)   #chain2

# CAUTION: this next line takes several hours to run!  Recommend using data(sim_data) instead!
#Out=hierarchical_DS(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Observers=Observers,Bin.length=Bin.length,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.formula=Hab.formula,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,pol.eff=NULL,point.ind=TRUE,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,misID=misID,Inits=Inits,grps=grps,M=M,Control=Control,Levels=Levels,adapt=TRUE,Prior.pars=Prior.pars,misID.mat=misID.mat,misID.models=misID.models)
data(simdata)
#look at some plots and summary statistics
plot_obs_pred(simdata)
#table(simdata$MCMC,a=0.025)
dens.plot(simdata$MCMC)
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

