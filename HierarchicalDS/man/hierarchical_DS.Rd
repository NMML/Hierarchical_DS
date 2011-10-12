\name{hierarchical_DS}
\alias{hierarchical_DS}
\title{Primary function for hierarchical, areal analysis of distance sampling data.}
\usage{hierarchical_DS(Dat, Adj, Area.hab=1, Mapping, Area.trans, Bin.length,
    Hab.cov, Hab.formula, Det.formula, Cov.prior.pdf, Cov.prior.parms,
    n.obs.cov=0, pol.eff=c(1:2), point.ind=TRUE, spat.ind=FALSE, Inits,
    grps=FALSE, M, Control, adapt=TRUE, Prior.pars)}
\description{Primary function for hierarchical, areal analysis of distance sampling data.  This function
pre-processes data and calls other functions to perform the analysis, and is the only function
the user needs to call themselves.}
\value{returns a list with the following slots: 
MCMC: A list object containing posterior samples;
Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)}
\keyword{areal}
\keyword{model,}
\keyword{data}
\keyword{augmentation,}
\keyword{distance}
\keyword{sampling,}
\keyword{mcmc,}
\keyword{reversible}
\keyword{jump}
\author{Paul B. Conn}
\arguments{\item{Dat}{A matrix or data frame with the following columns:
(1)transect ID; 
(2)match number  #currently, a maximum of 2 observers on each transect;
(3)(Observer ID);
(4)(Observation (0/1));
(5-x)(Observer covariates); 
(x+1)(Distance; if all integers, assumed to be discrete bins; if continuous, assumed standardized to (0,1) interval);
(x+2-??)(Group size and other individual covariates thought to influence detection; if group size is one of them, it's assumed to be column x+2);
Note that column names can be used to tag covariates}
\item{Adj}{Adjacency matrix for habitat cells (diagonal matrix implies spatial independence)}
\item{Area.hab}{A vector giving the area of each geographical strata (default is equal area)}
\item{Mapping}{A vector giving the habitat cell id # for each transect}
\item{Area.trans}{A vector giving the effective area covered by each transect as fraction of total area in the strata it is located}
\item{Bin.length}{If distances are binned, this vector gives the relative length of each distance bin (vector must sum to one)}
\item{n.obs.cov}{Number of observer covariates (e.g., seat position, visibility, etc.)}
\item{Hab.cov}{A data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates}
\item{Hab.formula}{A formula object giving the specific model for abundance intensity at the strata level (e.g., ~Vegetation+Latitude)}
\item{Det.formula}{A formula giving the model for detection probability (e.g. ~Distance+Group+Visibility+Observer). Note that
there are several "reserved" variable names.  "Distance", "Observer", and "Group" are reserved variable names.}
\item{Cov.prior.pdf}{If individual covariates are provided, this character vector gives the form of the prior pdfs for each covariate
current possibilities are "poisson", "pois1","uniform.disc", "uniform.cont", or "normal".
"pois1" is 1+x where x~poisson}
\item{Cov.prior.parms}{A (2 X n) matrix where n is the number of individual covariates (other than distance).  Each column
gives the parameters associated with the prior pdf of a given covariate (only the value in the first row is
used for "poisson"; for normal, first row gives mean, second row gives sd; for uniform, first row gives lower,
second row gives upper; for constant, the parameter entries are just placeholders (no parameters are required)
note that these priors are also used to propose covariate values during RJMCMC, so should be 
made to be biologically plausible (i.e., don't use 'vague' priors!)}
\item{pol.eff}{For continuous distance, which polynomial degrees to model (default is c(1:2); an intercept is always estimated when "Distance" is listed in "Det.formula")}
\item{point.ind}{Estimate a correlation parameter for detection probability that's an increasing function of distance?}
\item{spat.ind}{If TRUE, assumes spatial independence (no spatial random effects on abundance intensity) default is FALSE}
\item{grps}{If FALSE, detections are assumed to all be of individual animals}
\item{M}{Vector giving maximum possible value for number of groups present in each transect (in practice just set high enough that values at M and above are never sampled during MCMC)
and can be fine tuned as needed}
\item{Control}{A list object including the following slots:
"iter": number of MCMC iterations;
"burnin": number of MCMC burnin iterations;
"thin": if specified, how many iterations to skip between recorded posterior samples;
"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run; 
"MH.cor": Metropolis-hastings tuning parameter for updating the correlation parameter (if point.ind==TRUE);
"MH.nu": MH tuning parameter for Nu parameters (Langevin-Hastings multivariate update);
"MH.beta": A vector of tuning parameters for betas of the abundance process (dimension = number of columns of habitat DM);
"RJ.N"}{A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect}
\item{Inits}{An (optional) list object providing initial values for model parameters, with the following slots:
"Beta.hab": Initial values for habitat linear predictor parameters;
"Beta.det": Initial values for detection model (includes distance, observer, env. variables, and individual covariates);
"cor.par": If point.ind==TRUE, this is an initial value for the correlation parameter (which must be in (0,1));	
"Nu": Gives log(lambda) for each spatial strata;
"Eta": If spat.ind==FALSE, spatial random effects; one for each strata; 
"tau.eta": If spat.ind==FALSE, precision for spatial ICAR model;  
"tau.nu": Precision for Nu (overdispersion relative to the Poisson distribution)}
\item{adapt}{If adapt==TRUE, run an additional Control$adapt number of MCMC iterations to optimize MCMC proposal distributions prior to primary MCMC}
\item{Prior.pars}{A list object giving parameters of prior distribution.  Includes the following slots
"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
"b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
"beta.sd": standard deviation for regression coefficients (assumed Normal(0,beta.sd^2)}}
