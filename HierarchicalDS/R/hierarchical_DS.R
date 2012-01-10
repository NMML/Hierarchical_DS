#' Primary function for hierarchical, areal analysis of distance sampling data.  This function
#' pre-processes data and calls other functions to perform the analysis, and is the only function
#' the user needs to call themselves.
#'
#' @param Dat 	A matrix or data frame with the following columns:
#' 		(1)transect ID; 
#' 		(2)match number  #currently, a maximum of 2 observers on each transect;
#' 		(3)(Observer ID);
#' 		(4)(Observation (0/1));
#' 		(5-x)(Observer covariates); 
#' 		(x+1)(Distance; if all integers, assumed to be discrete bins; if continuous, assumed standardized to (0,1) interval);
#' 		(x+2-??)(Group size and other individual covariates thought to influence detection; if group size is one of them, it's assumed to be column x+2);
#' 		Note that column names can be used to tag covariates
#' @param Adj   Adjacency matrix for habitat cells (diagonal matrix implies spatial independence)
#' @param Area.hab   A vector giving the area of each geographical strata (default is equal area)
#' @param Mapping  A vector giving the habitat cell id # for each transect
#' @param Area.trans	A vector giving the effective area covered by each transect as fraction of total area in the strata it is located
#' @param Observers	A (2 x number of transects) matrix giving the observers IDs that were present for each transect (the 2nd row is to contain NAs if only 1 observer was present)
#' @param Bin.length	If distances are binned, this vector gives the relative length of each distance bin (vector must sum to one)
#' @param n.obs.cov	Number of observer covariates (e.g., seat position, visibility, etc.)
#' @param Hab.cov	A data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates
#' @param Obs.cov  A (max # of observers X # of transects X # of observer covariates) size array giving observer covariate values for each transect flown
#' @param Hab.formula	A formula object giving the specific model for abundance intensity at the strata level (e.g., ~Vegetation+Latitude)
#' @param Det.formula  A formula giving the model for detection probability (e.g. ~Distance+Group+Visibility+Observer). Note that
#'				there are several "reserved" variable names.  "Distance", "Observer", and "Group" are reserved variable names.
#' @param Cov.prior.pdf	If individual covariates are provided, this character vector gives the form of the prior pdfs for each covariate
#'		  current possibilities are "poisson", "pois1","poisson_ln","pois1_ln",uniform.disc","multinom","uniform.cont", or "normal".
#'		  "pois1" is 1+x where x~poisson; "poisson_ln" and "pois1_ln" are lognormal poisson models that incorporate overdispersion.
#' @param Cov.prior.parms	A (k X n) matrix where n is the number of individual covariates (other than distance), and
#' 		k is the maximum number of parameters considered for a single covariate (NAs can be used to fill this matrix
#'      out for covariate priors that have <k parameters).  If Cov.prior.fixed=1 for a given entry, the prior parameters supplied
#'      in each column apply to the prior pdf itself, and are treated as fixed.  If Cov.prior.fixed=0, the model will attempt
#'  	to estimate the posterior distribution of model parameters, given hyperpriors.  In this case, it is actually the hyperpriors
#'      that are being specified.  For "poisson", and "pois1", it is assumed that lambda~gamma(alpha,beta), so alpha
#' 		and beta must be supplied.  For "poisson_ln", and "pois1_ln", the model is lambda_i=exp(-sigma*Z_i+theta), so it is priors
#' 		for theta and sigma that are specified (in that order).  Theta is assumed to have a normal(mu,s^2) distribution,
#' 		and sigma is assumed to have a uniform(0,a) distribution; thus, priors are specified for these models as (mu,s, and a).
#' 		For the multinomial pdf, prior parameters of the dirichlet distribution must be specified if Cov.prior.fixed=1.
#' @param Cov.prior.fixed  An indicator vector specifying which (if any) individual covariate distributions should be fixed during estimation
#' @param pol.eff 	For continuous distance, which polynomial degrees to model (default is c(1:2); an intercept is always estimated when "Distance" is listed in "Det.formula")
#' @param point.ind  Estimate a correlation parameter for detection probability that's an increasing function of distance?
#' @param spat.ind	If TRUE, assumes spatial independence (no spatial random effects on abundance intensity) default is FALSE
#' @param fix.tau.nu  If TRUE, fixes tau.nu during estimation (the value to fix it to can be provided in "Inits")
#' @param srr  If TRUE, uses spatially retricted regression, where smoothing occurs on residuals and all spatial effects are orthogonal to the linear predictors (by default, analysis is limited to the highest 50 eigenvalues of the decomposition of the residual projection matrix to reduce computing time)
#' @param srr.tol Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation (default is 0.5)
#' @param grps 	If FALSE, detections are assumed to all be of individual animals
#' @param M		Vector giving maximum possible value for number of groups present in each transect (in practice just set high enough that values at M and above are never sampled during MCMC)
#' 			and can be fine tuned as needed
#' @param Levels An optional list object with slots corresponding to factor variable names - giving the possible levels for factors (if not included, the function attempts to ascertain from data)
#' @param Control	A list object including the following slots:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run; 
#'	"MH.cor": Metropolis-hastings tuning parameter for updating the correlation parameter (if point.ind==TRUE);
#'	"MH.nu": MH tuning parameter for Nu parameters (Langevin-Hastings multivariate update);
#'	"MH.beta": A vector of tuning parameters for betas of the abundance process (dimension = number of columns of habitat DM);
#'	"RJ.N"}{A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect
#' @param Inits	An (optional) list object providing initial values for model parameters, with the following slots:
#' "Beta.hab": Initial values for habitat linear predictor parameters;
#'	"Beta.det": Initial values for detection model (includes distance, observer, env. variables, and individual covariates);
#'	"cor.par": If point.ind==TRUE, this is an initial value for the correlation parameter (which must be in (0,1));	
#'	"Nu": Gives log(lambda) for each spatial strata;
#'	"Eta": If spat.ind==FALSE, spatial random effects; one for each strata; 
#'	"tau.eta": If spat.ind==FALSE, precision for spatial ICAR model;  
#'	"tau.nu": Precision for Nu (overdispersion relative to the Poisson distribution)
#'  One need not specify an initial value for all parameter types (if less are specified, the others are generated randomly)
#' @param adapt	If adapt==TRUE, run an additional Control$adapt number of MCMC iterations to optimize MCMC proposal distributions prior to primary MCMC
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following slots
#'	"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'  "b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'	"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
#'	"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
#'	"beta.sd": standard deviation for regression coefficients (assumed Normal(0,beta.sd^2)
#' @return returns a list with the following slots: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' @export
#' @import Matrix
#' @keywords areal model, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn
hierarchical_DS<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,Observers,Bin.length,Hab.cov,Obs.cov,Hab.formula,Det.formula,Cov.prior.pdf,Cov.prior.parms,Cov.prior.fixed,n.obs.cov=0,pol.eff=c(1:2),point.ind=TRUE,spat.ind=FALSE,fix.tau.nu=FALSE,srr=TRUE,srr.tol=0.5,Inits=NULL,Levels=NA,grps=FALSE,M,Control,adapt=TRUE,Prior.pars){
	require(mvtnorm)
	require(Matrix)
	require(truncnorm)
	require(mc2d)
	
	S=nrow(Adj)
	n.transects=length(Area.trans)
	G.transect=Dat
	DM.hab=model.matrix(Hab.formula,data=Hab.cov)

	Par=generate_inits(DM.hab=DM.hab,G.transect=G.transect,Area.trans=Area.trans,Area.hab=Area.hab,Mapping=Mapping,spat.ind=spat.ind)	
	if(is.null(Inits)==FALSE){  #replace random inits with user provided inits for all parameters specified
		I.init=names(Inits)
		for(ipar in 1:length(I.init)){
			eval(parse(text=paste("Par$",names(Inits)[ipar],"=Inits$",names(Inits[ipar]))))
		}
	}
	#start out at true value for now
	#Par$det=c(1.2,-.2,-.4,-.8,-1.4,-1.8,-2,.1,.2,-.4,-.2) 
	Par$Nu=DM.hab%*%Par$hab+rnorm(S,0,0.1)
	#get initial individual covariate parameter values
	n.hab.cov=ifelse(is.null(Hab.cov)==1,0,ncol(Hab.cov))
		
	i.Covered=c(1:S)%in%Mapping
	Covered.area=rep(0,S)
	for(i in 1:S){
		if(i.Covered[i]==1){
			Covered.area[i]=sum(Area.trans[which(Mapping==i)])
		}
	}
	
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q=Matrix(Q)	

	Meta=list(n.transects=n.transects,S=S,spat.ind=spat.ind,Area.hab=Area.hab,Area.trans=Area.trans,
			Adj=Adj,Mapping=Mapping,Covered.area=Covered.area,
			G.transect=G.transect,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol)
		
	if(adapt==TRUE){
		cat('\n Beginning adapt phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$adapt,adapt=1,Control=Control,DM.hab=DM.hab,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Out$Control,DM.hab=DM.hab,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	else{
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Control,DM.hab=DM.hab,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	Out	
}
