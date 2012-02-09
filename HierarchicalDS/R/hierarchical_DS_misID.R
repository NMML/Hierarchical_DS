#' Primary function for hierarchical, areal analysis of distance sampling data.  This function
#' pre-processes data and calls other functions to perform the analysis, and is the only function
#' the user needs to call themselves.
#'
#' @param Dat 	A matrix or data frame with the following columns:
#' 		(1)transect ID; 
#' 		(2)match number  #currently, a maximum of 2 observers on each transect;
#' 		(3)(Observer ID);
#' 		(4)(Observation (0/1));
#' 		(5) Observed species (integer - the max integer being 'unknown' if applicable)
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
#' @param Hab.formula	A formula vector giving the specific model for abundance intensity at the strata level (e.g., ~Vegetation+Latitude) for each species
#' @param Det.formula  A formula giving the model for detection probability (e.g. ~Distance+Group+Visibility+Observer). Note that
#'				there are several "reserved" variable names.  "Distance", "Observer", "Species", and "Group" are reserved variable names.
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
#' @param misID.mat With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#' @param misID.models A formula vector providing linar model-type formulas for each positive value of misID.mat.  If the same model is used in multiple columns it is assumed that all fixed effects (except the intercept) are shared
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
#'	"MH.beta": A matrix of tuning parameters for betas of the abundance process (nrows=number of species, ncol = max number of columns of habitat DM);
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
hierarchical_DS_misID<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,Observers,Bin.length,Hab.cov,Obs.cov,Hab.formula,Det.formula,Cov.prior.pdf,Cov.prior.parms,Cov.prior.fixed,n.obs.cov=0,pol.eff=c(1:2),point.ind=TRUE,spat.ind=FALSE,fix.tau.nu=FALSE,srr=TRUE,srr.tol=0.5,misID.models=NULL,misID.mat=NULL,Inits=NULL,Levels=NA,grps=FALSE,M,Control,adapt=TRUE,Prior.pars){
	require(mvtnorm)
	require(Matrix)
	require(truncnorm)
	require(mc2d)
	
	S=nrow(Adj)
	n.transects=length(Area.trans)
	n.ind.cov=ncol(Dat)-(6+n.obs.cov) #number of individual covariates 
	n.species=nrow(misID.mat)
	n.obs.species=ncol(misID.mat)
	i.unknown=n.obs.species-n.species

	Dat[,3]=as.factor(Dat[,3])  
	n.obs.max=ifelse(sum(is.na(Observers[2,]))==n.transects,1,2)
	cur.colnames=colnames(Dat)
	cur.colnames[6+n.obs.cov]="Distance"
	if(grps==TRUE)cur.colnames[7+n.obs.cov]="Group"
	colnames(Dat)=cur.colnames
	i.binned=0
	if(is.factor(Dat[1,"Distance"])==1){
		i.binned=1
		n.bins=length(unique(Dat[,"Distance"]))
	}

	factor.ind=sapply(Dat[1,],is.factor)
	
	Dat.num=Dat
	for(icol in 1:ncol(Dat)){
		Dat.num[,icol]=as.numeric(as.character(Dat[,icol]))
	}
	Dat.num=as.data.frame(Dat.num)
	
	n.Observers=apply(1-is.na(Observers),2,'sum')
	M=M*n.Observers	#actual dimension of M goes up if >1 observer
	max.M=max(M)
	
	#add an additional column for "True species" and fill
	True.sp=Dat.num[,"Species"]
	unk.ind=which(True.sp==(n.species+1))
	if(length(unk.ind)>0)True.sp[unk.ind]=sample(c(1:n.species),length(unk.ind),replace=TRUE)
	True.sp[which(duplicated(Dat.num[,"Match"])==TRUE)]=True.sp[which(duplicated(Dat.num[,"Match"])==TRUE)-1]
	Dat.num=cbind(Dat.num[1:5],True.sp,Dat.num[6:ncol(Dat.num)])
	
	#Initialize data augmentation multi-d array ("Data"), parameter vectors and matrices
	Data<-array(0,dim=c(n.species,n.transects,max.M,5+n.obs.cov+n.ind.cov)) #array cols are Obs ID,Y,Obs species, True species,Obs covariates,Distance,Ind covariates
	G.transect=matrix(0,n.species,n.transects)  #number of groups by transect; each row gives results for separate species
	n.Records=G.transect #number of records by transect (=G.transect*n.Observers)
	N.transect=G.transect #total abundance by transect

	for(isp in 1:n.species){
		for(itrans in 1:n.transects){
			cur.gt0=sum(Dat.num[,"Transect"]==itrans & Dat.num[,"True.sp"]==isp)
			if(cur.gt0>0){
				Cur.dat=Dat.num[which(Dat.num[,"Transect"]==itrans & Dat.num[,"True.sp"]==isp),3:ncol(Dat.num)]
				Data[isp,itrans,1:nrow(Cur.dat),1:(3+n.obs.cov)]=as.matrix(Cur.dat[,1:(3+n.obs.cov)])
				Data[isp,itrans,1:nrow(Cur.dat),5+n.obs.cov]=as.matrix(Cur.dat[,"Distance"])
				#fill distances for unobserved
				if(i.binned==1)Data[isp,itrans,(nrow(Cur.dat)+1):M[isp,itrans],5+n.obs.cov]=rep(sample(c(1:n.bins),size=(M[isp,itrans]-nrow(Cur.dat))/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
				else Data[isp,itrans,(nrow(Cur.dat)+1):M[isp,itrans],5+n.obs.cov]=rep(runif((M[isp,itrans]-nrow(Cur.dat))/n.Observers[itrans]),each=n.Observers[itrans])
				#fill individual covariate values for (potential) animals that weren't observed
				if(n.ind.cov>0){
					Data[isp,itrans,1:nrow(Cur.dat),(6+n.obs.cov):(6+n.obs.cov+n.ind.cov-1)]=as.matrix(Cur.dat[,(6+n.obs.cov):(6+n.obs.cov+n.ind.cov-1)])
					for(icov in 1:n.ind.cov){
						rsamp=switch_sample(n=(M[isp,itrans]-nrow(Cur.dat))/n.Observers[itrans],pdf=Cov.prior.pdf[isp,icov],cur.par=Cov.prior.parms[isp,,icov],RE=0) 
						Data[isp,itrans,(nrow(Cur.dat)+1):M[isp,itrans],5+n.obs.cov+icov]=rep(rsamp,each=n.Observers[itrans])
					}
				}
				n.Records[isp,itrans]=nrow(Cur.dat)
				G.transect[isp,itrans]=nrow(Cur.dat)/n.Observers[itrans]		#initialize abundance in each transect area to be = to total number of animals observed
				#fill species for latent animals				
				Data[isp,itrans,(n.Records[isp,itrans]+1):max.M,4]=rep(sample(c(1:n.species),(max.M-n.Records[isp,itrans])/n.Observers[itrans],replace=TRUE),each=n.Observers[itrans])
				if(grps==FALSE)N.transect[isp,itrans]=G.transect[isp,itrans]
				else N.transect[isp,itrans]=sum(Cur.dat[,"Group"])/n.Observers[itrans]
			}
			else{
				if(i.binned==1)Data[isp,itrans,1:M[isp,itrans],5+n.obs.cov]=rep(sample(c(1:n.bins),size=M[isp,itrans]/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
				else Data[isp,itrans,1:M[isp,itrans],5+n.obs.cov]=rep(runif(M[isp,itrans]/n.Observers[itrans]),each=n.Observers[itrans])
				if(n.ind.cov>0){
					for(icov in 1:n.ind.cov){
						rsamp=switch_sample(n=M[isp,itrans]/n.Observers[itrans],pdf=Cov.prior.pdf[isp,icov],cur.par=Cov.prior.parms[isp,,icov],RE=0)
						Data[isp,itrans,1:M[isp,itrans],5+n.obs.cov+icov]=rep(rsamp,each=n.Observers[itrans])
					}
				}
				#fill species
				Data[isp,itrans,1:max.M,4]=rep(sample(c(1:n.species),max.M/n.Observers[itrans],replace=TRUE),each=n.Observers[itrans])
				G.transect[isp,itrans]=0
				n.Records[isp,itrans]=0
				N.transect[isp,itrans]=0
			}
			#fill observer ids
			Data[isp,itrans,(n.Records[isp,itrans]+1):max.M,1]=rep(Observers[1:n.Observers[itrans],itrans],(max.M-n.Records[isp,itrans])/n.Observers[itrans])
			#fill observer covariates
			if(n.obs.cov>0){
				for(icov in 1:n.obs.cov){
					Data[isp,itrans,1:max.M,4+icov]=rep(Obs.cov[1:n.Observers[itrans],itrans,icov],max.M/n.Observers[itrans])
				}
			}
		}
	}
	#for debugging, set data aug to truth (keep 0's when simulating data)
	for(isp in 1:n.species){
		for(itrans in 1:n.transects){
			if(G.transect[isp,itrans]>0){
				Cur.dat=Data[isp,itrans,1:n.Records[isp,itrans],]
				if(n.Observers[itrans]==2){
					Obs.ind=matrix(Cur.dat[,2],2,G.transect[isp,itrans])
					Obs.ind[1,]=apply(Obs.ind,2,'max')
					Obs.ind[2,]=Obs.ind[1,]
					Obs.ind=as.vector(Obs.ind)
				}
				else Obs.ind=Cur.dat[,2]
				Data[isp,itrans,1:n.Records[isp,itrans],]=rbind(Cur.dat[which(Obs.ind==1),],Cur.dat[which(Obs.ind==0),])			
			}
		}
	}
	
	stacked.names=c(colnames(Dat)[3:4],"Obs.species","Species",colnames(Dat)[6:ncol(Dat)])
	Stacked=stack_data(Data[1,,,],G.transect[1,]*n.Observers,n.transects,stacked.names,factor.ind) #a stacked form of detection data for updating beta parameters
	for(isp in 2:n.species)Stacked=rbind(Stacked,stack_data(Data[isp,,,],G.transect[isp,]*n.Observers,n.transects,stacked.names,factor.ind))
	Stacked[,"Species"]=as.factor(Stacked[,"Species"])
	
	#determine levels for each factor variable to help in assembling compatible DMs for smaller datasets 
	# (stored in list object named 'Levels')
    factor.ind["Species"]=TRUE
	factor.cols=which(factor.ind[stacked.names]==TRUE) 
	if(length(factor.cols)>0 & is.na(Levels)[1]==1){
		Temp=Stacked[,factor.cols]
		Levels=eval(parse(text=paste('list(',colnames(Temp)[1],'=levels(Temp[,1]))',sep='')))
		if(length(factor.cols)>1){
			for(icol in 2:length(factor.cols)){
				eval(parse(text=paste('Levels$',colnames(Temp)[icol],'=levels(Temp[,icol])',sep='')))	
			}
		}		
	}

	N.hab.par=rep(0,n.species)
	DM.hab=list(sp1=model.matrix(Hab.formula[[1]],data=Hab.cov))
	N.hab.par[1]=ncol(DM.hab$sp1)
	for(i in 2:n.species){  #create design matrices for each species. e.g., name for first species will be DM.hab1
	  eval(parse(text=paste("DM.hab$sp",i,"=model.matrix(Hab.formula[[",i,"]],data=Hab.cov)",sep='')))
	  N.hab.par[i]=eval(parse(text=paste("ncol(DM.hab$sp",i,")",sep='')))
  	}
  	DM.det=get_mod_matrix(Cur.dat=Stacked,stacked.names,factor.ind,Det.formula,Levels)
	
	#now, deal with misID parameters
	N.par.misID=rep(0,max(misID.mat))
	for(i in 1:max(misID.mat)){
		N.par.misID[i]=ncol(model.matrix(misID.models[[i]],data=Stacked))
	}

	Par=generate_inits_misID(DM.hab=DM.hab,DM.det=DM.det,N.hab.par=N.hab.par,G.transect=G.transect,Area.trans=Area.trans,Area.hab=Area.hab,Mapping=Mapping,point.ind=point.ind,spat.ind=spat.ind,grp.mean=Cov.prior.parms[,1,1],misID.mat=misID.mat,N.par.misID=N.par.misID)	
	if(is.null(Inits)==FALSE){  #replace random inits with user provided inits for all parameters specified
		I.init=names(Inits)
		for(ipar in 1:length(I.init)){
			eval(parse(text=paste("Par$",names(Inits)[ipar],"=Inits$",names(Inits[ipar]))))
		}
	}
	#start out at true value for now
	Par$Nu[1,]=DM.hab$sp1%*%Par$hab[1,]
	Par$Nu[2,]=DM.hab$sp2%*%Par$hab[2,]
	#Par$det=c(1.2,-.2,-.4,-.8,-1.4,-1.8,-2,.1,.2,-.4,-.2) 
	#get initial individual covariate parameter values
	Par$Cov.par=Cov.prior.parms 
	for(i in 1:n.ind.cov){	
		for(j in 1:n.species){
			if(Cov.prior.fixed[j,i]==1)Par$Cov.par[j,,i]=Cov.prior.parms[j,,i]
			else{
				temp=switch_sample_prior(Cov.prior.pdf[j,i],Cov.prior.parms[j,,i])
				Par$Cov.par[j,1:length(temp),i]=temp
			}
		}
	}
	
	dist.pl=5+n.obs.cov
	
	n.hab.cov=ifelse(is.null(Hab.cov)==1,0,ncol(Hab.cov))
	
	#Check to make sure input values are internally consistent
	if(n.obs.max>2)cat("\n ERROR: Current max number of observers per transect is 2\n")
	#Later...
	
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

	Meta=list(n.transects=n.transects,n.species=n.species,S=S,spat.ind=spat.ind,Area.hab=Area.hab,Area.trans=Area.trans,
			Adj=Adj,Mapping=Mapping,Covered.area=Covered.area,n.Observers=n.Observers,M=M,stacked.names=stacked.names,
			factor.ind=factor.ind,Det.formula=Det.formula,Levels=Levels,i.binned=i.binned,dist.pl=dist.pl,
			G.transect=G.transect,N.transect=N.transect,grps=grps,n.bins=n.bins,Bin.length=Bin.length,n.ind.cov=n.ind.cov,
			Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,point.ind=point.ind,fix.tau.nu=fix.tau.nu,
			srr=srr,srr.tol=srr.tol,misID.models=misID.models,misID.mat=misID.mat,N.par.misID=N.par.misID,N.hab.par=N.hab.par)
		
	if(adapt==TRUE){
		cat('\n Beginning adapt phase \n')
		Out=mcmc_ds_misID(Par=Par,Data=Data,cur.iter=Control$adapt,adapt=1,Control=Control,DM.hab=DM.hab,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds_misID(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Out$Control,DM.hab=DM.hab,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	else{
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds_misID(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Control,DM.hab=DM.hab,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	Out	
}
