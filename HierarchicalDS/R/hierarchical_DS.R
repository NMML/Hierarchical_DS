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
#' @param Bin.length	If distances are binned, this vector gives the relative length of each distance bin (vector must sum to one)
#' @param n.obs.cov	Number of observer covariates (e.g., seat position, visibility, etc.)
#' @param Hab.cov	A data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates
#' @param Hab.formula	A formula object giving the specific model for abundance intensity at the strata level (e.g., ~Vegetation+Latitude)
#' @param Det.formula  A formula giving the model for detection probability (e.g. ~Distance+Group+Visibility+Observer). Note that
#'				there are several "reserved" variable names.  "Distance", "Observer", and "Group" are reserved variable names.
#' @param Cov.prior.pdf	If individual covariates are provided, this character vector gives the form of the prior pdfs for each covariate
#'		  current possibilities are "poisson", "pois1","uniform.disc", "uniform.cont", or "normal".
#'		  "pois1" is 1+x where x~poisson
#' @param Cov.prior.parms	A (2 X n) matrix where n is the number of individual covariates (other than distance).  Each column
#'		gives the parameters associated with the prior pdf of a given covariate (only the value in the first row is
#'		used for "poisson"; for normal, first row gives mean, second row gives sd; for uniform, first row gives lower,
#'		second row gives upper; for constant, the parameter entries are just placeholders (no parameters are required)
#'		note that these priors are also used to propose covariate values during RJMCMC, so should be 
#'		made to be biologically plausible (i.e., don't use 'vague' priors!)
#' @param pol.eff 	For continuous distance, which polynomial degrees to model (default is c(1:2); an intercept is always estimated when "Distance" is listed in "Det.formula")
#' @param point.ind  Estimate a correlation parameter for detection probability that's an increasing function of distance?
#' @param spat.ind	If TRUE, assumes spatial independence (no spatial random effects on abundance intensity) default is FALSE
#' @param grps 	If FALSE, detections are assumed to all be of individual animals
#' @param M		Vector giving maximum possible value for number of groups present in each transect (in practice just set high enough that values at M and above are never sampled during MCMC)
#' 			and can be fine tuned as needed
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
hierarchical_DS<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,Bin.length,Hab.cov,Hab.formula,Det.formula,Cov.prior.pdf,Cov.prior.parms,n.obs.cov=0,pol.eff=c(1:2),point.ind=TRUE,spat.ind=FALSE,Inits=NULL,grps=FALSE,M,Control,adapt=TRUE,Prior.pars){
	require(mvtnorm)
	require(Matrix)
	require(truncnorm)

	S=nrow(Adj)
	n.transects=length(unique(Dat[,1]))
	n.Observers=rep(0,n.transects)
	n.ind.cov=ncol(Dat)-(5+n.obs.cov) #number of individual covariates

	unique.observers=unique(Dat[,3])
	Dat[,3]=as.factor(Dat[,3])  
	n.observers=length(unique.observers)
	n.obs.max=0
	for(i in 1:max(Dat[,2])){
		n.obs.max=max(n.obs.max,length(which(Dat[,2]==i)))  #maximum number of observers on a transect
	}
	n.transects=length(Area.trans)
	cur.colnames=colnames(Dat)
	cur.colnames[5+n.obs.cov]="Distance"
	if(grps==TRUE)cur.colnames[6+n.obs.cov]="Group"
	colnames(Dat)=cur.colnames
	i.binned=0
	if(sum(Dat[,"Distance"]%%1)==0){
		i.binned=1
		Dat[,"Distance"]=as.factor(Dat[,"Distance"])
		n.bins=length(unique(Dat[,"Distance"]))
	}

	factor.ind=sapply(Dat[1,],is.factor)
	
	Dat.num=Dat
	for(icol in 1:ncol(Dat)){
		Dat.num[,icol]=as.numeric(Dat[,icol])
	}
	Dat.num=as.matrix(Dat.num)
	
    #compute number of observers in each transect
	for(itrans in 1:n.transects){
		Cur.dat=Dat.num[which(Dat.num[,1]==itrans),2:ncol(Dat.num)]
		first.match=Cur.dat[1,1]
		n.Observers[itrans]=length(which(Cur.dat[,1]==first.match))
	}
	M=M*n.Observers	#actual dimension of M goes up if >1 observer
	
	#Initialize data augmentation multi-d array ("Data"), parameter vectors and matrices
	Data<-array(0,dim=c(n.transects,max(M),3+n.obs.cov+n.ind.cov)) #array cols are Obs ID,Y,Obs covariates,Distance,Ind covariates
	G.transect=rep(0,n.transects)  #number of groups by transect
	n.Records=rep(0,n.transects) #number of records by transect (=G.transect*n.Observers)
	N.transect=G.transect #total abundance by transect

	for(itrans in 1:n.transects){
		Cur.dat=Dat.num[which(Dat.num[,1]==itrans),2:ncol(Dat.num)]
		Cur.dat=Cur.dat[,-1]
		Data[itrans,1:nrow(Cur.dat),1:(2+n.obs.cov)]=Cur.dat[,1:(2+n.obs.cov)]
		Data[itrans,1:nrow(Cur.dat),3+n.obs.cov]=Cur.dat[,"Distance"]
		#fill distances for unobserved
		if(i.binned==1)Data[itrans,(nrow(Cur.dat)+1):M[itrans],3+n.obs.cov]=rep(sample(c(1:n.bins),size=(M[itrans]-nrow(Cur.dat))/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
		else Data[itrans,(nrow(Cur.dat)+1):M[itrans],3+n.obs.cov]=rep(runif((M[itrans]-nrow(Cur.dat))/n.Observers[itrans]),each=n.Observers[itrans])
		#fill individual covariate values for (potential) animals that weren't observed
		if(n.ind.cov>0){
			Data[itrans,1:nrow(Cur.dat),(4+n.obs.cov):(4+n.obs.cov+n.ind.cov-1)]=Cur.dat[,(4+n.obs.cov):(4+n.obs.cov+n.ind.cov-1)]
			for(icov in 1:n.ind.cov){
				rsamp=switch_sample(n=(M[itrans]-nrow(Cur.dat))/n.Observers[itrans],pdf=Cov.prior.pdf[icov],cur.par=Cov.prior.parms[,icov])
				Data[itrans,(nrow(Cur.dat)+1):M[itrans],3+n.obs.cov+icov]=rep(rsamp,each=n.Observers[itrans])
			}
		}
		G.transect[itrans]=nrow(Cur.dat)/n.Observers[itrans]		#initialize abundance in each transect area to be = to total number of animals observed
		n.Records[itrans]=G.transect[itrans]*n.Observers[itrans]
		N.transect[itrans]=ifelse(grps==FALSE,G.transect[itrans],sum(Cur.dat[,which(colnames(Cur.dat)=="Group")])/n.Observers[itrans])
		#fill observer ids
		Data[itrans,(n.Records[itrans]+1):M[itrans],1]=rep(Data[itrans,1:n.Observers[itrans],1],(M[itrans]-n.Records[itrans])/n.Observers[itrans])
		curcol=3
		#fill observer covariates
		if(n.obs.cov>0){
			for(icol in curcol:(curcol+n.obs.cov-1)){
				Data[itrans,(n.Records[itrans]+1):M[itrans],curcol]=rep(Data[itrans,1:n.Observers[itrans],curcol],(M[itrans]-n.Records[itrans])/n.Observers[itrans])
			}
			curcol=curcol+n.obs.cov
		}
	}
	stacked.names=colnames(Dat)[3:ncol(Dat)]
	Stacked=stack_data(Data,G.transect,n.transects,stacked.names,factor.ind) #a stacked form of detection data for updating beta parameters
	
	#determine levels for each factor variable to help in assembling compatible DMs for smaller datasets 
	# (stored in list object named 'Levels')
	factor.cols=which(factor.ind[stacked.names]==TRUE) 
	if(length(factor.cols>0)){
		Temp=Stacked[,factor.cols]
		Levels=eval(parse(text=paste('list(',colnames(Temp)[1],'=levels(Temp[,1]))',sep='')))
		if(length(factor.cols)>1){
			for(icol in 2:length(factor.cols)){
				eval(parse(text=paste('Levels$',colnames(Temp)[icol],'=levels(Temp[,icol])',sep='')))	
			}
		}		
	}
	
	DM.hab=model.matrix(Hab.formula,data=Hab.cov)
	DM.det=model.matrix(Det.formula,data=Stacked)	
	
	Par=Inits
	if(is.null(Inits)==TRUE)Par=generate_inits(DM.hab=DM.hab,DM.det=DM.det,G.transect=G.transect,Area.trans=Area.trans,Area.hab=Area.hab,Mapping=Mapping,point.ind=point.ind,spat.ind=spat.ind,grp.mean=Cov.prior.parms[1,1])	
	#start out at true value for now
	Par$hab=c(log(100),1)
	Par$Nu=DM.hab%*%Par$hab
	
	dist.pl=3+n.obs.cov
	
	n.hab.cov=ifelse(is.null(Hab.cov)==1,0,ncol(Hab.cov))
	N.remain=Adj*0  #keep track of "remaining" abundance in each strata for portion of area not covered by transects
	N.total=N.remain #keep track of total abundance by strata
	
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

	Meta=list(n.transects=n.transects,S=S,spat.ind=spat.ind,Area.hab=Area.hab,Area.trans=Area.trans,
			Mapping=Mapping,Covered.area=Covered.area,n.Observers=n.Observers,M=M,stacked.names=stacked.names,
			factor.ind=factor.ind,Det.formula=Det.formula,Levels=Levels,i.binned=i.binned,dist.pl=dist.pl,
			G.transect=G.transect,N.transect=N.transect,grps=grps,n.bins=n.bins,Bin.length=Bin.length,n.ind.cov=n.ind.cov,
			Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,point.ind=point.ind)
		
	if(adapt==TRUE){
		cat('\n Beginning adapt phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$adapt,adapt=1,Control=Control,DM.hab=DM.hab,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Out$Control,DM.hab=DM.hab,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	else{
		cat('\n Beginning MCMC phase \n')
		Out=mcmc_ds(Par=Par,Data=Data,cur.iter=Control$iter,adapt=0,Control=Control,DM.hab=DM.hab,DM.det=DM.det,Q=Q,Prior.pars=Prior.pars,Meta=Meta)
	}
	
	Out	
}

	