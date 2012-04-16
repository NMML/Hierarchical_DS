#' Function for MCMC analysis 
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		"det": a vector giving the current iteration's linear model parameters for the detection model;
#' 		"hab": a vector giving the current iteration's linear model parameters for abundance intensity;
#' 		"cor": a correlation parameter for detections that's an increasing function of distance (correlation at the maximum distance);
#' 		"Nu": a vector giving the log of the abundance intensity for each strata;
#' 		"G": a vector giving the number of groups of animals in each strata; 
#' 		"N": a vector giving the number of animals in each strata
#' 		"Cov.par": an (n.species X n X n.ind.cov)  array holding parameters of individual covariate distributions.
#' @param Data   A four dimensional array; the first dimension gives species, the second gives the transect, the third indexes a (possible) observation, 
#' 			and the fourth dimension gives observations and covariates associated with a given animal.
#' 			These final columns are: Observer ID,Y(observation=0/1),Observed species,Obs covariates,Distance,Ind covariates
#' @param cur.iter   Number of iterations to run
#' @param adapt	If adapt==TRUE, run MCMC in adapt mode, optimizing MCMC proposal distributions prior to primary MCMC
#' @param Control	A list object including the following slots:
#'	"iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'	"adapt": if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal 
#' 				ranges prior to final MCMC run; 
#'	"MH.cor": Metropolis-hastings tuning parameter for updating the correlation parameter (if Meta$point.ind==TRUE);
#'	"MH.nu": MH tuning parameter for Nu parameters (Langevin-Hastings multivariate update);
#'	"MH.beta": A matrix of tuning parameters for betas of the abundance process (nrows=number of species, ncol = max number of columns of habitat DM);
#'	"RJ.N": A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect
#' @param DM.hab	A design matrix for the log of abundance intensity
#' @param DM.det	A design matrix for the probit of detection probability
#' @param Q			An inverse precision matrix for the spatial ICAR process
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following slots
#'	"a.eta": alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'  "b.eta": beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))
#'	"a.nu": alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))
#'	"b.nu": beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu)) 
#'	"beta.sd": standard deviation for regression coefficients (assumed Normal(0,beta.sd^2)
#'  "misID.mu": a list vector, each entry gives normal prior means for misID regression coefficients for the corresponding model in Meta$misID.mat (can be set to null if no misID)
#'  "misID.sd": a list vector, each entry gives normal prior sd for misID regression coefficients for the corresponding model in Meta$misID.mat (can be set to null if no misID)
#' @param Meta	A list object giving a number of other features of the dataset, including:
#' 	"n.transects"	Number of transects
#'  "n.species"     Number of species
#' 	"S"				Number of strata cells
#'  "spat.ind"		Indicator for spatial dependence
#'  "Area.hab"		Vector giving relative area covered by each strata
#'  "Area.trans"	Vector giving fraction of area of relevant strata covered by each transect
#' 	"Adj"			Adjacency matrix giving connectivity of spatial grid cells
#'  "Mapping" 		Vector mapping each transect into a parent strata
#'  "Covered.area"	Vector giving the fraction of each strata covered by transects
#' 	"n.Observers"	Vector giving the number of observers that operated on each transect
#'  "M"				Vector giving maximum number of groups in each transect
#'  "stacked.names" Character vector giving column names for the dataset
#'  "factor.ind"	Indicator vector specifying whether data columns are factors (1) or continuous (0)
#'  "Det.formula"	a formula object specifying the model for the detection process
#'  "Levels"		a list object, where slot names are comprised of detection model names; each slot gives total # of levels in the combined dataset
#'  "i.binned"		indicator for whether distances are recorded in bins (1) or are continuous (0)
#'  "dist.pl"		gives the column in Data where distances are located	
#'  "G.transect"	vector holding current number of groups of animals present in area covered by each transect		
#'  "N.transect"    vector holding current number of animals present in covered area by each transect
#'  "grps"			indicator for whether observations are for groups rather than individuals
#'  "n.bins"		number of distance bins (provided i.binned=1)
#'  "Bin.length"	vector giving relative size of distance bins
#'  "n.ind.cov" 	Number of individual covariates (distance is not included in this total, but group size is)
#'  "Cov.prior.pdf" character vector giving the probability density function associated with each individual covariate (type ? hierarchical_DS for more info)
#'  "Cov.prior.parms"	An (n.species X n X n.ind.cov) array providing "pseudo-prior" parameters for individual covarate distributions (only the first row used if a signle parameter distribution)
#'  "Cov.prior.fixed" indicator vector for whether parameters of each covariate distribution should be fixed within estimation routine
#'  "Cov.prior.n" 	(#species X #covariates) Matrix giving number of parameters in each covariate pdf 
#'  "point.ind"		Indicator for whether point independence assumed (if no, then no correlation modeled b/w multiple observers as function of distance)
#'  "fix.tau.nu"	Indicator for whether tau.nu should be fixed (1) or estimated(0)
#'  "srr"			Indicator for whether a spatially restricted regression model should be employed (1) or not (0)
#'  "srr.tol"		Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation
#'  "misID"			If TRUE, misidentification of species is modeled
#'  "misID.mat"     With true state on rows and assigned state on column, each positive entry provides an index to misID.models (i.e. what model to assume on multinomial logit space); a 0 indicates an impossible assigment; a negative number designates which column is to be obtained via subtraction
#'  "misID.models"  A formula vector providing linar model-type formulas for each positive value of misID.mat.  If the same model is used in multiple columns it is assumed that all fixed effects (except the intercept) are shared
#'  "N.par.misID"   A vector specifying the number of parameters needed for each misID model
#'  "N.hab.par"	    A vector specifying the number of parameters needed for each species' habitat model
#' @return returns a list with the following slots: 
#' 	"MCMC": An 'mcmc' object (see 'coda' R package) containing posterior samples;
#'  "Accept": A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  "Control": A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used) 
#'  "Obs.N":  Records latent abundance in each transect; dimension is (n.species X # samples X # transects)
#'  "Pred.N": Posterior predictive distribution for abundance in each transect; obtained by sampling a Poisson distribution given current parameter values
#'  "Post": Holds posterior samples for strata specific group sizes ("Post$G") and abundance ("Post$N")
#' @export
#' @import Matrix
#' @keywords areal, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn

mcmc_ds<-function(Par,Data,cur.iter,adapt,Control,DM.hab,DM.det,Q,Prior.pars,Meta){	
	#require(mvtnorm)
	#require(Matrix)
	#require(truncnorm)
	require(compiler)
	switch_sample<-cmpfun(switch_sample)
	switch_sample_prior<-cmpfun(switch_sample_prior)
	switch_pdf<-cmpfun(switch_pdf)
	get_mod_matrix<-cmpfun(get_mod_matrix)
	log_lambda_gradient<-cmpfun(log_lambda_gradient)
	log_lambda_log_likelihood<-cmpfun(log_lambda_log_likelihood)

	
	Lam.index=c(1:Meta$S)
	if(Meta$i.binned==0)dist.mult=1
	if(Meta$i.binned==1)dist.mult=1/(Meta$n.bins-1)
	n.beta.det=ncol(DM.det)
	n.Records=Meta$G.transect*Meta$n.Observers
	
	#initialize G.obs (number of groups observed per transect)
	G.obs=Meta$G.transect
	for(isp in 1:Meta$n.species){
		for(itrans in 1:Meta$n.transects){
			Tmp=matrix(Data[isp,itrans,1:n.Records[isp,itrans],2],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
			G.obs[isp,itrans]=sum(apply(Tmp,1,'sum')>0)
		}
	}
	n.samp.misID=max(1,round(0.05*sum(G.obs)))  #currently only updating species for 1/10 of population at each iteration
	
	#initialize Y.tilde (temporarily pretending correlation between observations is zero)
	Y.tilde=array(0,dim=c(Meta$n.species,max(Meta$M),Meta$n.transects))
	for(isp in 1:Meta$n.species){
		for(itrans in 1:Meta$n.transects){
			X=get_mod_matrix(Cur.dat=Data[isp,itrans,,],Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
			ExpY=X%*%Par$det
			Y.tilde[isp,,itrans]=rtruncnorm(max(Meta$M), a=ifelse(Data[isp,itrans,,2]==0,-Inf,0), b=ifelse(Data[isp,itrans,,2]==0,0,Inf), ExpY, 1)		
		}
	}

	#initialize lambda
	Lambda=matrix(0,Meta$n.species,Meta$S)
	Lambda.trans=matrix(0,Meta$n.species,Meta$n.transects)
	for(isp in 1:Meta$n.species){
		Lambda[isp,]=exp(Par$Nu[isp,])*Meta$Area.hab
		Lambda.trans[isp,]=Lambda[isp,Meta$Mapping]*Meta$Area.trans
	}
	grp.lam=rep(0,Meta$n.species)
	
	
	#initialize statistics/matrices needed for MCMC updates
	XpXinv.hab=vector('list',Meta$n.species)
	XpXinvXp.hab=XpXinv.hab
	for(isp in 1:Meta$n.species){
		XpXinv.hab[[isp]]=solve(crossprod(DM.hab[[isp]]))
		XpXinvXp.hab[[isp]]=XpXinv.hab[[isp]]%*%t(DM.hab[[isp]])
	}

	if(Meta$srr){
		L.t=XpXinv.hab
		L=L.t
		Qt=L.t
		cross.L=L.t
		Theta=L.t
		N.theta=rep(0,Meta$n.species)
		for(isp in 1:Meta$n.species){
			P.c=diag(Meta$S)-DM.hab[[isp]]%*%solve(crossprod(DM.hab[[isp]]),t(DM.hab[[isp]]))
			Omega=(P.c%*%Meta$Adj%*%P.c)*(Meta$S/sum(Meta$Adj))
			Eigen=eigen(Omega)
			if(max(Eigen$values)<Meta$srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
			Ind=which(Eigen$values>Meta$srr.tol)
			L.t[[isp]]=Eigen$vectors[,Ind]
			cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
			L[[isp]]=t(L.t[[isp]])
			Qt[[isp]]=L[[isp]]%*%Q%*%L.t[[isp]]
			cross.L[[isp]]=L[[isp]]%*%L.t[[isp]]
			N.theta[isp]=nrow(Qt[[isp]])
			Theta[[isp]]=rnorm(N.theta[isp],0,sqrt(1/Par$tau.eta[isp]))
		}
	}
	Sampled=unique(Meta$Mapping)
	n.unique=length(Sampled)
	Sampled.area.by.strata=rep(0,n.unique)
	for(i in 1:Meta$n.transects)Sampled.area.by.strata[Sampled==Meta$Mapping[i]]=Sampled.area.by.strata[which(Sampled==Meta$Mapping[i])]+Meta$Area.trans[i]
	
	#initialize MCMC, Acceptance rate matrices
	mcmc.length=(Control$iter-Control$burnin)/Control$thin
	MCMC=list(MisID=array(0,dim=c(mcmc.length,dim(Par$MisID))),N.tot=matrix(0,Meta$n.species,mcmc.length),N=array(0,dim=c(Meta$n.species,mcmc.length,Meta$S)),G=array(0,dim=c(Meta$n.species,mcmc.length,Meta$S)),Hab=array(0,dim=c(Meta$n.species,mcmc.length,ncol(Par$hab))),Det=data.frame(matrix(0,mcmc.length,length(Par$det))),cor=rep(0,mcmc.length),tau.eta=matrix(0,Meta$n.species,mcmc.length),tau.nu=matrix(0,Meta$n.species,mcmc.length),Cov.par=array(0,dim=c(Meta$n.species,mcmc.length,length(Par$Cov.par[1,,]))))
	#colnames(MCMC$Hab)=colnames(DM.hab)
	colnames(MCMC$Det)=colnames(DM.det)
	if(Meta$misID==TRUE)Accept=list(cor=0,N=matrix(0,Meta$n.species,Meta$n.transects),Nu=rep(0,Meta$n.species),MisID=matrix(0,dim(Par$MisID)[1],dim(Par$MisID)[2]))
	if(Meta$misID==FALSE)Accept=list(cor=0,N=matrix(0,Meta$n.species,Meta$n.transects),Nu=rep(0,Meta$n.species))
	Pred.N=array(0,dim=c(Meta$n.species,mcmc.length,Meta$n.transects))
	Obs.N=Pred.N

	#check that intial misID parameters are plausible
	if(Meta$misID){
		for(isp in 1:Meta$n.species){
			Cur.dat=matrix(0,sum(G.obs[isp,]*Meta$n.Observers),dim(Data)[4])
			ipl=1
			for(itrans in 1:Meta$n.transects){
				if(G.obs[isp,itrans]>0)Cur.dat[ipl:(ipl+G.obs[isp,itrans]*Meta$n.Observers[itrans]-1),]=Data[isp,itrans,1:(G.obs[isp,itrans]*Meta$n.Observers[itrans]),]
				ipl=ipl+G.obs[isp,itrans]*Meta$n.Observers[itrans]
			}
			Cur.dat=Cur.dat[-which(Cur.dat[,3]==0),]
			DM=vector('list',ncol(Meta$misID.mat))
			XBeta=matrix(0,nrow(Cur.dat),ncol(Meta$misID.mat))	
			for(ipl in 1:ncol(Meta$misID.mat)){  #set up initial parameter values, cell probabilities
				if(Meta$misID.mat[isp,ipl]>0){
					DM[[ipl]]=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[Meta$misID.mat[isp,ipl]]],Levels=Meta$Levels)
					XBeta[,ipl]=exp(DM[[ipl]]%*%Par$MisID[Meta$misID.mat[isp,ipl],1:Meta$N.par.misID[Meta$misID.mat[isp,ipl]]])
				}
			}
			if(sum(XBeta[,isp]<0 | XBeta[,isp]<apply(XBeta[,-isp],1,'max'))!=0){ 
				cat('\n Error: Initial misidentifiation parameters are not consistent with assumption that correct ID probability is > misID probability for all animals \n')#implies XB bigger than other ID possibilities for all animals
			}
		}
	}
	
	#initialize random effect matrices for individual covariates if required
	if(sum(1-Meta$Cov.prior.fixed)>0)RE.cov=array(0,dim=c(Meta$n.species,Meta$n.transects,max(Meta$M),Meta$n.ind.cov))
		
	PROFILE=FALSE
	st <- Sys.time()
	##################################################
	############   Begin MCMC Algorithm ##############
	##################################################
	for(iiter in 1:cur.iter){
		#cat(paste('\n ', iiter))
		for(isp in 1:Meta$n.species){		
		########## update abundance parameters at the strata scale   ################
		
		#update nu parameters (log lambda)
		#1) for sampled cells
			Hab=Par$hab[isp,1:Meta$N.hab.par[isp]]
			Eta=Par$Eta[isp,]
			Mu=DM.hab[[isp]]%*%Hab+Eta
			G.sampled=rep(0,n.unique) #total number of groups currently in each sampled strata
			for(i in 1:Meta$n.transects)G.sampled[Sampled==Meta$Mapping[i]]=G.sampled[Sampled==Meta$Mapping[i]]+Meta$G.transect[isp,i]
			Grad1=log_lambda_gradient(Mu=Mu,Nu=Par$Nu[isp,],Sampled=Sampled,Area=Sampled.area.by.strata,N=G.sampled,var.nu=1/Par$tau.nu[isp])
			Prop=Par$Nu[isp,]
			Prop[Sampled]=Par$Nu[isp,Sampled]+Control$MH.nu[isp]^2*0.5*Grad1+Control$MH.nu[isp]*rnorm(n.unique)
			new.post=log_lambda_log_likelihood(Log.lambda=Prop[Sampled],DM=DM.hab[[isp]],Beta=Hab,SD=sqrt(1/Par$tau.nu[isp]),N=G.sampled,Sampled=Sampled,Area=Sampled.area.by.strata)
			old.post=log_lambda_log_likelihood(Log.lambda=Par$Nu[isp,Sampled],DM=DM.hab[[isp]],Beta=Hab,SD=sqrt(1/Par$tau.nu[isp]),N=G.sampled,Sampled=Sampled,Area=Sampled.area.by.strata)
			Grad2=log_lambda_gradient(Mu=Mu,Nu=Prop,Sampled=Sampled,Area=Sampled.area.by.strata,N=G.sampled,var.nu=1/Par$tau.nu[isp])
			diff1=as.vector(Par$Nu[isp,Sampled]-Prop[Sampled]-0.5*Control$MH.nu[isp]^2*Grad2)	
			diff2=as.vector(Prop[Sampled]-Par$Nu[isp,Sampled]-0.5*Control$MH.nu[isp]^2*Grad1)
			log.jump=0.5/Control$MH.nu[isp]^2*(sqrt(crossprod(diff1,diff1))-sqrt(crossprod(diff2,diff2))) #ratio of jumping distributions using e.g. Robert and Casella 2004 p. 319	
			#cat(paste("\n iter=",iiter," new=",new.post," old=",old.post," jump=",log.jump))
			if(runif(1)<exp(new.post-old.post+log.jump)){
				Par$Nu[isp,]=Prop
				Accept$Nu[isp]=Accept$Nu[isp]+1	
			}
			#2) simulate nu for areas not sampled
#			for(ipar in 1:n.unique){
#				cur=Par$Nu[isp,Sampled[ipar]]
#				prop=cur+runif(1,-1,1)
#				old.logL=dnorm(cur,Mu[Sampled[ipar]],sqrt(1/Par$tau.nu[isp]),log=TRUE)+dpois(Meta$G.transect[1,ipar],0.25*exp(cur),log=1)
#				new.logL=dnorm(prop,Mu[Sampled[ipar]],sqrt(1/Par$tau.nu[isp]),log=TRUE)+dpois(Meta$G.transect[1,ipar],0.25*exp(prop),log=1)
#				if(runif(1)<exp(new.logL-old.logL))Par$Nu[isp,Sampled[ipar]]=prop
#			}
			Par$Nu[isp,-Sampled]=rnorm(Meta$S-n.unique,Mu[-Sampled],1/sqrt(Par$tau.nu[isp]))
		
		
			if(Meta$spat.ind==FALSE){
				if(Meta$srr==FALSE){
					#update eta parameters (spatial random effects)
					V.eta.inv <- Par$tau.nu[isp]*diag(Meta$S) + Par$tau.eta[isp]*Q
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*(Par$Nu[isp,]-DM.hab[[isp]]%*%Hab))		
					Par$Eta[isp,]<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
					Par$Eta[isp,]=Par$Eta[isp,]-mean(Par$Eta[isp,])  #centering
					
					#update tau_eta  (precision of spatial process)
					Par$tau.eta[isp] <- rgamma(1, (Meta$S-1)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta[isp,], Q %*% Par$Eta[isp,])*0.5) + Prior.pars$b.eta)
				}
				else{
					#Update Theta
					Dat.minus.Exp=Par$Nu[isp,]-DM.hab[[isp]]%*%Hab
					V.eta.inv <- cross.L[[isp]]*Par$tau.nu[isp] + Par$tau.eta[isp]*Qt[[isp]]
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*L[[isp]]%*%Dat.minus.Exp)
					Theta <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(N.theta[isp],0,1))
					Par$Eta[isp,]=as.numeric(L.t[[isp]]%*%Theta)
					
					#update tau.eta
					Par$tau.eta[isp] <- rgamma(1, N.theta[isp]*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta, Qt[[isp]] %*% Theta)*0.5) + Prior.pars$b.eta)
				}
			}
			if(PROFILE==TRUE){
				cat(paste("Nu: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			#update tau_nu	 (precision for Poisson overdispersion)
			Mu=DM.hab[[isp]]%*%Hab+Par$Eta[isp,]
			if(Meta$fix.tau.nu==FALSE){
				Diff=Par$Nu[isp,Sampled]-Mu[Sampled]
				Par$tau.nu[isp] <- rgamma(1,n.unique*0.5 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.nu)
			}
			if(PROFILE==TRUE){
				cat(paste("Tau nu: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			
			#translate to lambda scale
			Lambda[isp,]=exp(Par$Nu[isp,])*Meta$Area.hab
			Lambda.trans[isp,]=Lambda[isp,Meta$Mapping]*Meta$Area.trans

			#update Betas for habitat relationships
			Hab=rmvnorm(1,XpXinvXp.hab[[isp]]%*%(Par$Nu[isp,]-Par$Eta[isp,]),XpXinv.hab[[isp]]/Par$tau.nu[isp])
			Par$hab[isp,1:Meta$N.hab.par[isp]]=Hab

			
			########## update group abundance at strata level
			Par$G[isp,]=rpois(Meta$S,Lambda[isp,]*(1-Meta$Covered.area))
			grp.lam[isp]=ifelse(Meta$Cov.prior.pdf[isp,1] %in% c("pois1_ln","poisson_ln"),exp(Par$Cov.par[isp,1,1]+(Par$Cov.par[isp,2,1])^2/2),Par$Cov.par[isp,1,1])
			Par$N[isp,]=Par$G[isp,]+rpois(Meta$S,grp.lam[isp]*Par$G[isp,]) #add the Par$G since number in group > 0 
			Par$G[isp,Meta$Mapping]=Par$G[isp,Meta$Mapping]+Meta$G.transect[isp,]
			Par$N[isp,Meta$Mapping]=Par$N[isp,Meta$Mapping]+Meta$N.transect[isp,]
		}
		
		if(adapt==TRUE){
			if(iiter%%100==0){
				if(Accept$cor<30)Control$MH.cor=Control$MH.cor*.95
				if(Accept$cor>40)Control$MH.cor=Control$MH.cor*1.053
				if(Meta$misID){
					for(irow in 1:nrow(Par$MisID)){
						for(icol in 1:ncol(Par$MisID)){
							if(Accept$MisID[irow,icol]<30)Control$MH.misID[irow,icol]=Control$MH.misID[irow,icol]*0.95
							if(Accept$MisID[irow,icol]>40)Control$MH.misID[irow,icol]=Control$MH.misID[irow,icol]*1.053
						}
					}
				}
				for(ipar in 1:Meta$n.species){
					if(Accept$Nu[ipar]<55)Control$MH.nu[ipar]=Control$MH.nu[ipar]*.95
					if(Accept$Nu[ipar]>60)Control$MH.nu[ipar]=Control$MH.nu[ipar]*1.053
				}
				for(ipar in 1:length(Accept$N)){
					if(Accept$N[ipar]<30)Control$RJ.N[ipar]=max(1,Control$RJ.N[ipar]-1)
					if(Accept$N[ipar]>40)Control$RJ.N[ipar]=Control$RJ.N[ipar]+1
				}
				Accept$cor=0
				Accept$Hab=Accept$Hab*0
				Accept$Nu=Accept$Nu*0
				Accept$N=Accept$N*0
				Accept$MisID=Accept$MisID*0
			}
		}
		
		
		#store results if applicable
		if(iiter>Control$burnin & iiter%%Control$thin==0){
			MCMC$cor[(iiter-Control$burnin)/Control$thin]=Par$cor
			MCMC$Det[(iiter-Control$burnin)/Control$thin,]=Par$det
			if(Meta$misID)MCMC$MisID[(iiter-Control$burnin)/Control$thin,,]=Par$MisID
			for(isp in 1:Meta$n.species){
				MCMC$G[isp,(iiter-Control$burnin)/Control$thin,]=Par$G[isp,]
				MCMC$N[isp,(iiter-Control$burnin)/Control$thin,]=Par$N[isp,]
				MCMC$N.tot[isp,(iiter-Control$burnin)/Control$thin]=sum(Par$N[isp,])
				MCMC$Hab[isp,(iiter-Control$burnin)/Control$thin,]=Par$hab[isp,]
				MCMC$tau.eta[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.eta[isp]
				MCMC$tau.nu[isp,(iiter-Control$burnin)/Control$thin]=Par$tau.nu[isp]
				MCMC$Cov.par[isp,(iiter-Control$burnin)/Control$thin,]=Par$Cov.par[isp,,]
				Obs.N[isp,(iiter-Control$burnin)/Control$thin,]=Meta$N.transect[isp,]
				Temp.G=Meta$Area.hab[Meta$Mapping]*Meta$Area.trans*exp(rnorm(Meta$n.transects,(DM.hab[[isp]]%*%Par$hab[isp,1:Meta$N.hab.par[isp]]+Par$Eta[isp,])[Meta$Mapping],sqrt(1/Par$tau.nu[isp])))
				Pred.N[isp,(iiter-Control$burnin)/Control$thin,]=Temp.G+rpois(Meta$n.transects,grp.lam[isp]*Temp.G)	
			}
		}
		
		if(iiter==100){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/100
			ttc <- round((cur.iter-100)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}	
		if((iiter%%1000)==1)cat(paste('iteration ', iiter,' of ',cur.iter,' completed \n'))
	}
	cat(paste('\n total elapsed time: ',difftime(Sys.time(),st,units="mins"),' minutes \n'))

	Post=list(N=MCMC$N,G=MCMC$G)
	#convert Out$MCMC into mcmc object for use with coda, S3 methods
	Hab.names=vector("list",Meta$n.species)
	for(isp in 1:Meta$n.species)Hab.names[[isp]]=colnames(DM.hab[[isp]])
	Det.names=colnames(get_mod_matrix(Cur.dat=matrix(Data[1,1,1,],nrow=1),Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels))
	Cov.names=vector("list",Meta$n.ind.cov)
	Cov.par.n=0
	if(Meta$n.ind.cov>0){
		for(icov in 1:Meta$n.ind.cov){
			Par.name=switch(Meta$Cov.prior.pdf[icov],pois1_ln=c("mean.minus.1","sd"),poisson_ln=c("mean","sd"),multinom=paste("prop.cell.",c(1:(Meta$Cov.prior.n[isp,icov]-1)),sep=''),normal="mean",pois1="mean.minus.1",poisson="mean")
			Cov.names[[icov]]=paste(Meta$stacked.names[Meta$dist.pl+icov],".",Par.name,sep='')
			Cov.par.n=Cov.par.n+length(Cov.names[[icov]])
		}
	}
	MisID.names=NULL
	if(Meta$misID==TRUE){
		MisID.names=vector("list",max(Meta$misID.mat))
		for(imod in 1:max(Meta$misID.mat)){
			MisID.names[[imod]]=colnames(get_mod_matrix(Cur.dat=matrix(Data[1,1,1,],nrow=1),stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[imod]],Levels=Meta$Levels))		
		}
	}																																				
	MCMC=convert.HDS.to.mcmc(MCMC=MCMC,N.hab.par=Meta$N.hab.par,Cov.par.n=Cov.par.n,Hab.names=Hab.names,Det.names=Det.names,Cov.names=Cov.names,MisID.names=MisID.names,N.par.misID=Meta$N.par.misID,misID.mat=Meta$misID.mat,fix.tau.nu=Meta$fix.tau.nu,misID=Meta$misID,spat.ind=Meta$spat.ind,point.ind=Meta$point.ind)
	Out=list(Post=Post,MCMC=MCMC,Accept=Accept,Control=Control,Obs.N=Obs.N,Pred.N=Pred.N)
	Out
}

