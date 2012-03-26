#' Function for MCMC analysis 
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		"det": a vector giving the current iteration's linear model parameters for the detection model;
#' 		"hab": a vector giving the current iteration's linear model parameters for abundance intensity;
#' 		"cor": a correlation parameter for detections that's an increasing function of distance (correlation at the maximum distance);
#' 		"Nu": a vector giving the log of the abundance intensity for each strata;
#' 		"G": a vector giving the number of groups of animals in each strata; 
#' 		"N": a vector giving the number of animals in each strata
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
#'  "Cov.prior.parms"	(n x n.ind.cov) matrix providing "pseudo-prior" parameters for individual covarate distributions (only the first row used if a signle parameter distribution)
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
			Par$Nu[isp,-Sampled]=rnorm(Meta$S-n.unique,Mu[-Sampled],1/sqrt(Par$tau.nu[isp]))
		
		
			if(Meta$spat.ind==FALSE){
				if(Meta$srr==FALSE){
					#update eta parameters (spatial random effects)
					V.eta.inv <- Par$tau.nu[isp]*diag(Meta$S) + Par$tau.eta[isp]*Q
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*(Par$Nu[isp,]-DM.hab[[isp]]%*%Hab))		
					Par$Eta[isp,]<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
					Par$Eta[isp,]=Par$Eta[isp,]-mean(Par$Eta[isp,])  #centering
					
					#update tau_eta  (precision of spatial process)
					Par$tau.eta[isp] <- rgamma(1, (Meta$S-1)/2 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta[isp,], Q %*% Par$Eta[isp,])/2) + Prior.pars$b.eta)
				}
				else{
					#Update Theta
					Dat.minus.Exp=Par$Nu[isp,]-DM.hab[[isp]]%*%Hab
					V.eta.inv <- cross.L[[isp]]*Par$tau.nu[isp] + Par$tau.eta[isp]*Qt[[isp]]
					M.eta <- solve(V.eta.inv, Par$tau.nu[isp]*L[[isp]]%*%Dat.minus.Exp)
					Theta <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(N.theta[isp],0,1))
					Par$Eta[isp,]=as.numeric(L.t[[isp]]%*%Theta)
					
					#update tau.eta
					Par$tau.eta[isp] <- rgamma(1, N.theta[isp]/2 + Prior.pars$a.eta, as.numeric(crossprod(Theta, Qt[[isp]] %*% Theta)/2) + Prior.pars$b.eta)
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
				Par$tau.nu[isp] <- rgamma(1,n.unique/2 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))/2 + Prior.pars$b.nu)
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

			if(PROFILE==TRUE){
				cat(paste("Hab, etc.: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			########## update abundance, distances, ind. covariates for observed transects using RJMCMC  #############
			
			for(itrans in 1:Meta$n.transects){
				Sample=c(-Control$RJ.N[isp,itrans]:Control$RJ.N[isp,itrans])
				Sample=Sample[-(Control$RJ.N[isp,itrans]+1)] #don't make proposals where pop size stays the same
				a=sample(Sample,1)
				Sigma=matrix(Par$cor,Meta$n.Observers[itrans],Meta$n.Observers[itrans])
				diag(Sigma)=1
				offdiag=which(Sigma!=1)
				
				if(a>0){ # proposal an addition
					if(((Meta$G.transect[isp,itrans]+a)*Meta$n.Observers[itrans])>=Meta$M[isp,itrans]){
						if(adapt==FALSE)cat(paste('\n Warning: proposed abundance for transect ',itrans,' species ',isp, '> M; consider increasing M value! \n'))
						else{
							temp=floor(Meta$M[isp,itrans]*1.25)
							if(temp%%2==1)temp=temp+1
							Meta$M[isp,itrans]=min(temp,max(Meta$M[isp,]))
						}
					}
					else{
						Cur.dat=Data[isp,itrans,(n.Records[isp,itrans]+1):(n.Records[isp,itrans]+a*Meta$n.Observers[itrans]),]
						if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
						X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$Det.formula,Levels=Meta$Levels)
						ExpY=X%*%Par$det
						P=c(1:a)
						for(i in 1:a){
							Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
							P[i]=pmvnorm(upper=rep(0,Meta$n.Observers[itrans]),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
						}
						tmp.sum=0
						for(i in 1:a)tmp.sum=tmp.sum-log(Meta$G.transect[isp,itrans]-G.obs[isp,itrans]+i)+log(P[i])
						MH.prob=exp(a*log(Lambda.trans[isp,itrans])+tmp.sum)
						if(runif(1)<MH.prob){
							Meta$G.transect[isp,itrans]=Meta$G.transect[isp,itrans]+a
							n.Records[isp,itrans]=Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]
							if(Meta$grps==FALSE)Meta$N.transect[isp,itrans]=Meta$G.transect[isp,itrans]
							else Meta$N.transect[isp,itrans]=sum(Data[isp,itrans,1:n.Records[isp,itrans],Meta$stacked.names=="Group"])/Meta$n.Observers[itrans]
							Accept$N[isp,itrans]=Accept$N[isp,itrans]+1
							#generate Y-tilde values
							Tmp=matrix(ExpY,a,Meta$n.Observers[itrans],byrow=TRUE)
							Y.tilde.temp=rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,1],sd=1)
							if(Meta$n.Observers[itrans]>1){
								Dist=matrix(Cur.dat[,Meta$dist.pl],a,Meta$n.Observers[itrans],byrow=TRUE)
								Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
								Y.tilde.temp=cbind(Y.tilde.temp,rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,2]+Cor*(Y.tilde.temp-Tmp[,1]),sd=sqrt(1-Cor^2)))
							}
							Y.tilde[isp,(n.Records[isp,itrans]-a*Meta$n.Observers[itrans]+1):n.Records[isp,itrans],itrans]=as.vector(t(Y.tilde.temp))
						}	
					}
				}
				else{ #proposal a deletion
					if((Meta$G.transect[isp,itrans]+a)>=G.obs[isp,itrans]){ #can only delete records where an animal isn't observed
						Cur.dat=Data[isp,itrans,n.Records[isp,itrans]:(n.Records[isp,itrans]+a*Meta$n.Observers[itrans]+1),]
						if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
						X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
						ExpY=X%*%Par$det
						P=c(1:-a)
						for(i in 1:-a){
							Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
							P[i]=pmvnorm(upper=rep(0,Meta$n.Observers[itrans]),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
						}
						tmp.sum=0
						for(i in 1:-a){
							tmp.sum=tmp.sum+log(Meta$G.transect[isp,itrans]-G.obs[isp,itrans]-i+1)-log(P[i])
						}
						MH.prob=exp(a*log(Lambda.trans[isp,itrans])+tmp.sum)
						if(runif(1)<MH.prob){
							Meta$G.transect[isp,itrans]=Meta$G.transect[isp,itrans]+a
							n.Records[isp,itrans]=Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]
							Accept$N[isp,itrans]=Accept$N[isp,itrans]+1
						}					
						
					}
				}
				
				#update distances, individual covariates for latent animals not currently in the population
				#fill distances for unobserved
				if(Meta$i.binned==1)Data[isp,itrans,(n.Records[isp,itrans]+1):Meta$M[isp,itrans],Meta$dist.pl]=rep(sample(c(1:Meta$n.bins),size=(Meta$M[isp,itrans]-n.Records[isp,itrans])/Meta$n.Observers[itrans],replace=TRUE,prob=Meta$Bin.length),each=Meta$n.Observers[itrans])
				else Data[isp,itrans,(n.Records[isp,itrans]+1):Meta$M[isp,itrans],Meta$dist.pl]=rep(runif((Meta$M[isp,itrans]-n.Records[isp,itrans])/Meta$n.Observers[itrans]),each=Meta$n.Observers[itrans])
				#fill individual covariate values for (potential) animals that weren't observed
				if(Meta$n.ind.cov>0){
					for(icov in 1:Meta$n.ind.cov){
						if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=RE.cov[isp,itrans,(Meta$G.transect[isp,itrans]+1):(Meta$M[isp,itrans]/Meta$n.Observers[itrans]),icov]
						else cur.RE=0
						rsamp=switch_sample(n=(Meta$M[isp,itrans]-n.Records[isp,itrans])/Meta$n.Observers[itrans],pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
						Data[isp,itrans,(n.Records[isp,itrans]+1):Meta$M[isp,itrans],Meta$dist.pl+icov]=rep(rsamp,each=Meta$n.Observers[itrans])
					}
				}
				
				#update distances, individual covariates for animals that ARE in the population but never observed
				cur.G=Meta$G.transect[isp,itrans]-G.obs[isp,itrans]
				if(cur.G>0){
					#distance
					if(Meta$i.binned==1)dist.star=sample(c(1:Meta$n.bins),cur.G,replace=TRUE,prob=Meta$Bin.length)
					else dist.star=runif(cur.G)
					Cur.dat=Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),]
					if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
					cur.dist=Cur.dat[,Meta$dist.pl]
					X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY=X%*%Par$det
					L.old=c(1:length(dist.star))
					Tmp.Y.tilde=Y.tilde[isp,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),itrans]
					for(i in 1:length(dist.star)){
						Sigma[offdiag]=Par$cor*(Meta$i.binned*(cur.dist[i]-1)*dist.mult+(1-Meta$i.binned)*cur.dist[i])
						L.old[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
					}	
					Cur.dat[,Meta$dist.pl]=rep(dist.star,each=Meta$n.Observers[itrans])
					X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY=X%*%Par$det
					L.star=L.old
					for(i in 1:length(dist.star)){
						Sigma[offdiag]=Par$cor*(Meta$i.binned*(dist.star[i]-1)*dist.mult+(1-Meta$i.binned)*dist.star[i])
						L.star[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
					}	
					Acc=(runif(length(L.star))<(L.star/L.old))
					Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl]+rep(Acc,each=Meta$n.Observers[itrans])*rep(dist.star,each=Meta$n.Observers[itrans])
					
					#individual covariates
					if(Meta$n.ind.cov>0){
						for(icov in 1:Meta$n.ind.cov){
							if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=RE.cov[isp,itrans,(G.obs[isp,itrans]+1):Meta$G.transect[isp,itrans],icov]
							else cur.RE=0
							Cov.star=switch_sample(n=cur.G,pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
							Cur.dat=Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),]
							if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
							cur.dist=Cur.dat[,Meta$dist.pl]
							X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
							ExpY=X%*%Par$det
							L.old=c(1:length(Cov.star))
							for(i in 1:length(Cov.star)){
								Sigma[offdiag]=Par$cor*(Meta$i.binned*(cur.dist[i]-1)*dist.mult+(1-Meta$i.binned)*cur.dist[i])
								L.old[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
							}	
							Cur.dat[,Meta$dist.pl+icov]=rep(Cov.star,each=Meta$n.Observers[itrans])
							X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
							ExpY=X%*%Par$det
							L.star=L.old
							for(i in 1:length(Cov.star)){
								Sigma[offdiag]=Par$cor*(Meta$i.binned*(cur.dist[i]-1)*dist.mult+(1-Meta$i.binned)*cur.dist[i])
								L.star[i]=dmvnorm(Tmp.Y.tilde[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
							}	
							Acc=(runif(length(L.star))<(L.star/L.old))
							Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[isp,itrans,(G.obs[isp,itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[isp,itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]+rep(Acc,each=Meta$n.Observers[itrans])*rep(Cov.star,each=Meta$n.Observers[itrans])						
						}
					}
				}
				if(Meta$grps==FALSE)Meta$N.transect[isp,itrans]=Meta$G.transect[isp,itrans]
				else Meta$N.transect[isp,itrans]=sum(Data[isp,itrans,1:n.Records[isp,itrans],Meta$stacked.names=="Group"])/Meta$n.Observers[itrans]

			
				#update Y.tilde
				if(n.Records[isp,itrans]>0){
					Cur.dat=Data[isp,itrans,1:n.Records[isp,itrans],]
					if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
					X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY=matrix(X%*%Par$det,Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
					Temp.Y.tilde=matrix(Y.tilde[isp,1:n.Records[isp,itrans],itrans],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
					if(Meta$n.Observers[itrans]==1){
						Temp.Y.tilde<-rtruncnorm(Meta$G.transect[isp,itrans],a=ifelse(Cur.dat[,2]==0,-Inf,0),b=ifelse(Cur.dat[,2]==0,0,Inf),ExpY,1)
					}
					else{
						Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
						Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
						Resp=matrix(Cur.dat[,2],Meta$G.transect[isp,itrans],Meta$n.Observers[itrans],byrow=TRUE)
						EY1=ExpY[,1]+Cor*(Temp.Y.tilde[,2]-ExpY[,2])
						Temp.Y.tilde[,1] <- rtruncnorm(Meta$G.transect[isp,itrans], a=ifelse(Resp[,1]==0,-Inf,0), b=ifelse(Resp[,1]==0,0,Inf), EY1, sqrt(1-Cor^2))
						EY2=ExpY[,2]+Cor*(Temp.Y.tilde[,1]-ExpY[,1])
						Temp.Y.tilde[,2] <- rtruncnorm(Meta$G.transect[isp,itrans], a=ifelse(Resp[,2]==0,-Inf,0), b=ifelse(Resp[,2]==0,0,Inf), EY2, sqrt(1-Cor^2))
					}
					Y.tilde[isp,1:n.Records[isp,itrans],itrans]=as.vector(t(Temp.Y.tilde))
				}
			}
			if(PROFILE==TRUE){
				cat(paste("Data aug: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			
		}
		
		##### update true species for observed animals ######
		if(Meta$misID){
			for(isamp in 1:n.samp.misID){
				isp=sample(Meta$n.species,1,prob=apply(G.obs,1,'sum'))
				itrans=sample(Meta$n.transects,1,prob=G.obs[isp,])
				iind=sample(G.obs[isp,itrans],1)
				Cur.dat=Data[isp,itrans,iind*Meta$n.Observers[itrans]+(1-Meta$n.Observers[itrans]):0,]
				if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
				dist=Cur.dat[1,Meta$dist.pl]
				cor=Par$cor*(Meta$i.binned*(dist-1)*dist.mult+(1-Meta$i.binned)*dist)
				Cur.Y.tilde=Y.tilde[isp,iind*Meta$n.Observers[itrans]+(1-Meta$n.Observers[itrans]):0,itrans]
				X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
				ExpY=X%*%Par$det
				Other.sp=c(1:Meta$n.species)[-isp]
				if(length(Other.sp)==1)prop.sp=Other.sp
				else prop.sp=sample(Other.sp,1)
				New.dat=Cur.dat
				New.dat[,4]=prop.sp
				X=get_mod_matrix(Cur.dat=New.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
				ExpY.prop=X%*%Par$det
				Obs.species=Cur.dat[,3]
				#abundance component
				mh=log(Lambda.trans[prop.sp,itrans])+log(G.obs[isp,itrans])-log(Lambda.trans[isp,itrans])-log(G.obs[prop.sp,itrans]+1)  
				#detection (Y-tilde) component
				if(Meta$n.Observers[itrans]==1){
					Y.tilde.llik.old=dnorm(Cur.Y.tilde,ExpY,1,log=TRUE)
					Y.tilde.llik.new=dnorm(Cur.Y.tilde,ExpY.prop,1,log=TRUE)
				}
				else{
					Sigma[offdiag]=cor
					Y.tilde.llik.old=dmvnorm(Cur.Y.tilde,ExpY,Sigma,log=TRUE)
					Y.tilde.llik.new=dmvnorm(Cur.Y.tilde,ExpY.prop,Sigma,log=TRUE)						
				}
				mh=mh+Y.tilde.llik.new-Y.tilde.llik.old
				#individual covariate component
				for(icov in 1:Meta$n.ind.cov){
					if(Meta$Cov.prior.pdf[isp,icov]=='poisson_ln' | Meta$Cov.prior.pdf[isp,icov]=='pois1_ln')cur.RE=RE.cov[c(isp,prop.sp),itrans,iind*Meta$n.Observers[itrans],icov]
					else cur.RE=c(0,0)
					cur.cov=Cur.dat[1,Meta$dist.pl+icov]
					mh=mh+switch_pdf(x=cur.cov,pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Meta$Cov.prior.parms[prop.sp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)								
					mh=mh-switch_pdf(x=cur.cov,pdf=Meta$Cov.prior.pdf[isp,icov],cur.par=Meta$Cov.prior.parms[isp,1:Meta$Cov.prior.n[isp,icov],icov],RE=cur.RE)
				}							
				#misID component
				#old species
				Conf=matrix(0,Meta$n.Observers[itrans],ncol(Meta$misID.mat))
				for(ipl in 1:ncol(Meta$misID.mat)){
					if(Meta$misID.mat[isp,ipl]>0){
						DM=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[Meta$misID.mat[isp,ipl]]],Levels=Meta$Levels)
						Conf[,ipl]=exp(DM%*%Par$MisID[Meta$misID.mat[isp,ipl],1:ncol(DM)])
					}
				}
				Sum=apply(Conf,1,'sum')
				Conf=Conf/(1+Sum)
				Conf[,Meta$misID.mat[isp,]<0]=1-apply(Conf,1,'sum')	
				for(iobs in 1:Meta$n.Observers[itrans])if(Obs.species[iobs]>0)mh=mh-log(Conf[iobs,Obs.species[iobs]])
				#new species
				Conf=matrix(0,Meta$n.Observers[itrans],ncol(Meta$misID.mat))
				Cur.dat[,4]=prop.sp
				for(ipl in 1:ncol(Meta$misID.mat)){
					if(Meta$misID.mat[prop.sp,ipl]>0){
						DM=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[Meta$misID.mat[prop.sp,ipl]]],Levels=Meta$Levels)
						Conf[,ipl]=exp(DM%*%Par$MisID[Meta$misID.mat[prop.sp,ipl],1:ncol(DM)])
					}
				}
				Sum=apply(Conf,1,'sum')
				Conf=Conf/(1+Sum)
				Conf[,Meta$misID.mat[prop.sp,]<0]=1-apply(Conf,1,'sum')	
				for(iobs in 1:Meta$n.Observers[itrans])if(Obs.species[iobs]>0)mh=mh+log(Conf[iobs,Obs.species[iobs]])
				
				if(runif(1)<exp(mh)){  #accept species change!
					if((n.Records[prop.sp,itrans]+Meta$n.Observers[itrans])<Meta$M[prop.sp,itrans]){						
						#cat(paste("species changed; ind ",iind))
						#1) add data to new target species' data aug matrix
						Data[prop.sp,itrans,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+Meta$n.Observers[itrans]+1):Meta$M[prop.sp,itrans],]=Data[prop.sp,itrans,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1):(Meta$M[prop.sp,itrans]-Meta$n.Observers[itrans]),]	
						Data[prop.sp,itrans,G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1:Meta$n.Observers[itrans],]=Cur.dat
						#2) remove entry from old species; shift data aug array down
						Data[isp,itrans,((iind-1)*(Meta$n.Observers[itrans])+1):(Meta$M[isp,itrans]-Meta$n.Observers[itrans]),]=Data[isp,itrans,((iind-1)*(Meta$n.Observers[itrans])+1+Meta$n.Observers[itrans]):Meta$M[isp,itrans],]
						#3) add Y.tilde to new target species Y.tilde matrix
						Y.tilde[prop.sp,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+Meta$n.Observers[itrans]+1):Meta$M[prop.sp,itrans],itrans]=Y.tilde[prop.sp,(G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1):(Meta$M[prop.sp,itrans]-Meta$n.Observers[itrans]),itrans]	
						Y.tilde[prop.sp,G.obs[prop.sp,itrans]*Meta$n.Observers[itrans]+1:Meta$n.Observers[itrans],itrans]=Cur.Y.tilde	
						#4) remove Y.tilde from current species' Y.tilde matrix
						Y.tilde[isp,((iind-1)*(Meta$n.Observers[itrans])+1):(Meta$M[isp,itrans]-Meta$n.Observers[itrans]),itrans]=Y.tilde[isp,((iind-1)*(Meta$n.Observers[itrans])+1+Meta$n.Observers[itrans]):Meta$M[isp,itrans],itrans]
						#5) adjust dimensions of things
						G.obs[isp,itrans]=G.obs[isp,itrans]-1
						Meta$G.transect[isp,itrans]=Meta$G.transect[isp,itrans]-1
						n.Records[isp,itrans]=n.Records[isp,itrans]-Meta$n.Observers[itrans]
						G.obs[prop.sp,itrans]=G.obs[prop.sp,itrans]+1
						Meta$G.transect[prop.sp,itrans]=Meta$G.transect[prop.sp,itrans]+1
						n.Records[prop.sp,itrans]=n.Records[prop.sp,itrans]+Meta$n.Observers[itrans]
					}
				}		
			}
			if(PROFILE==TRUE){
				cat(paste("True species: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			
			##### update misID parameters
			for(isp in 1:Meta$n.species){
				Cur.dat=matrix(0,sum(G.obs[isp,]*Meta$n.Observers),dim(Data)[4])
				ipl=1
				for(itrans in 1:Meta$n.transects){
					if(G.obs[isp,itrans]>0)Cur.dat[ipl:(ipl+G.obs[isp,itrans]*Meta$n.Observers[itrans]-1),]=Data[isp,itrans,1:(G.obs[isp,itrans]*Meta$n.Observers[itrans]),]
					ipl=ipl+G.obs[isp,itrans]*Meta$n.Observers[itrans]
				}
				Cur.dat=Cur.dat[-which(Cur.dat[,3]==0),]
				DM=vector('list',ncol(Meta$misID.mat))
				Conf=matrix(0,nrow(Cur.dat),ncol(Meta$misID.mat))	
				XBeta=Conf
				for(ipl in 1:ncol(Meta$misID.mat)){  #set up initial parameter values, cell probabilities
					if(Meta$misID.mat[isp,ipl]>0){
						DM[[ipl]]=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$misID.models[[Meta$misID.mat[isp,ipl]]],Levels=Meta$Levels)
						XBeta[,ipl]=exp(DM[[ipl]]%*%Par$MisID[Meta$misID.mat[isp,ipl],1:Meta$N.par.misID[Meta$misID.mat[isp,ipl]]])
					}
				}
				Conf=XBeta/(1+apply(XBeta,1,'sum'))
				Conf[,which(Meta$misID.mat[isp,]==-1)]=1-apply(Conf,1,'sum')
				
				Row.num=c(1:nrow(Cur.dat))
				logL.old=sum(log(Conf[(Cur.dat[,3]-1)*nrow(Cur.dat)+Row.num]))  #categorical distribution
				
				for(ipl in 1:ncol(Meta$misID.mat)){  #now, update parameter values
					if(Meta$misID.mat[isp,ipl]>0){
						Temp.par=Par$MisID[Meta$misID.mat[isp,ipl],1:Meta$N.par.misID[Meta$misID.mat[isp,ipl]]]
						for(ipar in 1:Meta$N.par.misID[Meta$misID.mat[isp,ipl]]){
							beta.old=Par$MisID[Meta$misID.mat[isp,ipl],ipar]
							beta.star=beta.old+rnorm(1,0,Control$MH.misID[Meta$misID.mat[isp,ipl],ipar])
							Temp.par[ipar]=beta.star
							XBeta[,ipl]=exp(DM[[ipl]]%*%Temp.par)
							if(sum(XBeta[,isp]<0 | XBeta[,isp]<apply(XBeta[,-isp],1,'max'))==0){ #implies XB bigger than other ID possibilities for all animals
								Conf=XBeta/(1+apply(XBeta,1,'sum'))
								Conf[,which(Meta$misID.mat[isp,]==-1)]=1-apply(Conf,1,'sum')
								logL.new=sum(log(Conf[(Cur.dat[,3]-1)*nrow(Cur.dat)+Row.num]))  #categorical distribution
								if(runif(1)<exp(logL.new-logL.old+dnorm(beta.star,0,100,log=1)-dnorm(beta.old,0,100,log=1))){
									Par$MisID[Meta$misID.mat[isp,ipl],ipar]=beta.star
									Accept$MisID[[Meta$misID.mat[isp,ipl],ipar]]=Accept$MisID[[Meta$misID.mat[isp,ipl],ipar]]+1
									logL.old=logL.new
								}
								else Temp.par[ipar]=beta.old
							}
							else Temp.par[ipar]=beta.old
						}
						XBeta[,ipl]=exp(DM[[ipl]]%*%Temp.par)  #in case last parameter wasn't accepted
					}
				}			
			}
			if(PROFILE==TRUE){
				cat(paste("misID pars: ", (Sys.time()-st),'\n'))
				st=Sys.time()
			}
			
		}
		###############       update detection process parameters       ##############
		# First, assemble stacked adjusted Response, X matrices across all transects; 
		#basic form of response is Ytilde[obs1]-cor*Ytilde[obs2]
		#adjusted DM rows are of form X[obs1]-cor*X[obs2]
		X.beta=matrix(NA,1,n.beta.det)
		Y.beta=NA
		Cor=NA
		for(isp in 1:Meta$n.species){
			GT0=which(Meta$G.transect[isp,]>0)
			n.gt0=length(GT0)
			if(n.gt0>0){
				for(itrans in 1:n.gt0){
					Cur.dat=Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],]
					if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
					Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1]
					X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,Meta$n.Observers[GT0[itrans]],Meta$G.transect[isp,GT0[itrans]]))
					if(Meta$n.Observers[GT0[itrans]]==2){
						Tmp.cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
						Cor=c(Cor,rep(Tmp.cor,2))  #assemble vector of correlation parameters for each observation
						X.beta=rbind(X.beta,t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
						Y.temp=matrix(Y.tilde[isp,1:n.Records[isp,GT0[itrans]],GT0[itrans]],Meta$G.transect[isp,GT0[itrans]],2,byrow=TRUE)
						Y.beta=c(Y.beta,Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
					}
					else{
						X.beta=rbind(X.beta,t(X.temp[,1,]))
						Y.beta=c(Y.beta,Y.tilde[isp,1:n.Records[isp,GT0[itrans]],GT0[itrans]])
						Cor=c(Cor,rep(0,n.Records[isp,GT0[itrans]]))
					}
				}
			}
		}
		X.beta=X.beta[-1,]
		Y.beta=Y.beta[-1]
		Cor=Cor[-1]
		#now use basic matrix equations from Gelman '04 book (14.11 14.12) to update beta parms
		Sig.inv=1/(1-Cor^2)  #for use in eq. 14.11, 14.12 of Gelman et al.; don't need a matrix since it would be diagonal
		V.inv <- crossprod(X.beta,Sig.inv*X.beta) 	
		M.z <- solve(V.inv, crossprod(X.beta,Sig.inv*Y.beta))
		Par$det <- M.z + solve(chol(V.inv), rnorm(n.beta.det,0,1))

		if(PROFILE==TRUE){
			cat(paste("Detection pars: ", (Sys.time()-st),'\n'))
			st=Sys.time()
		}
		
		
		#update correlation parameter for detection process (if applicable)
		if(Meta$point.ind==1){
			cor.star=Par$cor+runif(1,-Control$MH.cor,Control$MH.cor)
			if(cor.star>0 & cor.star<1){
				Delta1=rep(NA,sum(Meta$G.transect))
				Delta2=Delta1
				Dist=Delta1
				counter=1
				for(isp in 1:Meta$n.species){
					I.gt.one=which(Meta$n.Observers>1 & Meta$G.transect[isp,]>0)
					n.gt.one=length(I.gt.one)
					if(n.gt.one>0){
						for(itrans in 1:n.gt.one){
							Cur.dat=Data[isp,I.gt.one[itrans],1:n.Records[isp,I.gt.one[itrans]],]
							if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
							Dist[counter:(counter+Meta$G.transect[isp,I.gt.one[itrans]]-1)]=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[isp,I.gt.one[itrans]],Meta$n.Observers[I.gt.one[itrans]],byrow=TRUE)[,1]
							X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[isp,I.gt.one[itrans]]))
							Y.temp=matrix(Y.tilde[isp,1:n.Records[isp,I.gt.one[itrans]],I.gt.one[itrans]],Meta$G.transect[isp,I.gt.one[itrans]],2,byrow=TRUE)			
							Delta1[counter:(counter+Meta$G.transect[isp,I.gt.one[itrans]]-1)]=Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det)
							Delta2[counter:(counter+Meta$G.transect[isp,I.gt.one[itrans]]-1)]=Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det)
							counter=counter+Meta$G.transect[isp,I.gt.one[itrans]]
						}
					}
				}
				Cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
				logP.old=-.5*(sum(log(1-Cor^2))+sum((Delta1^2+Delta2^2-2*Cor*Delta1*Delta2)/(1-Cor^2)))
				Cor=cor.star*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
				logP.new=-.5*(sum(log(1-Cor^2))+sum((Delta1^2+Delta2^2-2*Cor*Delta1*Delta2)/(1-Cor^2)))
				if(runif(1)<exp(logP.new-logP.old)){
					Par$cor=cor.star
					Accept$cor=Accept$cor+1
				}				
			}
		}

		if(PROFILE==TRUE){
			cat(paste("Correlation: ", (Sys.time()-st),'\n'))
			st=Sys.time()
		}
		
		#update parameters of individual covariate distributions (if fixed=0)
		for(isp in 1:Meta$n.species){
			GT0=which(Meta$G.transect[isp,]>0)
			n.gt0=length(GT0)
			for(icov in 1:Meta$n.ind.cov){
				if(Meta$Cov.prior.fixed[isp,icov]==0){
					if(Meta$Cov.prior.pdf[isp,icov]=="normal")cat("\n Warning: hyper-priors not yet implemented for normal dist. \n")
					if(Meta$Cov.prior.pdf[isp,icov]=="poisson"){
						Cur.cov=matrix(Data[isp,GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
						if(n.gt0>1){
							for(itrans in 2:n.gt0){
								Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
							}
						}
						Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
					}
					if(Meta$Cov.prior.pdf[isp,icov]=="pois1"){
						Cur.cov=matrix(Data[isp,GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
						if(n.gt0>1){
							for(itrans in 2:n.gt0){
								Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
							}
						}
						Par$Cov.par[isp,1,icov]=rgamma(1,sum(Cur.cov)-length(Cur.cov)+Meta$Cov.prior.parms[isp,1,icov],length(Cur.cov)+Meta$Cov.prior.parms[isp,2,icov])
					}
					if(Meta$Cov.prior.pdf[isp,icov]=="poisson_ln" | Meta$Cov.prior.pdf[isp,icov]=="pois1_ln"){
						Cur.cov=matrix(Data[isp,GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
						Cur.RE=RE.cov[isp,GT0[1],1:Meta$G.transect[isp,GT0[1]],icov]
						if(n.gt0>1){
							for(itrans in 2:n.gt0){
								Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
								Cur.RE=c(Cur.RE,RE.cov[isp,GT0[itrans],1:Meta$G.transect[isp,GT0[itrans]],icov])
							}
						}
						Cur.cov=Cur.cov-(Meta$Cov.prior.pdf[isp,icov]=="pois1_ln")
						#1) update theta
						par.star=Par$Cov.par[isp,1,icov]+runif(1,-0.05,0.05)
						sum.y=sum(Cur.cov)
						sum.yZ=sum(Cur.cov*Cur.RE)
						log.post.new=par.star*sum.y-sum(exp(par.star+Par$Cov.par[isp,2,icov]*Cur.RE))+dnorm(par.star,Meta$Cov.prior.parms[isp,1,icov],Meta$Cov.prior.parms[isp,2,icov],log=1)
						log.post.old=Par$Cov.par[isp,1,icov]*sum.y-sum(exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Cur.RE))+dnorm(Par$Cov.par[isp,1,icov],Meta$Cov.prior.parms[isp,1,icov],Meta$Cov.prior.parms[isp,2,icov],log=1)
						if(runif(1)<exp(log.post.new-log.post.old))Par$Cov.par[isp,1,icov]=par.star
						#2) update sigma
						par.star=Par$Cov.par[isp,2,icov]+runif(1,-.01,.01)
						if(par.star>0 & par.star<Meta$Cov.prior.parms[isp,3,icov]){
							log.post.new=par.star*sum.yZ-sum(exp(Par$Cov.par[isp,1,icov]+par.star*Cur.RE))
							log.post.old=Par$Cov.par[isp,2,icov]*sum.yZ-sum(exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Cur.RE))
							if(runif(1)<exp(log.post.new-log.post.old))Par$Cov.par[isp,2,icov]=par.star
						}
						#3) update random effects		
						for(itrans in 1:n.gt0){
							#animals currently in population
							Cur.cov=matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1]-(Meta$Cov.prior.pdf[isp,icov]=="pois1_ln")						
							Cur.RE=RE.cov[isp,GT0[itrans],1:Meta$G.transect[isp,GT0[itrans]],icov]
							Prop=Cur.RE+runif(length(Cur.RE),-.1,.1)
							LogPost.new=Cur.cov*Par$Cov.par[isp,2,icov]*Prop-exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Prop)+dnorm(Prop,0,1,log=1)
							LogPost.old=Cur.cov*Par$Cov.par[isp,2,icov]*Cur.RE-exp(Par$Cov.par[isp,1,icov]+Par$Cov.par[isp,2,icov]*Cur.RE)+dnorm(Cur.RE,0,1,log=1)
							Acc=(runif(length(Cur.RE))<(exp(LogPost.new-LogPost.old)))
							RE.cov[isp,GT0[itrans],1:Meta$G.transect[isp,GT0[itrans]],icov]=Acc*Prop+(1-Acc)*Cur.RE
						}
						#animals currently not in population
						for(itrans in 1:Meta$n.transects){
							RE.cov[isp,itrans,(Meta$G.transect[isp,itrans]+1):Meta$M[isp,itrans],icov]=rnorm(Meta$M[isp,itrans]-Meta$G.transect[isp,itrans],0,1)
						}
					}
					if(Meta$Cov.prior.pdf[isp,icov]=="multinom"){
						Cur.cov=matrix(Data[isp,1:GT0[1],1:n.Records[isp,GT0[1]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
						if(n.gt0>1){
							for(itrans in 2:n.gt0){
								Cur.cov=c(Cur.cov,matrix(Data[isp,GT0[itrans],1:n.Records[isp,GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[isp,GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
							}
						}
						Par$Cov.par[isp,1:Meta$Cov.prior.n[isp,icov],icov]=rdirichlet(1,Meta$Cov.prior.parms[isp,1:Meta$Cov.prior.n[isp,icov],icov]+tabulate(factor(Cur.cov)))
					}
				}
			}
		}
		if(PROFILE==TRUE){
			cat(paste("Ind cov pars: ", (Sys.time()-st),'\n'))
			st=Sys.time()
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

