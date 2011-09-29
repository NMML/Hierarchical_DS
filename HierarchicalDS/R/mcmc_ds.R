#' Function for MCMC analysis 
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		"det": a vector giving the current iteration's linear model parameters for the detection model;
#' 		"hab": a vector giving the current iteration's linear model parameters for abundance intensity;
#' 		"cor": a correlation parameter for detections that's an increasing function of distance (correlation at the maximum distance);
#' 		"Nu": a vector giving the log of the abundance intensity for each strata;
#' 		"G": a vector giving the number of groups of animals in each strata; 
#' 		"N": a vector giving the number of animals in each strata
#' @param Data   A three dimensional array; the first dimension gives the transect, the second dimension indexes a (possible) observation, 
#' 			and the third dimension gives observations and covariates associated with a given animal.
#' 			These final columns are: Observer ID,Y(observation=0/1),Obs covariates,Distance,Ind covariates
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
#'	"MH.beta": A vector of tuning parameters for betas of the abundance process (dimension = number of columns of habitat DM);
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
#' 	"S"				Number of strata cells
#'  "spat.ind"		Indicator for spatial dependence
#'  "Area.hab"		Vector giving relative area covered by each strata
#'  "Area.trans"	Vector giving fraction of area of relevant strata covered by each transect
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
#'  "Cov.prior.pdf" character vector giving the probability density function associated with each individual covariate (current choices: "pois1","poisson","normal","unif.disc","unif.cont")
#'  "Cov.prior.parms"	(2xn.ind.cov) matrix providing "pseudo-prior" parameters for individual covarate distributions (only the first row used if a signle parameter distribution)
#'  "point.ind"		Indicator for whether point independence assumed (if no, then no correlation modeled b/w multiple observers as function of distance)
#' @return returns a list with the following slots: 
#' 	"MCMC": A list object containing posterior samples;
#'  "Accept": A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  "Control": A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used) 
#' @export
#' @import Matrix
#' @keywords areal, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn


mcmc_ds<-function(Par,Data,cur.iter,adapt,Control,DM.hab,DM.det,Q,Prior.pars,Meta){	
	require(mvtnorm)
	require(Matrix)
	require(truncnorm)
	
	Lam.index=c(1:Meta$S)
	if(Meta$i.binned==0)dist.mult=1
	if(Meta$i.binned==1)dist.mult=1/(Meta$n.bins-1)
	n.beta.det=ncol(DM.det)
	n.Records=Meta$G.transect*Meta$n.Observers
	
	#initialize G.obs (number of groups observed per transect)
	G.obs=Meta$G.transect
	for(itrans in 1:Meta$n.transects){
		Tmp=matrix(Data[itrans,1:n.Records[itrans],2],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
		G.obs[itrans]=sum(apply(Tmp,1,'sum')>0)
	}
	
	#initialize Y.tilde (temporarily pretending correlation between observations is zero)
	Y.tilde=matrix(0,max(Meta$M),Meta$n.transects)
	for(itrans in 1:Meta$n.transects){
		X=get_mod_matrix(Cur.dat=Data[itrans,,],Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
		ExpY=X%*%Par$det
		Y.tilde[,itrans]=rtruncnorm(max(Meta$M), a=ifelse(Data[itrans,,2]==0,-Inf,0), b=ifelse(Data[itrans,,2]==0,0,Inf), ExpY, 1)		
	}
	
	#initialize MCMC, Acceptance rate matrices
	mcmc.length=(Control$iter-Control$burnin)/Control$thin
	MCMC=list(N.tot=rep(0,mcmc.length),N=matrix(0,mcmc.length,Meta$S),G=matrix(0,mcmc.length,Meta$S),Hab=data.frame(matrix(0,mcmc.length,length(Par$hab))),Det=data.frame(matrix(0,mcmc.length,length(Par$det))),cor=rep(0,mcmc.length),tau.eta=rep(0,mcmc.length),tau.nu=rep(0,mcmc.length))
	colnames(MCMC$Hab)=colnames(DM.hab)
	colnames(MCMC$Det)=colnames(DM.det)
	Accept=list(cor=0,N=rep(0,Meta$n.transects),Nu=0,Hab=rep(0,length(Par$hab)))
		
	st <- Sys.time()
	##################################################
	############   Begin MCMC Algorithm ##############
	##################################################
	for(iiter in 1:cur.iter){
		#if((iiter%%1000)==1)cat(paste('\n ', iiter))
		#cat(paste('\n ', iiter))
		
		########## update abundance parameters at the strata scale   ################
		
		#update nu parameters (log lambda)
		Mu=DM.hab%*%Par$hab+Par$Eta
		Grad1=sapply(Lam.index,'log_lambda_gradient',Mu=Mu,Nu=Par$Nu,N=Par$G,var.nu=1/Par$tau.nu)
		Prop=Par$Nu+Control$MH.nu^2*0.5*Grad1+Control$MH.nu*rnorm(Meta$S)
		new.post=log_lambda_log_likelihood(Log.lambda=Prop,DM=DM.hab,Beta=Par$hab,SD=sqrt(1/Par$tau.nu),N=Par$G)
		old.post=log_lambda_log_likelihood(Log.lambda=Par$Nu,DM=DM.hab,Beta=Par$hab,SD=sqrt(1/Par$tau.nu),N=Par$G)
		Grad2=sapply(Lam.index,'log_lambda_gradient',Mu=Mu,Nu=Prop,N=Par$G,var.nu=1/Par$tau.nu)
		diff1=as.vector(Par$Nu-Prop-0.5*Control$MH.nu^2*Grad2)	
		diff2=as.vector(Prop-Par$Nu-0.5*Control$MH.nu^2*Grad1)
		log.jump=0.5/Control$MH.nu^2*(sqrt(crossprod(diff1,diff1))-sqrt(crossprod(diff2,diff2))) #ratio of jumping distributions using e.g. Robert and Casella 2004 p. 319	
		if(runif(1)<exp(new.post-old.post+log.jump)){
			Par$Nu=Prop
			Accept$Nu=Accept$Nu+1	
		}
		
		if(Meta$spat.ind==0){
			#update eta parameters (spatial random effects)
			V.eta.inv <- Par$tau.nu*diag(Meta$S) + Par$tau.eta*Q
			M.eta <- solve(V.eta.inv, Par$tau.nu*(Par$Nu-DM.hab%*%Par$hab))		
			Par$Eta<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
			Par$Eta=Par$Eta-mean(Par$Eta)  #centering
			
			#update tau_eta  (precision of spatial process)
			Par$tau.eta <- rgamma(1, (Meta$S-1)/2 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta, Q %*% Par$Eta)/2) + Prior.pars$b.eta)
		}
		
		#update tau_nu	 (precision for Poisson overdispersion)
		Mu=DM.hab%*%Par$hab+Par$Eta
		Diff=Par$Nu-Mu
		Par$tau.nu <- rgamma(1,Meta$S/2 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))/2 + Prior.pars$b.nu)
		
		
		#translate to lambda scale
		Lambda=Meta$Area.hab*exp(Par$Nu)
		Lambda.trans=Lambda[Meta$Mapping]
		
		#update Betas for habitat relationships
		for(ipar in 1:length(Par$hab)){
			Prop=Par$hab
			Prop[ipar]=Par$hab[ipar]+Control$MH.beta[ipar]*rnorm(1,0,1)
			old.post=sum(dnorm(Par$Nu,Mu,sqrt(1/Par$tau.nu),log=1))+sum(dnorm(Par$hab,0,Prior.pars$beta.sd,log=1))
			Mu=DM.hab%*%Prop+Par$Eta			
			new.post=sum(dnorm(Par$Nu,Mu,sqrt(1/Par$tau.nu),log=1))+sum(dnorm(Prop,0,Prior.pars$beta.sd,log=1))
			if(runif(1)<exp(new.post-old.post)){
				Par$hab[ipar]=Prop[ipar]
				Accept$Hab[ipar]=Accept$Hab[ipar]+1
			}
		}
		
		########## update group abundance at strata level
		rsamp=switch_sample(n=(Meta$M[itrans]-n.Records[itrans])/Meta$n.Observers[itrans],pdf=Meta$Cov.prior.pdf[1],cur.par=Meta$Cov.prior.parms[,1])
		
		Par$G=rpois(Meta$S,Lambda*(1-Meta$Covered.area))
		Par$N=Par$G+rpois(Meta$S,Meta$Cov.prior.parms[1,1]*Par$G) #add the Par$G since number in group > 0 
		Par$G[Meta$Mapping]=Par$G[Meta$Mapping]+Meta$G.transect
		Par$N[Meta$Mapping]=Par$N[Meta$Mapping]+Meta$N.transect
		
		########## update abundance, distances, ind. covariates for observed transects using RJMCMC  #############
		for(itrans in 1:Meta$n.transects){
			Sample=c(-Control$RJ.N[itrans]:Control$RJ.N[itrans])
			Sample=Sample[-(Control$RJ.N[itrans]+1)] #don't make proposals where pop size stays the same
			a=sample(Sample,1)
			Sigma=matrix(Par$cor,Meta$n.Observers[itrans],Meta$n.Observers[itrans])
			diag(Sigma)=1
			offdiag=which(Sigma!=1)
			
			if(a>0){ # proposal an addition
				if(((Meta$G.transect[itrans]+a)*Meta$n.Observers[itrans])>Meta$M[itrans])cat('\n Error: proposed abundance > M; increase M value! \n')
				else{
					Cur.dat=Data[itrans,(n.Records[itrans]+1):(n.Records[itrans]+a*Meta$n.Observers[itrans]),]
					X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY=X%*%Par$det
					P=c(1:a)
					for(i in 1:a){
						Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
						P[i]=pmvnorm(upper=c(0,0),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:a){
						tmp.sum=tmp.sum-log(Meta$G.transect[itrans]-G.obs[itrans]+i)+log(P[i])
					}
					MH.prob=exp(a*log(Lambda.trans[itrans]*Meta$Area.trans[itrans])+tmp.sum)
					if(runif(1)<MH.prob){
						Meta$G.transect[itrans]=Meta$G.transect[itrans]+a
						n.Records[itrans]=Meta$G.transect[itrans]*Meta$n.Observers[itrans]
						if(Meta$grps==FALSE)Meta$N.transect[itrans]=Meta$G.transect[itrans]
						else Meta$N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(Meta$stacked.names=="Group")])/Meta$n.Observers[itrans]
						Accept$N[itrans]=Accept$N[itrans]+1
						#generate Y-tilde values
						Tmp=matrix(ExpY,a,Meta$n.Observers[itrans],byrow=TRUE)
						Y.tilde.temp=rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,1],sd=1)
						if(Meta$n.Observers[itrans]>1){
							Dist=matrix(Cur.dat[,Meta$dist.pl],a,Meta$n.Observers[itrans],byrow=TRUE)
							Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
							Y.tilde.temp=cbind(Y.tilde.temp,rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,2]+Cor*(Y.tilde.temp-Tmp[,2]),sd=sqrt(1-Cor^2)))
						}
						Y.tilde[(n.Records[itrans]-a*Meta$n.Observers[itrans]+1):n.Records[itrans],itrans]=as.vector(t(Y.tilde.temp))
					}	
				}
			}
			else{ #proposal a deletion
				if((Meta$G.transect[itrans]+a)>=G.obs[itrans]){ #can only delete records where an animal isn't observed
					Cur.dat=Data[itrans,n.Records[itrans]:(n.Records[itrans]+a*Meta$n.Observers[itrans]+1),]
					X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
					ExpY=X%*%Par$det
					P=c(1:-a)
					for(i in 1:-a){
						Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
						P[i]=pmvnorm(upper=rep(0,Meta$n.Observers[itrans]),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:-a){
						tmp.sum=tmp.sum+log(Meta$G.transect[itrans]-G.obs[itrans]-i+1)-log(P[i])
					}
					MH.prob=exp(a*log(Lambda.trans[itrans]*Meta$Area.trans[itrans])+tmp.sum)
					if(runif(1)<MH.prob){
						Meta$G.transect[itrans]=Meta$G.transect[itrans]+a
						n.Records[itrans]=Meta$G.transect[itrans]*Meta$n.Observers[itrans]
						Accept$N[itrans]=Accept$N[itrans]+1
					}					
					
				}
			}
			
			#update distances, individual covariates for latent animals not currently in the population
			#fill distances for unobserved
			if(Meta$i.binned==1)Data[itrans,(n.Records[itrans]+1):Meta$M[itrans],Meta$dist.pl]=rep(sample(c(1:Meta$n.bins),size=(Meta$M[itrans]-n.Records[itrans])/Meta$n.Observers[itrans],replace=TRUE,prob=Meta$Bin.length),each=Meta$n.Observers[itrans])
			else Data[itrans,(n.Records[itrans]+1):Meta$M[itrans],Meta$dist.pl]=rep(runif((Meta$M[itrans]-n.Records[itrans])/Meta$n.Observers[itrans]),each=Meta$n.Observers[itrans])
			#fill individual covariate values for (potential) animals that weren't observed
			if(Meta$n.ind.cov>0){
				for(icov in 1:Meta$n.ind.cov){
					rsamp=switch_sample(n=(Meta$M[itrans]-n.Records[itrans])/Meta$n.Observers[itrans],pdf=Meta$Cov.prior.pdf[icov],cur.par=Meta$Cov.prior.parms[,icov])
					Data[itrans,(n.Records[itrans]+1):Meta$M[itrans],Meta$dist.pl+icov]=rep(rsamp,each=Meta$n.Observers[itrans])
				}
			}
			
			#update distances, individual covariates for animals that ARE in the population but never observed
			cur.G=Meta$G.transect[itrans]-G.obs[itrans]
			if(Meta$G.transect[itrans]>G.obs[itrans]){
				#distance
				if(Meta$i.binned==1)dist.star=sample(c(1:Meta$n.bins),cur.G,replace=TRUE,prob=Meta$Bin.length)
				else dist.star=runif(cur.G)
				Cur.dat=Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),]
				cur.dist=Cur.dat[,Meta$dist.pl]
				X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
				ExpY=X%*%Par$det
				L.old=c(1:length(dist.star))
				Tmp.Y.tilde=Y.tilde[(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),itrans]
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
				Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),Meta$dist.pl]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),Meta$dist.pl]+rep(Acc,each=Meta$n.Observers[itrans])*rep(dist.star,each=Meta$n.Observers[itrans])
				
				#individual covariates
				if(Meta$n.ind.cov>0){
					for(icov in 1:Meta$n.ind.cov){
						Cov.star=switch_sample(n=cur.G,pdf=Meta$Cov.prior.pdf[icov],cur.par=Meta$Cov.prior.parms[,icov])
						Cur.dat=Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),]
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
						Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]+rep(Acc,each=Meta$n.Observers[itrans])*rep(Cov.star,each=Meta$n.Observers[itrans])						
					}
				}
			}
			if(Meta$grps==FALSE)Meta$N.transect[itrans]=Meta$G.transect[itrans]
			else Meta$N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(Meta$stacked.names=="Group")])/Meta$n.Observers[itrans]
			
			
			#update Y.tilde
			Cur.dat=Data[itrans,1:n.Records[itrans],]
			X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
			ExpY=matrix(X%*%Par$det,Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
			Temp.Y.tilde=matrix(Y.tilde[1:n.Records[itrans],itrans],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
			Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
			Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
			Resp=matrix(Cur.dat[,2],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
			EY1=ExpY[,1]+Cor*(Temp.Y.tilde[,2]-ExpY[,2])
			Temp.Y.tilde[,1] <- rtruncnorm(Meta$G.transect[itrans], a=ifelse(Resp[,1]==0,-Inf,0), b=ifelse(Resp[,1]==0,0,Inf), EY1, sqrt(1-Cor^2))
			EY2=ExpY[,2]+Cor*(Temp.Y.tilde[,1]-ExpY[,1])
			Temp.Y.tilde[,2] <- rtruncnorm(Meta$G.transect[itrans], a=ifelse(Resp[,2]==0,-Inf,0), b=ifelse(Resp[,2]==0,0,Inf), EY2, sqrt(1-Cor^2))
			Y.tilde[1:n.Records[itrans],itrans]=as.vector(t(Temp.Y.tilde))
		}
		
		###############       update detection process parameters       ##############
		# First, assemble stacked adjusted Response, X matrices across all transects; 
		#basic form of response is Ytilde[obs1]-cor*Ytilde[obs2]
		#adjusted DM rows are of form X[obs1]-cor*X[obs2]
		Cur.dat=Data[1,1:n.Records[1],]
		Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[1],Meta$n.Observers[1],byrow=TRUE)[,1]
		X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[1]))
		if(Meta$n.Observers[1]==2){
			Tmp.cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
			Cor=rep(Tmp.cor,2)  #assemble vector of correlation parameters for each observation
			X.beta=rbind(t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
			Y.temp=matrix(Y.tilde[1:n.Records[1],1],Meta$G.transect[1],2,byrow=TRUE)
			Y.beta=c(Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
		}
		else{
			X.beta=t(X.temp[,1,])
			Y.beta=Y.tilde[1:n.Records[1],1]
			Cor=rep(1,n.Records[1])
		}
		for(itrans in 2:Meta$n.transects){
			Cur.dat=Data[itrans,1:n.Records[itrans],]
			Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)[,1]
			X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[itrans]))
			if(Meta$n.Observers[1]==2){
				Tmp.cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
				Cor=c(Cor,rep(Tmp.cor,2))  #assemble vector of correlation parameters for each observation
				X.beta=rbind(X.beta,t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
				Y.temp=matrix(Y.tilde[1:n.Records[itrans],itrans],Meta$G.transect[itrans],2,byrow=TRUE)
				Y.beta=c(Y.beta,Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
			}
			else{
				X.beta=rbind(X.beta,t(X.temp[,1,]))
				Y.beta=c(Y.beta,Y.tilde[1:n.Records[itrans],1])
				Cor=c(Cor,rep(1,n.Records[itrans]))
			}
		}
		#now use basic matrix equations from Gelman '04 book (14.11 14.12) to update beta parms
		Sig.inv=sqrt(1-Cor^2)  #for use in eq. 14.11, 14.12 of Gelman et al.; don't need a matrix since it would be diagonal
		V.inv <- crossprod(X.beta,Sig.inv*X.beta) 	
		M.z <- solve(V.inv, crossprod(X.beta,Sig.inv*Y.beta))
		Par$det <- M.z + solve(chol(V.inv), rnorm(n.beta.det,0,1))
		
		
		#update correlation parameter for detection process (if applicable)
		if(Meta$point.ind==1){
			cor.star=Par$cor+runif(1,-Control$MH.cor,Control$MH.cor)
			if(cor.star>0 & cor.star<1){
				I.gt.one=which(Meta$n.Observers>1)
				n.gt.one=length(I.gt.one)
				
				Cur.dat=Data[I.gt.one[1],1:n.Records[I.gt.one[1]],]
				Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[I.gt.one[1]],Meta$n.Observers[I.gt.one[1]],byrow=TRUE)[,1]
				X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[I.gt.one[1]]))
				Y.temp=matrix(Y.tilde[1:n.Records[I.gt.one[1]],I.gt.one[1]],Meta$G.transect[I.gt.one[1]],2,byrow=TRUE)			
				Delta1=Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det)
				Delta2=Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det)
				
				if(n.gt.one>1){
					for(itrans in 2:n.gt.one){
						Cur.dat=Data[I.gt.one[itrans],1:n.Records[I.gt.one[itrans]],]
						Dist=c(Dist,matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[I.gt.one[itrans]],Meta$n.Observers[I.gt.one[itrans]],byrow=TRUE)[,1])
						X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[I.gt.one[itrans]]))
						Y.temp=matrix(Y.tilde[1:n.Records[I.gt.one[itrans]],I.gt.one[itrans]],Meta$G.transect[I.gt.one[itrans]],2,byrow=TRUE)			
						Delta1=c(Delta1,Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det))
						Delta2=c(Delta2,Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det))
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
		
		if(adapt==TRUE){
			if(iiter%%100==0){
				if(Accept$cor<30)Control$MH.cor=Control$MH.cor*.95
				if(Accept$cor>40)Control$MH.cor=Control$MH.cor*1.053
				for(ipar in 1:length(Control$MH.beta)){
					if(Accept$Hab[ipar]<30)Control$MH.beta[ipar]=Control$MH.beta[ipar]*.95
					if(Accept$Hab[ipar]>40)Control$MH.beta[ipar]=Control$MH.beta[ipar]*1.053
				}
				if(Accept$Nu<55)Control$MH.nu=Control$MH.nu*.95
				if(Accept$Nu>60)Control$MH.nu=Control$MH.nu*1.053
				for(ipar in 1:Meta$n.transects){
					if(Accept$N[ipar]<30)Control$RJ.N[ipar]=max(1,Control$RJ.N[ipar]-1)
					if(Accept$N[ipar]>40)Control$RJ.N[ipar]=Control$RJ.N[ipar]+1
				}
				Accept$cor=0
				Accept$Hab=Accept$Hab*0
				Accept$Nu=0
				Accept$N=Accept$N*0
			}
		}
		
		#store results if applicable
		if(iiter>Control$burnin & iiter%%Control$thin==0){
			MCMC$G[(iiter-Control$burnin)/Control$thin,]=Par$G
			MCMC$N[(iiter-Control$burnin)/Control$thin,]=Par$N
			MCMC$N.tot[(iiter-Control$burnin)/Control$thin]=sum(Par$N)
			MCMC$cor[(iiter-Control$burnin)/Control$thin]=Par$cor
			MCMC$Hab[(iiter-Control$burnin)/Control$thin,]=Par$hab
			MCMC$Det[(iiter-Control$burnin)/Control$thin,]=Par$det
			MCMC$tau.eta[(iiter-Control$burnin)/Control$thin]=Par$tau.eta
			MCMC$tau.nu[(iiter-Control$burnin)/Control$thin]=Par$tau.nu
			
		}
		
		if(iiter==15){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/15
			ttc <- round((cur.iter-15)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}
		
	}
	Out=list(MCMC=MCMC,Accept=Accept,Control=Control)
	Out
}

