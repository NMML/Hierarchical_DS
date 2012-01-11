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
#'  "point.ind"		Indicator for whether point independence assumed (if no, then no correlation modeled b/w multiple observers as function of distance)
#'  "fix.tau.nu"	Indicator for whether tau.nu should be fixed (1) or estimated(0)
#'  "srr"			Indicator for whether a spatially restricted regression model should be employed (1) or not (0)
#'  "srr.tol"		Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation
#' @return returns a list with the following slots: 
#' 	"MCMC": A list object containing posterior samples;
#'  "Accept": A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms;
#'  "Control": A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used) 
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

	if(Meta$srr==TRUE){
		P.c=diag(Meta$S)-DM.hab%*%solve(crossprod(DM.hab,DM.hab),t(DM.hab))
		Omega=(P.c%*%Meta$Adj%*%P.c)*(Meta$S/sum(Adj))
		Eigen=eigen(Omega)
		if(max(Eigen$values)<srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < srr.tol; decrease srr.tol"))
		Ind=which(Eigen$values>Meta$srr.tol)
		L.t=Eigen$vectors[,Ind]
		cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
		L=t(L.t)
		Qt=L%*%Q%*%L.t
		cross.L=L%*%L.t	
		n.theta=nrow(Qt)
		Theta=rnorm(n.theta,0,sqrt(1/Par$tau.eta))
	}
	
	Sampled=unique(Meta$Mapping)
	n.unique=length(Sampled)
	Sampled.area.by.strata=rep(0,n.unique)
	for(i in 1:Meta$n.transects)Sampled.area.by.strata[which(Sampled==Meta$Mapping[i])]=Sampled.area.by.strata[which(Sampled==Meta$Mapping[i])]+Meta$Area.trans[i]
	
	#initialize MCMC, Acceptance rate matrices
	mcmc.length=(Control$iter-Control$burnin)/Control$thin
	MCMC=list(N.tot=rep(0,mcmc.length),N=matrix(0,mcmc.length,Meta$S),G=matrix(0,mcmc.length,Meta$S),Hab=data.frame(matrix(0,mcmc.length,length(Par$hab))),Det=data.frame(matrix(0,mcmc.length,length(Par$det))),cor=rep(0,mcmc.length),tau.eta=rep(0,mcmc.length),tau.nu=rep(0,mcmc.length),Cov.par=matrix(0,mcmc.length,length(Par$Cov.par)))
	colnames(MCMC$Hab)=colnames(DM.hab)
	colnames(MCMC$Det)=colnames(DM.det)
	Accept=list(cor=0,N=rep(0,Meta$n.transects),Nu=0,Hab=rep(0,length(Par$hab)))
	Pred.N=matrix(0,mcmc.length,Meta$n.transects)
	Obs.N=Pred.N
	
	#initialize random effect matrices for individual covariates if required
	if(sum(1-Meta$Cov.prior.fixed)>0)RE.cov=array(0,dim=c(Meta$n.transects,max(Meta$M),Meta$n.ind.cov))
		
	st <- Sys.time()
	##################################################
	############   Begin MCMC Algorithm ##############
	##################################################
	for(iiter in 1:cur.iter){
		#cat(paste('\n ', iiter))
		
		########## update abundance parameters at the strata scale   ################
		
		#update nu parameters (log lambda)
		#1) for sampled cells
		Mu=DM.hab%*%Par$hab+Par$Eta
		#cat(paste("\n Mu",Mu))
		#cat(paste("\n Nu",Par$Nu))
		G.sampled=rep(0,n.unique) #total number of groups currently in each sampled strata
		for(i in 1:Meta$n.transects)G.sampled[which(Sampled==Meta$Mapping[i])]=G.sampled[which(Sampled==Meta$Mapping[i])]+Meta$G.transect[i]
		Grad1=log_lambda_gradient(Mu=Mu,Nu=Par$Nu,Sampled=Sampled,Area=Sampled.area.by.strata,N=G.sampled,var.nu=1/Par$tau.nu)
		Prop=Par$Nu
		Prop[Sampled]=Par$Nu[Sampled]+Control$MH.nu^2*0.5*Grad1+Control$MH.nu*rnorm(n.unique)
		new.post=log_lambda_log_likelihood(Log.lambda=Prop[Sampled],DM=DM.hab,Beta=Par$hab,SD=sqrt(1/Par$tau.nu),N=G.sampled,Sampled=Sampled,Area=Sampled.area.by.strata)
		old.post=log_lambda_log_likelihood(Log.lambda=Par$Nu[Sampled],DM=DM.hab,Beta=Par$hab,SD=sqrt(1/Par$tau.nu),N=G.sampled,Sampled=Sampled,Area=Sampled.area.by.strata)
		Grad2=log_lambda_gradient(Mu=Mu,Nu=Prop,Sampled=Sampled,Area=Sampled.area.by.strata,N=G.sampled,var.nu=1/Par$tau.nu)
		diff1=as.vector(Par$Nu[Sampled]-Prop[Sampled]-0.5*Control$MH.nu^2*Grad2)	
		diff2=as.vector(Prop[Sampled]-Par$Nu[Sampled]-0.5*Control$MH.nu^2*Grad1)
		log.jump=0.5/Control$MH.nu^2*(sqrt(crossprod(diff1,diff1))-sqrt(crossprod(diff2,diff2))) #ratio of jumping distributions using e.g. Robert and Casella 2004 p. 319	
		#cat(paste("\n iter=",iiter," new=",new.post," old=",old.post," jump=",log.jump))
		if(runif(1)<exp(new.post-old.post+log.jump)){
			Par$Nu=Prop
			Accept$Nu=Accept$Nu+1	
		}
		#2) simulate nu for areas not sampled
		Par$Nu[-Sampled]=rnorm(Meta$S-n.unique,Mu[-Sampled],1/sqrt(Par$tau.nu))
		
		
		if(Meta$spat.ind==0){
			if(Meta$srr==FALSE){
				#update eta parameters (spatial random effects)
				V.eta.inv <- Par$tau.nu*diag(Meta$S) + Par$tau.eta*Q
				M.eta <- solve(V.eta.inv, Par$tau.nu*(Par$Nu-DM.hab%*%Par$hab))		
				Par$Eta<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(Meta$S,0,1)))
				Par$Eta=Par$Eta-mean(Par$Eta)  #centering
				
				#update tau_eta  (precision of spatial process)
				Par$tau.eta <- rgamma(1, (Meta$S-1)/2 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta, Q %*% Par$Eta)/2) + Prior.pars$b.eta)
			}
			else{
				#Update Theta
				Dat.minus.Exp=Par$Nu-DM.hab%*%Par$hab
				V.eta.inv <- cross.L + Par$tau.eta*Qt
				M.eta <- solve(V.eta.inv, L%*%Dat.minus.Exp)
				Theta <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(n.theta,0,1))
				Par$Eta=as.numeric(L.t%*%Theta)
				
				#update tau.eta
				Par$tau.eta <- rgamma(1, n.theta/2 + Prior.pars$a.eta, as.numeric(crossprod(Theta, Qt %*% Theta)/2) + Prior.pars$b.eta)
			}
		}
		
		#update tau_nu	 (precision for Poisson overdispersion)
		Mu=DM.hab%*%Par$hab+Par$Eta
		if(Meta$fix.tau.nu==FALSE){
			Diff=Par$Nu[Sampled]-Mu[Sampled]
			Par$tau.nu <- rgamma(1,n.unique/2 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))/2 + Prior.pars$b.nu)
		}
		
		#translate to lambda scale
		Lambda=Meta$Area.hab*exp(Par$Nu)
		Lambda.trans=Lambda[Meta$Mapping]*Meta$Area.trans
		
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
		Par$G=rpois(Meta$S,Lambda*(1-Meta$Covered.area))
		grp.lam=ifelse(Meta$Cov.prior.pdf[1] %in% c("pois1_ln","poisson_ln"),exp(Par$Cov.par[1,1]+(Par$Cov.par[2,1])^2/2),Par$Cov.par[1,1])
		Par$N=Par$G+rpois(Meta$S,grp.lam*Par$G) #add the Par$G since number in group > 0 
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
				if(((Meta$G.transect[itrans]+a)*Meta$n.Observers[itrans])>=Meta$M[itrans])cat(paste('\n Warning: proposed abundance for transect ',itrans,' > M; consider increasing M value! \n'))
				else{
					Cur.dat=Data[itrans,(n.Records[itrans]+1):(n.Records[itrans]+a*Meta$n.Observers[itrans]),]
					if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
					X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names=Meta$stacked.names,factor.ind=Meta$factor.ind,Det.formula=Meta$Det.formula,Levels=Meta$Levels)
					ExpY=X%*%Par$det
					P=c(1:a)
					for(i in 1:a){
						Sigma[offdiag]=Par$cor*(Meta$i.binned*(Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl]-1)*dist.mult+(1-Meta$i.binned)*Cur.dat[i*Meta$n.Observers[itrans],Meta$dist.pl])
						P[i]=pmvnorm(upper=rep(0,Meta$n.Observers[itrans]),mean=ExpY[(i*Meta$n.Observers[itrans]-Meta$n.Observers[itrans]+1):(i*Meta$n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:a){
						tmp.sum=tmp.sum-log(Meta$G.transect[itrans]-G.obs[itrans]+i)+log(P[i])
					}
					MH.prob=exp(a*log(Lambda.trans[itrans])+tmp.sum)
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
						tmp.sum=tmp.sum+log(Meta$G.transect[itrans]-G.obs[itrans]-i+1)-log(P[i])
					}
					MH.prob=exp(a*log(Lambda.trans[itrans])+tmp.sum)
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
					if(Meta$Cov.prior.pdf[icov]=='poisson_ln' | Meta$Cov.prior.pdf[icov]=='pois1_ln')cur.RE=RE.cov[itrans,(Meta$G.transect[itrans]+1):(Meta$M[itrans]/Meta$n.Observers[itrans]),icov]
					else cur.RE=0
					rsamp=switch_sample(n=(Meta$M[itrans]-n.Records[itrans])/Meta$n.Observers[itrans],pdf=Meta$Cov.prior.pdf[icov],cur.par=Par$Cov.par[,icov],RE=cur.RE)
					Data[itrans,(n.Records[itrans]+1):Meta$M[itrans],Meta$dist.pl+icov]=rep(rsamp,each=Meta$n.Observers[itrans])
				}
			}
			
			#update distances, individual covariates for animals that ARE in the population but never observed
			cur.G=Meta$G.transect[itrans]-G.obs[itrans]
			if(cur.G>0){
				#distance
				if(Meta$i.binned==1)dist.star=sample(c(1:Meta$n.bins),cur.G,replace=TRUE,prob=Meta$Bin.length)
				else dist.star=runif(cur.G)
				Cur.dat=Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),]
				if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
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
						if(Meta$Cov.prior.pdf[icov]=='poisson_ln' | Meta$Cov.prior.pdf[icov]=='pois1_ln')cur.RE=RE.cov[itrans,(G.obs[itrans]+1):Meta$G.transect[itrans],icov]
						else cur.RE=0
						Cov.star=switch_sample(n=cur.G,pdf=Meta$Cov.prior.pdf[icov],cur.par=Par$Cov.par[,icov],RE=cur.RE)
						Cur.dat=Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),]
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
						Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]=(1-rep(Acc,each=Meta$n.Observers[itrans]))*Data[itrans,(G.obs[itrans]*Meta$n.Observers[itrans]+1):(Meta$G.transect[itrans]*Meta$n.Observers[itrans]),Meta$dist.pl+icov]+rep(Acc,each=Meta$n.Observers[itrans])*rep(Cov.star,each=Meta$n.Observers[itrans])						
					}
				}
			}
			if(Meta$grps==FALSE)Meta$N.transect[itrans]=Meta$G.transect[itrans]
			else Meta$N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(Meta$stacked.names=="Group")])/Meta$n.Observers[itrans]
			
			
			#update Y.tilde
			if(n.Records[itrans]>0){
				Cur.dat=Data[itrans,1:n.Records[itrans],]
				if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
				X=get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)
				ExpY=matrix(X%*%Par$det,Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
				Temp.Y.tilde=matrix(Y.tilde[1:n.Records[itrans],itrans],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
				if(Meta$n.Observers[itrans]==1){
					Temp.Y.tilde<-rtruncnorm(Meta$G.transect[itrans],a=ifelse(Cur.dat[,2]==0,-Inf,0),b=ifelse(Cur.dat[,2]==0,0,Inf),ExpY,1)
				}
				else{
					Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
					Cor=Par$cor*(Meta$i.binned*(Dist[,1]-1)*dist.mult+(1-Meta$i.binned)*Dist[,1])
					Resp=matrix(Cur.dat[,2],Meta$G.transect[itrans],Meta$n.Observers[itrans],byrow=TRUE)
					EY1=ExpY[,1]+Cor*(Temp.Y.tilde[,2]-ExpY[,2])
					Temp.Y.tilde[,1] <- rtruncnorm(Meta$G.transect[itrans], a=ifelse(Resp[,1]==0,-Inf,0), b=ifelse(Resp[,1]==0,0,Inf), EY1, sqrt(1-Cor^2))
					EY2=ExpY[,2]+Cor*(Temp.Y.tilde[,1]-ExpY[,1])
					Temp.Y.tilde[,2] <- rtruncnorm(Meta$G.transect[itrans], a=ifelse(Resp[,2]==0,-Inf,0), b=ifelse(Resp[,2]==0,0,Inf), EY2, sqrt(1-Cor^2))
				}
				Y.tilde[1:n.Records[itrans],itrans]=as.vector(t(Temp.Y.tilde))
			}
		}
		
		###############       update detection process parameters       ##############
		# First, assemble stacked adjusted Response, X matrices across all transects; 
		#basic form of response is Ytilde[obs1]-cor*Ytilde[obs2]
		#adjusted DM rows are of form X[obs1]-cor*X[obs2]
		GT0=which(Meta$G.transect>0)
		n.gt0=length(GT0)
		Cur.dat=Data[GT0[1],1:n.Records[GT0[1]],]
		Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
		X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,Meta$n.Observers[GT0[1]],Meta$G.transect[GT0[1]]))
		if(Meta$n.Observers[GT0[1]]==2){
			Tmp.cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
			Cor=rep(Tmp.cor,2)  #assemble vector of correlation parameters for each observation
			X.beta=rbind(t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
			Y.temp=matrix(Y.tilde[1:n.Records[GT0[1]],1],Meta$G.transect[GT0[1]],2,byrow=TRUE)
			Y.beta=c(Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
		}
		else{
			X.beta=t(X.temp[,1,])
			Y.beta=Y.tilde[1:n.Records[GT0[1]],1]
			Cor=rep(0,n.Records[GT0[1]])
		}
		if(n.gt0>1){
			for(itrans in 2:n.gt0){
				Cur.dat=Data[itrans,1:n.Records[GT0[itrans]],]
				if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
				Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1]
				X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,Meta$n.Observers[GT0[itrans]],Meta$G.transect[GT0[itrans]]))
				if(Meta$n.Observers[GT0[itrans]]==2){
					Tmp.cor=Par$cor*(Meta$i.binned*(Dist-1)*dist.mult+(1-Meta$i.binned)*Dist)
					Cor=c(Cor,rep(Tmp.cor,2))  #assemble vector of correlation parameters for each observation
					X.beta=rbind(X.beta,t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
					Y.temp=matrix(Y.tilde[1:n.Records[GT0[itrans]],GT0[itrans]],Meta$G.transect[GT0[itrans]],2,byrow=TRUE)
					Y.beta=c(Y.beta,Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
				}
				else{
					X.beta=rbind(X.beta,t(X.temp[,1,]))
					Y.beta=c(Y.beta,Y.tilde[1:n.Records[GT0[itrans]],GT0[itrans]])
					Cor=c(Cor,rep(0,n.Records[GT0[itrans]]))
				}
			}
		}
		#now use basic matrix equations from Gelman '04 book (14.11 14.12) to update beta parms
		Sig.inv=1/(1-Cor^2)  #for use in eq. 14.11, 14.12 of Gelman et al.; don't need a matrix since it would be diagonal
		V.inv <- crossprod(X.beta,Sig.inv*X.beta) 	
		M.z <- solve(V.inv, crossprod(X.beta,Sig.inv*Y.beta))
		Par$det <- M.z + solve(chol(V.inv), rnorm(n.beta.det,0,1))
		
		
		#update correlation parameter for detection process (if applicable)
		if(Meta$point.ind==1){
			cor.star=Par$cor+runif(1,-Control$MH.cor,Control$MH.cor)
			if(cor.star>0 & cor.star<1){
				I.gt.one=which(Meta$n.Observers>1 & Meta$G.transect>0)
				n.gt.one=length(I.gt.one)
				
				Cur.dat=Data[I.gt.one[1],1:n.Records[I.gt.one[1]],]
				if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
				Dist=matrix(Cur.dat[,Meta$dist.pl],Meta$G.transect[I.gt.one[1]],2,byrow=TRUE)[,1]
				X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,Meta$stacked.names,Meta$factor.ind,Meta$Det.formula,Meta$Levels)),dim=c(n.beta.det,2,Meta$G.transect[I.gt.one[1]]))
				Y.temp=matrix(Y.tilde[1:n.Records[I.gt.one[1]],I.gt.one[1]],Meta$G.transect[I.gt.one[1]],2,byrow=TRUE)			
				Delta1=Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det)
				Delta2=Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det)
				
				if(n.gt.one>1){
					for(itrans in 2:n.gt.one){
						Cur.dat=Data[I.gt.one[itrans],1:n.Records[I.gt.one[itrans]],]
						if(is.vector(Cur.dat))Cur.dat=matrix(Cur.dat,1,length(Cur.dat))
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
		
		#update parameters of individual covariate distributions (if fixed=0)
		for(icov in 1:Meta$n.ind.cov){
			if(Meta$Cov.prior.fixed[icov]==0){
				if(Meta$Cov.prior.pdf[icov]=="normal")cat("\n Warning: hyper-priors not yet implemented for normal dist. \n")
				if(Meta$Cov.prior.pdf[icov]=="poisson"){
					Cur.cov=matrix(Data[GT0[1],1:n.Records[GT0[1]],Meta$dist.pl+icov],Meta$G.transect[GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
					if(n.gt0>1){
						for(itrans in 2:n.gt0){
							Cur.cov=c(Cur.cov,matrix(Data[GT0[itrans],1:n.Records[GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
						}
					}
					Par$Cov.par[1,icov]=rgamma(1,sum(Cur.cov)+Meta$Cov.prior.parms[1,icov],length(Cur.cov)+Meta$Cov.prior.parms[2,icov])
				}
				if(Meta$Cov.prior.pdf[icov]=="pois1"){
					Cur.cov=matrix(Data[GT0[1],1:n.Records[GT0[1]],Meta$dist.pl+icov],Meta$G.transect[GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
					if(n.gt0>1){
						for(itrans in 2:n.gt0){
							Cur.cov=c(Cur.cov,matrix(Data[GT0[itrans],1:n.Records[GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
						}
					}
					Par$Cov.par[1,icov]=rgamma(1,sum(Cur.cov)-length(Cur.cov)+Meta$Cov.prior.parms[1,icov],length(Cur.cov)+Meta$Cov.prior.parms[2,icov])
				}
				if(Meta$Cov.prior.pdf[icov]=="poisson_ln" | Meta$Cov.prior.pdf[icov]=="pois1_ln"){
					Cur.cov=matrix(Data[GT0[1],1:n.Records[GT0[1]],Meta$dist.pl+icov],Meta$G.transect[GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
					Cur.RE=RE.cov[GT0[1],1:Meta$G.transect[GT0[1]],icov]
					if(n.gt0>1){
						for(itrans in 2:n.gt0){
							Cur.cov=c(Cur.cov,matrix(Data[GT0[itrans],1:n.Records[GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
							Cur.RE=c(Cur.RE,RE.cov[GT0[itrans],1:Meta$G.transect[GT0[itrans]],icov])
						}
					}
					Cur.cov=Cur.cov-(Meta$Cov.prior.pdf[icov]=="pois1_ln")
					#1) update theta
					par.star=Par$Cov.par[1,icov]+runif(1,-0.05,0.05)
					sum.y=sum(Cur.cov)
					sum.yZ=sum(Cur.cov*Cur.RE)
					log.post.new=par.star*sum.y-sum(exp(par.star+Par$Cov.par[2,icov]*Cur.RE))+dnorm(par.star,Meta$Cov.prior.parms[1,icov],Meta$Cov.prior.parms[2,icov],log=1)
					log.post.old=Par$Cov.par[1,icov]*sum.y-sum(exp(Par$Cov.par[1,icov]+Par$Cov.par[2,icov]*Cur.RE))+dnorm(Par$Cov.par[1,icov],Meta$Cov.prior.parms[1,icov],Meta$Cov.prior.parms[2,icov],log=1)
				    if(runif(1)<exp(log.post.new-log.post.old))Par$Cov.par[1,icov]=par.star
					#2) update sigma
					par.star=Par$Cov.par[2,icov]+runif(1,-.01,.01)
					if(par.star>0 & par.star<Meta$Cov.prior.parms[3,icov]){
						log.post.new=par.star*sum.yZ-sum(exp(Par$Cov.par[1,icov]+par.star*Cur.RE))
						log.post.old=Par$Cov.par[2,icov]*sum.yZ-sum(exp(Par$Cov.par[1,icov]+Par$Cov.par[2,icov]*Cur.RE))
						if(runif(1)<exp(log.post.new-log.post.old))Par$Cov.par[2,icov]=par.star
					}
					#3) update random effects		
					for(itrans in 1:n.gt0){
						#animals currently in population
						Cur.cov=matrix(Data[GT0[itrans],1:n.Records[GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1]-(Meta$Cov.prior.pdf[icov]=="pois1_ln")						
						Cur.RE=RE.cov[GT0[itrans],1:Meta$G.transect[GT0[itrans]],icov]
						Prop=Cur.RE+runif(length(Cur.RE),-.1,.1)
						LogPost.new=Cur.cov*Par$Cov.par[2,icov]*Prop-exp(Par$Cov.par[1,icov]+Par$Cov.par[2,icov]*Prop)+dnorm(Prop,0,1,log=1)
						LogPost.old=Cur.cov*Par$Cov.par[2,icov]*Cur.RE-exp(Par$Cov.par[1,icov]+Par$Cov.par[2,icov]*Cur.RE)+dnorm(Cur.RE,0,1,log=1)
						Acc=(runif(length(Cur.RE))<(exp(LogPost.new-LogPost.old)))
						RE.cov[GT0[itrans],1:Meta$G.transect[GT0[itrans]],icov]=Acc*Prop+(1-Acc)*Cur.RE
					}
					#animals currently not in population
					for(itrans in 1:Meta$n.transects){
					    RE.cov[itrans,(Meta$G.transect[itrans]+1):Meta$M[itrans],icov]=rnorm(Meta$M[itrans]-Meta$G.transect[itrans],0,1)
					}
				}
				if(Meta$Cov.prior.pdf[icov]=="multinom"){
					Cur.cov=matrix(Data[1:GT0[1],1:n.Records[GT0[1]],Meta$dist.pl+icov],Meta$G.transect[GT0[1]],Meta$n.Observers[GT0[1]],byrow=TRUE)[,1]
					if(n.gt0>1){
						for(itrans in 2:n.gt0){
							Cur.cov=c(Cur.cov,matrix(Data[GT0[itrans],1:n.Records[GT0[itrans]],Meta$dist.pl+icov],Meta$G.transect[GT0[itrans]],Meta$n.Observers[GT0[itrans]],byrow=TRUE)[,1])
						}
					}
					Par$Cov.par[,icov]=rdirichlet(1,Meta$Cov.prior.parms[,icov]+tabulate(factor(Cur.cov)))
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
			MCMC$Cov.par[(iiter-Control$burnin)/Control$thin,]=Par$Cov.par
			Obs.N[(iiter-Control$burnin)/Control$thin,]=Meta$N.transect
			Temp.G=Meta$Area.hab[Mapping]*Meta$Area.trans*exp(rnorm(Meta$n.transects,(DM.hab%*%Par$hab+Par$Eta)[Mapping],sqrt(1/Par$tau.nu)))
			Pred.N[(iiter-Control$burnin)/Control$thin,]=Temp.G+rpois(Meta$n.transects,grp.lam*Temp.G)			
		}
		
		if(iiter==100){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/100
			ttc <- round((cur.iter-100)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}	
		if((iiter%%1000)==1)cat(paste('iteration ', iiter,' of ',cur.iter,' completed \n'))
	}
	cat(paste('\n total elapsed time: ',difftime(Sys.time(),st,units="mins"),' minutes \n'))
	Out=list(MCMC=MCMC,Accept=Accept,Control=Control,Obs.N=Obs.N,Pred.N=Pred.N)
	Out
}

