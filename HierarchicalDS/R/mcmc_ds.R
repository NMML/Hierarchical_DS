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
#' 	"Adj"			Adjacency matrix giving connectivity of spatial grid cells
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
#'  "a.linex": Optimal parameter for linex loss function
#' @export
#' @import Matrix
#' @keywords areal, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn


mcmc_ds<-function(Par,Data,cur.iter,adapt,Control,DM.hab,Q,Prior.pars,Meta){	
	#require(mvtnorm)
	#require(Matrix)
	#require(truncnorm)
	
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
	Pred.G=matrix(0,mcmc.length,Meta$n.transects)
	Obs.G=Pred.G
	
	colnames(MCMC$Hab)=colnames(DM.hab)
	Accept=list(cor=0,N=rep(0,Meta$n.transects),Nu=0,Hab=rep(0,length(Par$hab)))
	
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
		#2) ssimulate nu for areas not sampled
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
		
		########## update group abundance at strata level
		Par$G=rpois(Meta$S,Lambda*(1-Meta$Covered.area))
		Par$G[Meta$Mapping]=Par$G[Meta$Mapping]+Meta$G.transect

		######update abundance at transect level
		#(assumed known with certainty)
		
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
		if(adapt==TRUE){
			if(iiter%%100==0){
				if(Accept$Nu<55)Control$MH.nu=Control$MH.nu*.95
				if(Accept$Nu>60)Control$MH.nu=Control$MH.nu*1.053
				Accept$cor=0
				Accept$Hab=Accept$Hab*0
				Accept$Nu=0
				Accept$N=Accept$N*0
			}
		}
		
		#store results if applicable
		if(iiter>Control$burnin & iiter%%Control$thin==0){
			MCMC$G[(iiter-Control$burnin)/Control$thin,]=Par$G
			MCMC$N.tot[(iiter-Control$burnin)/Control$thin]=sum(Par$G)
			MCMC$Hab[(iiter-Control$burnin)/Control$thin,]=Par$hab
			MCMC$tau.eta[(iiter-Control$burnin)/Control$thin]=Par$tau.eta
			MCMC$tau.nu[(iiter-Control$burnin)/Control$thin]=Par$tau.nu
			
			Obs.G[(iiter-Control$burnin)/Control$thin,]=Meta$G.transect
			Pred.G[(iiter-Control$burnin)/Control$thin,]=Meta$Area.hab[Mapping]*Meta$Area.trans*exp(rnorm(Meta$n.transects,(DM.hab%*%Par$hab+Par$Eta)[Mapping],sqrt(1/Par$tau.nu)))
			
		}
		
		if(iiter==100){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/100
			ttc <- round((cur.iter-100)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}	
		if((iiter%%1000)==1)cat(paste('iteration ', iiter,' of ',cur.iter,' completed \n'))
	}
	a.linex=calc_linex_a(Pred.G,Obs.G)$minimum
	cat(paste('\n total elapsed time: ',difftime(Sys.time(),st,units="mins"),' minutes \n'))
	Out=list(MCMC=MCMC,Accept=Accept,Control=Control,Obs.G=Obs.G,Pred.G=Pred.G)
	Out
}

