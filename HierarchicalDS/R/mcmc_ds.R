#' Function for MCMC analysis 
#' 
#' @param Par 	A list comprised of the following parameters:
#' 		\item{"det"}{a vector giving the current iteration's linear model parameters for the detection model}
#' 		\item{"hab"}{a vector giving the current iteration's linear model parameters for abundance intensity}
#' 		\item{"cor"}{a correlation parameter for detections that's an increasing function of distance (correlation at the maximum distance)}
#' 		\item{"Nu"}{a vector giving the log of the abundance intensity for each strata}
#' 		\item{"G"}{a vector giving the number of groups of animals in each strata } 
#' 		\item{"N"}{a vector giving the number of animals in each strata }
#' @param Data   A three dimensional array; the first dimension gives the transect, the second dimension indexes a (possible) observation, and the third dimension gives observations and covariates associated with a given animal.
#' 			Each of these columns is as for "Dat"
#' @param cur.iter   Number of iterations to run
#' @param adapt	If adapt==TRUE, run MCMC in adapt mode, optimizing MCMC proposal distributions prior to primary MCMC
#' @param Control	A list object including the following slots:
#'	\item{"iter"}{number of MCMC iterations}
#'  \item{"burnin"}{number of MCMC burnin iterations}
#'	\item{"thin"}{if specified, how many iterations to skip between recorded posterior samples}
#'	\item{"adapt"}{if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run} 
#'	\item{"MH.cor"}{Metropolis-hastings tuning parameter for updating the correlation parameter (if point.ind==TRUE)}
#'	\item{"MH.nu"}{MH tuning parameter for Nu parameters (Langevin-Hastings multivariate update)}
#'	\item{"MH.beta"}{A vector of tuning parameters for betas of the abundance process (dimension = number of columns of habitat DM)}
#'	\item{"grp.mean"}{Expected group size minus 1 ; used for prior & psedo-prior for RJMCMC updates}
#'	\item{"RJ.N"}{A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect}
#' @param DM.hab	A design matrix for the log of abundance intensity
#' @param DM.det	A design matrix for the probit of detection probability
#' @param Q			An inverse precision matrix for the spatial ICAR process
#' @param n.transects	Number of transects
#' @param S				Number of strata cells
#' @param Accept		A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms
#' @param Prior.pars	A list object giving parameters of prior distribution.  Includes the following slots
#'	\item{"a.eta"}{alpha parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))}
#'  \item{"b.eta"}{beta parameter for prior precision of spatial process (assumed Gamma(a.eta,b.eta))}
#'	\item{"a.nu"}{alpha parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))}
#'	\item{"b.nu"}{beta parameter for prior precision of overdispersion process (assumed Gamma(a.nu,b.nu))} 
#'	\item{"beta.sd"}{standard deviation for regression coefficients (assumed Normal(0,beta.sd^2)}
#' @return returns a list with the following slots: 
#' 	\item{MCMC}{A list object containing posterior samples}
#'  \item{Accept}{A list object indicating the number of proposals that were accepted for parameters updated via Metropolis- or Langevin-Hastings algorithms
#'  \item{Control}{A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' #' @export
#' @import Matrix
#' @keywords areal model, data augmentation, distance sampling, mcmc, reversible jump
#' @author Paul B. Conn
require(mvtnorm)
require(Matrix)
require(truncnorm)

##################################################
############   MCMC Algorithm ####################
##################################################
mcmc_ds<-function(Par,Data,cur.iter,adapt,Control,DM.hab,DM.det,Q,n.transects,S,Accept,Prior.pars){	
	for(iiter in 1:cur.iter){
		#if((iiter%%1000)==1)cat(paste('\n ', iiter))
		#cat(paste('\n ', iiter))
		
		########## update abundance parameters at the strata scale   ################
		
		#update nu parameters (log lambda)
		Mu=DM.hab%*%Par$hab+Par$Eta
		Grad1=sapply(Lam.index,'log_lambda_gradient',Mu=Mu,Nu=Par$Nu,N=Par$G,var.nu=1/Par$tau.nu)
		Prop=Par$Nu+Control$MH.nu^2*0.5*Grad1+Control$MH.nu*rnorm(S)
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
		
		if(spat.ind==0){
			#update eta parameters (spatial random effects)
			V.eta.inv <- Par$tau.nu*diag(S) + Par$tau.eta*Q
			M.eta <- solve(V.eta.inv, Par$tau.nu*(Par$Nu-DM.hab%*%Par$hab))		
			Par$Eta<-as.vector(M.eta+solve(chol(V.eta.inv),rnorm(S,0,1)))
			Par$Eta=Par$Eta-mean(Par$Eta)  #centering
			
			#update tau_eta  (precision of spatial process)
			Par$tau.eta <- rgamma(1, (S-1)/2 + Prior.pars$a.eta, as.numeric(crossprod(Par$Eta, Q %*% Par$Eta)/2) + Prior.pars$b.eta)
		}
		
		#update tau_nu	 (precision for Poisson overdispersion)
		Mu=DM.hab%*%Par$hab+Par$Eta
		Diff=Par$Nu-Mu
		Par$tau.nu <- rgamma(1,S/2 + Prior.pars$a.nu, as.numeric(crossprod(Diff,Diff))/2 + Prior.pars$b.nu)
		
		
		#translate to lambda scale
		Lambda=Area.hab*exp(Par$Nu)
		Lambda.trans=Lambda[Mapping]
		
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
		Par$G=rpois(S,Lambda*(1-Covered.area))
		Par$N=Par$G+rpois(S,Control$grp.mean*Par$G) #add the Par$G since number in group > 0 
		Par$G[Mapping]=Par$G[Mapping]+G.transect
		Par$N[Mapping]=Par$N[Mapping]+N.transect
		
		########## update abundance, distances, ind. covariates for observed transects using RJMCMC  #############
		for(itrans in 1:n.transects){
			Sample=c(-Control$RJ.N[itrans]:Control$RJ.N[itrans])
			Sample=Sample[-(Control$RJ.N[itrans]+1)] #don't make proposals where pop size stays the same
			a=sample(Sample,1)
			Sigma=matrix(Par$cor,n.Observers[itrans],n.Observers[itrans])
			diag(Sigma)=1
			offdiag=which(Sigma!=1)
			
			if(a>0){ # proposal an addition
				if(((G.transect[itrans]+a)*n.Observers[itrans])>M)cat('\n Error: proposed abundance > M; increase M value! \n')
				else{
					Cur.dat=Data[itrans,(n.Records[itrans]+1):(n.Records[itrans]+a*n.Observers[itrans]),]
					X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
					ExpY=X%*%Par$det
					P=c(1:a)
					for(i in 1:a){
						Sigma[offdiag]=Par$cor*(i.binned*(Cur.dat[i*n.Observers[itrans],dist.pl]-1)*dist.mult+(1-i.binned)*Cur.dat[i*n.Observers[itrans],dist.pl])
						P[i]=pmvnorm(upper=c(0,0),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:a){
						tmp.sum=tmp.sum-log(G.transect[itrans]-G.obs[itrans]+i)+log(P[i])
					}
					MH.prob=exp(a*log(Lambda.trans[itrans]*Area.trans[itrans])+tmp.sum)
					if(runif(1)<MH.prob){
						G.transect[itrans]=G.transect[itrans]+a
						n.Records[itrans]=G.transect[itrans]*n.Observers[itrans]
						if(grps==FALSE)N.transect[itrans]=G.transect[itrans]
						else N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(stacked.names=="Group")])/n.Observers[itrans]
						Accept$N[itrans]=Accept$N[itrans]+1
						#generate Y-tilde values
						Tmp=matrix(ExpY,a,n.Observers[itrans],byrow=TRUE)
						Y.tilde.temp=rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,1],sd=1)
						if(n.Observers[itrans]>1){
							Dist=matrix(Cur.dat[,dist.pl],a,n.Observers[itrans],byrow=TRUE)
							Cor=Par$cor*(i.binned*(Dist[,1]-1)*dist.mult+(1-i.binned)*Dist[,1])
							Y.tilde.temp=cbind(Y.tilde.temp,rtruncnorm(n=a,a=-Inf,b=0,mean=Tmp[,2]+Cor*(Y.tilde.temp-Tmp[,2]),sd=sqrt(1-Cor^2)))
						}
						Y.tilde[(n.Records[itrans]-a*n.Observers[itrans]+1):n.Records[itrans],itrans]=as.vector(t(Y.tilde.temp))
					}	
				}
			}
			else{ #proposal a deletion
				if((G.transect[itrans]+a)>=G.obs[itrans]){ #can only delete records where an animal isn't observed
					Cur.dat=Data[itrans,n.Records[itrans]:(n.Records[itrans]+a*n.Observers[itrans]+1),]
					X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
					ExpY=X%*%Par$det
					P=c(1:-a)
					for(i in 1:-a){
						Sigma[offdiag]=Par$cor*(i.binned*(Cur.dat[i*n.Observers[itrans],dist.pl]-1)*dist.mult+(1-i.binned)*Cur.dat[i*n.Observers[itrans],dist.pl])
						P[i]=pmvnorm(upper=rep(0,n.Observers[itrans]),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:-a){
						tmp.sum=tmp.sum+log(G.transect[itrans]-G.obs[itrans]-i+1)-log(P[i])
					}
					MH.prob=exp(a*log(Lambda.trans[itrans]*Area.trans[itrans])+tmp.sum)
					if(runif(1)<MH.prob){
						G.transect[itrans]=G.transect[itrans]+a
						n.Records[itrans]=G.transect[itrans]*n.Observers[itrans]
						Accept$N[itrans]=Accept$N[itrans]+1
					}					
					
				}
			}
			
			#update distances, individual covariates for latent animals not currently in the population
			#fill distances for unobserved
			if(i.binned==1)Data[itrans,(n.Records[itrans]+1):M,dist.pl]=rep(sample(c(1:n.bins),size=(M-n.Records[itrans])/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
			else Data[itrans,(n.Records[itrans]+1):M,dist.pl]=rep(runif((M-n.Records[itrans])/n.Observers[itrans]),each=n.Observers[itrans])
			#fill individual covariate values for (potential) animals that weren't observed
			if(n.ind.cov>0){
				for(icov in 1:n.ind.cov){
					rsamp=switch_sample(n=(M-n.Records[itrans])/n.Observers[itrans],pdf=Cov.prior.pdf[icov],cur.par=Cov.prior.parms[,icov])
					Data[itrans,(n.Records[itrans]+1):M,dist.pl+icov]=rep(rsamp,each=n.Observers[itrans])
				}
			}
			
			#update distances, individual covariates for animals that ARE in the population but never observed
			cur.G=G.transect[itrans]-G.obs[itrans]
			if(G.transect[itrans]>G.obs[itrans]){
				#distance
				if(i.binned==1)dist.star=sample(c(1:n.bins),cur.G,replace=TRUE)
				else dist.star=runif(cur.G)
				Cur.dat=Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),]
				cur.dist=Cur.dat[,dist.pl]
				X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
				ExpY=X%*%Par$det
				L.old=c(1:length(dist.star))
				Tmp.Y.tilde=Y.tilde[(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),itrans]
				for(i in 1:length(dist.star)){
					Sigma[offdiag]=Par$cor*(i.binned*(cur.dist[i]-1)*dist.mult+(1-i.binned)*cur.dist[i])
					L.old[i]=dmvnorm(Tmp.Y.tilde[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
				}	
				Cur.dat[,dist.pl]=rep(dist.star,each=n.Observers[itrans])
				X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
				ExpY=X%*%Par$det
				L.star=L.old
				for(i in 1:length(dist.star)){
					Sigma[offdiag]=Par$cor*(i.binned*(dist.star[i]-1)*dist.mult+(1-i.binned)*dist.star[i])
					L.star[i]=dmvnorm(Tmp.Y.tilde[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
				}	
				Acc=(runif(length(L.star))<(L.star/L.old))
				Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),dist.pl]=(1-rep(Acc,each=n.Observers[itrans]))*Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),dist.pl]+rep(Acc,each=n.Observers[itrans])*rep(dist.star,each=n.Observers[itrans])
				
				#individual covariates
				if(n.ind.cov>0){
					for(icov in 1:n.ind.cov){
						Cov.star=switch_sample(n=cur.G,pdf=Cov.prior.pdf[icov],cur.par=Cov.prior.parms[,icov])
						Cur.dat=Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),]
						cur.dist=Cur.dat[,dist.pl]
						X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
						ExpY=X%*%Par$det
						L.old=c(1:length(Cov.star))
						for(i in 1:length(Cov.star)){
							Sigma[offdiag]=Par$cor*(i.binned*(cur.dist[i]-1)*dist.mult+(1-i.binned)*cur.dist[i])
							L.old[i]=dmvnorm(Tmp.Y.tilde[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
						}	
						Cur.dat[,dist.pl+icov]=rep(Cov.star,each=n.Observers[itrans])
						X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
						ExpY=X%*%Par$det
						L.star=L.old
						for(i in 1:length(Cov.star)){
							Sigma[offdiag]=Par$cor*(i.binned*(cur.dist[i]-1)*dist.mult+(1-i.binned)*cur.dist[i])
							L.star[i]=dmvnorm(Tmp.Y.tilde[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
						}	
						Acc=(runif(length(L.star))<(L.star/L.old))
						Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),dist.pl+icov]=(1-rep(Acc,each=n.Observers[itrans]))*Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),dist.pl+icov]+rep(Acc,each=n.Observers[itrans])*rep(Cov.star,each=n.Observers[itrans])						
					}
				}
			}
			if(grps==FALSE)N.transect[itrans]=G.transect[itrans]
			else N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(stacked.names=="Group")])/n.Observers[itrans]
			
			
			#update Y.tilde
			Cur.dat=Data[itrans,1:n.Records[itrans],]
			X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
			ExpY=matrix(X%*%Par$det,G.transect[itrans],n.Observers[itrans],byrow=TRUE)
			Temp.Y.tilde=matrix(Y.tilde[1:n.Records[itrans],itrans],G.transect[itrans],n.Observers[itrans],byrow=TRUE)
			Dist=matrix(Cur.dat[,dist.pl],G.transect[itrans],n.Observers[itrans],byrow=TRUE)
			Cor=Par$cor*(i.binned*(Dist[,1]-1)*dist.mult+(1-i.binned)*Dist[,1])
			Resp=matrix(Cur.dat[,2],G.transect[itrans],n.Observers[itrans],byrow=TRUE)
			EY1=ExpY[,1]+Cor*(Temp.Y.tilde[,2]-ExpY[,2])
			Temp.Y.tilde[,1] <- rtruncnorm(G.transect[itrans], a=ifelse(Resp[,1]==0,-Inf,0), b=ifelse(Resp[,1]==0,0,Inf), EY1, sqrt(1-Cor^2))
			EY2=ExpY[,2]+Cor*(Temp.Y.tilde[,1]-ExpY[,1])
			Temp.Y.tilde[,2] <- rtruncnorm(G.transect[itrans], a=ifelse(Resp[,2]==0,-Inf,0), b=ifelse(Resp[,2]==0,0,Inf), EY2, sqrt(1-Cor^2))
			Y.tilde[1:n.Records[itrans],itrans]=as.vector(t(Temp.Y.tilde))
		}
		
		###############       update detection process parameters       ##############
		# First, assemble stacked adjusted Response, X matrices across all transects; 
		#basic form of response is Ytilde[obs1]-cor*Ytilde[obs2]
		#adjusted DM rows are of form X[obs1]-cor*X[obs2]
		Cur.dat=Data[1,1:n.Records[1],]
		Dist=matrix(Cur.dat[,dist.pl],G.transect[1],n.Observers[1],byrow=TRUE)[,1]
		X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)),dim=c(n.beta.det,2,G.transect[1]))
		if(n.Observers[1]==2){
			Tmp.cor=Par$cor*(i.binned*(Dist-1)*dist.mult+(1-i.binned)*Dist)
			Cor=rep(Tmp.cor,2)  #assemble vector of correlation parameters for each observation
			X.beta=rbind(t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
			Y.temp=matrix(Y.tilde[1:n.Records[1],1],G.transect[1],2,byrow=TRUE)
			Y.beta=c(Y.temp[,1]-Tmp.cor*Y.temp[,2],Y.temp[,2]-Tmp.cor*Y.temp[,1])
		}
		else{
			X.beta=t(X.temp[,1,])
			Y.beta=Y.tilde[1:n.Records[1],1]
			Cor=rep(1,n.Records[1])
		}
		for(itrans in 2:n.transects){
			Cur.dat=Data[itrans,1:n.Records[itrans],]
			Dist=matrix(Cur.dat[,dist.pl],G.transect[itrans],n.Observers[itrans],byrow=TRUE)[,1]
			X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)),dim=c(n.beta.det,2,G.transect[itrans]))
			if(n.Observers[1]==2){
				Tmp.cor=Par$cor*(i.binned*(Dist-1)*dist.mult+(1-i.binned)*Dist)
				Cor=c(Cor,rep(Tmp.cor,2))  #assemble vector of correlation parameters for each observation
				X.beta=rbind(X.beta,t(X.temp[,1,])-Tmp.cor*t(X.temp[,2,]),t(X.temp[,2,])-Tmp.cor*t(X.temp[,1,]))
				Y.temp=matrix(Y.tilde[1:n.Records[itrans],itrans],G.transect[itrans],2,byrow=TRUE)
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
		if(point.ind==1){
			cor.star=Par$cor+runif(1,-Control$MH.cor,Control$MH.cor)
			if(cor.star>0 & cor.star<1){
				I.gt.one=which(n.Observers>1)
				n.gt.one=length(I.gt.one)
				
				Cur.dat=Data[I.gt.one[1],1:n.Records[I.gt.one[1]],]
				Dist=matrix(Cur.dat[,dist.pl],G.transect[I.gt.one[1]],n.Observers[I.gt.one[1]],byrow=TRUE)[,1]
				X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)),dim=c(n.beta.det,2,G.transect[I.gt.one[1]]))
				Y.temp=matrix(Y.tilde[1:n.Records[I.gt.one[1]],I.gt.one[1]],G.transect[I.gt.one[1]],2,byrow=TRUE)			
				Delta1=Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det)
				Delta2=Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det)
				
				if(n.gt.one>1){
					for(itrans in 2:n.gt.one){
						Cur.dat=Data[I.gt.one[itrans],1:n.Records[I.gt.one[itrans]],]
						Dist=c(Dist,matrix(Cur.dat[,dist.pl],G.transect[I.gt.one[itrans]],n.Observers[I.gt.one[itrans]],byrow=TRUE)[,1])
						X.temp=array(t(get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)),dim=c(n.beta.det,2,G.transect[I.gt.one[itrans]]))
						Y.temp=matrix(Y.tilde[1:n.Records[I.gt.one[itrans]],I.gt.one[itrans]],G.transect[I.gt.one[itrans]],2,byrow=TRUE)			
						Delta1=c(Delta1,Y.temp[,1]-(t(X.temp[,1,]) %*% Par$det))
						Delta2=c(Delta2,Y.temp[,2]-(t(X.temp[,2,]) %*% Par$det))
					}
				}
				Cor=Par$cor*(i.binned*(Dist-1)*dist.mult+(1-i.binned)*Dist)
				logP.old=-.5*(sum(log(1-Cor^2))+sum((Delta1^2+Delta2^2-2*Cor*Delta1*Delta2)/(1-Cor^2)))
				Cor=cor.star*(i.binned*(Dist-1)*dist.mult+(1-i.binned)*Dist)
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
				for(ipar in 1:n.transects){
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
			ttc <- round((iiter-15)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}
		
	}
	Out=list(MCMC=MCMC,Accept=Accept,Control=Control)
	Out
}

