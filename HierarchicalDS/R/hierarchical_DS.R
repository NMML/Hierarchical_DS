
#hierarchical_DS<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,Bin.length,Hab.cov,n.obs.cov=0,Det.cov=NULL,pol.eff=c(1:2),point.ind=TRUE,Inits=NULL,grps=FALSE,M,Control,adapt=TRUE)
	#Dat - matrix or data frame with the following columns
		#1-transect ID
		#2-match number
		#3-Observer ID
		#4-Observation (0/1)
		#5-x - Observer covariates (these can vary among transects but are assumed constant within transects)
		#x+1-Distance (if all integers, assumed to be discrete bins; if continuous, assumed standardized to (0,1) interval)
		#x+2-?? Group size and other individual covariates thought to influence detection; if group size is one of them, it's assumed to be column x+2
			#Note that covariates at transect level need to be specified in "Det.cov"
			#Note: column names are used to tag covariates
	#Adj - adjacency matrix for habitat cells (diagonal matrix implies spatial independence)
	#Area.hab - a vector giving the area of each geographical strata (default is equal area)
	#Mapping - A vector giving the habitat cell id # for each transect
	#Area.trans - a vector giving the effective area covered by each transect as fraction of total area in the strata it is located
	#Bin.length - If distances are binned, this vector gives the relative length of each distance bin (vector must sum to one)
	#n.obs.cov - Number of observer covariates (e.g., seat position, visibility, etc.)
	#Hab.cov - a data.frame object giving covariates thought to influence abundance intensity at strata level; column names index individual covariates
	#Hab.formula - a formula object giving the specific model for abundance intensity at the strata level (e.g., ~Vegetation+Latitude)
	#Det.formula - a formula giving the model for detection probability (e.g. ~Distance+Group+Visibility+Observer). Note that
		#there are several "reserved" variable names.  "Distance", "Observer", and "Group" are reserved variable names.
	#Cov.prior.pdf - if individual covariates are provided, this character vector gives the form of the prior pdfs for each covariate
				#current possibilities are "poisson", "pois1","uniform.disc", "uniform.cont", or "normal".
				#"pois1" is the 1+x where x~poisson
	#Cov.prior.parms - a (2 X n) matrix where n is the number of individual covariates (other than distance).  Each column
			#gives the parameters associated with the prior pdf of a given covariate (only the value in the first row is
			#used for "poisson"; for normal, first row gives mean, second row gives sd; for uniform, first row gives lower,
 			#second row gives upper; for constant, the parameter entries are just placeholders (no parameters are required)
		    #note that these priors are also used to propose covariate values during RJMCMC, so should be 
			#made to be biologically plausible (i.e., don't use 'vague' priors!)
	#pol.eff - for continuous distance, which polynomial degrees to model (default is c(1:2); an intercept is always estimated when "Distance" is listed in "Det.formula")
	#point.ind - estimate a correlation parameter for detection probability that's an increasing function of distance?
	#grps - if FALSE, detections are assumed to all be of individual animals
	#M - maximum possible value for number of groups present in any one transect (in practice just set high enough that values at M and above are never sampled during MCMC)
	#			and can be fine tuned as needed
	#Control - a list object including the following slots:
		#"iter" - number of MCMC iterations
		#"burnin" - number of MCMC burnin iterations
		#"adapt" - if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run 
		#"MH.cor" - Metropolis-hastings tuning parameter for updating the correlation parameter (if point.ind==TRUE)
		#"RJ.N" - A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect
	#Inits - an (optional) list object providing initial values for model parameters
		#"Beta.hab" - initial values for habitat linear predictor parameters
		#"Beta.det" - initial values for detection model (includes distance, observer, env. variables, and individual covariates)
		#"cor.par" - if point.ind==TRUE, this is an initial value for the correlation parameter (which must be in (0,1))	
	#adapt - if adapt==TRUE, run an additional Control$adapt number of MCMC iterations to optimize MCMC proposal distributions prior to primary MCMC

	####start by defining stuff directly (just for debugging...switch to outside function later)
	require(mvtnorm)
	source("z:/git/New folder/hierarchicalDS/R/simulate_data.R")
	S=10 #number of sites (# transects=# habitat areas to start with)
	Dat=simulate_data(S=S) 
	n.bins=length(unique(Dat[,6]))
	Adj=diag(S)
	Area.hab=1
	Mapping=c(1:S)
	Area.trans=rep(1,S)
	Bin.length=rep(1,n.bins)
	Hab.cov=data.frame(matrix(c(1:S)/S,S,1))
	colnames(Hab.cov)="Cov1"
	n.obs.cov=1  #dummy 'seat' variable
	Hab.formula=~Cov1
	Det.formula=~Observer+Seat+Distance+Group
	Cov.prior.pdf="pois1"
	Cov.prior.parms=matrix(c(3,NA),2,1)
	colnames(Cov.prior.parms)="Group"
	pol.eff=2 #not currently used since using distance bins
	point.ind=TRUE
	grps=TRUE
	M=200
	Control=list(iter=15000,burnin=0,MH.cor=0.05,MH.dist=0.1,grp.mean=3,RJ.N=rep(5,S),adapt=1000)
	Inits=list(hab=c(0,100),det=c(1.2,-.2,-.4,0,-.8,-1.4,-1.8,-2,.2),cor=0.5) #start at true values during debugging
	#note, these differ from beta in simulate data because of the way DM is implemented to sim data

	adapt=FALSE
	
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
	M=M*n.obs.max
	n.hab.cov=ifelse(is.null(Hab.cov)==1,0,ncol(Hab.cov))
	n.ind.cov=ncol(Dat)-(5+n.obs.cov) #number of individual covariates
	N.remain=Adj*0  #keep track of "remaining" abundance in each strata for portion of area not covered by transects
	N.total=N.remain #keep track of total abundance by strata

	var.names=colnames(Dat)
	factor.ind=sapply(Dat[1,],is.factor)
	#Check to make sure input values are internally consistent
	#Later...
	

	switch_sample<-function(n,pdf,cur.par){
		switch(pdf,
			pois1=rpois(n,cur.par[1])+1,
			poisson=rpois(n,cur.par[1]),
			normal=rnorm(n,cur.par[1],cur.par[2]),
			unif.disc=sample(cur.par[1]:cur.par[2],n,replace=TRUE),
			unif.cont=runif(n,cur.par[1],cur.par[2])
		)
	}

	#Initialize data augmentation multi-d array ("Data"), parameter vectors and matrices
	Data<-array(0,dim=c(n.transects,M,3+n.obs.cov+n.ind.cov)) #array cols are Obs ID,Y,Obs covariates,Distance,Ind covariates
	n.Observers=rep(0,n.transects)	#holds number of observers in each transect
	G.transect=rep(0,n.transects)  #number of groups by transect
	n.Records=rep(0,n.transects) #number of records by transect (=G.transect*n.Observers)
	N.transect=G.transect #total abundance by transect
	Dat.num=Dat
	for(icol in 1:ncol(Dat)){
		Dat.num[,icol]=as.numeric(Dat[,icol])
	}
	Dat.num=as.matrix(Dat.num)
	for(itrans in 1:n.transects){
		Cur.dat=Dat.num[which(Dat.num[,1]==itrans),2:ncol(Dat.num)]
		first.match=Cur.dat[1,1]
		n.Observers[itrans]=length(which(Cur.dat[,1]==first.match))
		Cur.dat=Cur.dat[,-1]
		Data[itrans,1:nrow(Cur.dat),1:(2+n.obs.cov)]=Cur.dat[,1:(2+n.obs.cov)]
		Data[itrans,1:nrow(Cur.dat),3+n.obs.cov]=Cur.dat[,"Distance"]
		#fill distances for unobserved
		if(i.binned==1)Data[itrans,(nrow(Cur.dat)+1):M,3+n.obs.cov]=rep(sample(c(1:n.bins),size=(M-nrow(Cur.dat))/n.Observers[itrans],replace=TRUE,prob=Bin.length),each=n.Observers[itrans])
		else Data[itrans,(nrow(Cur.dat)+1):M,3+n.obs.cov]=rep(runif((M-nrow(Cur.dat))/n.Observers[itrans]),each=n.Observers[itrans])
		#fill individual covariate values for (potential) animals that weren't observed
		if(n.ind.cov>0){
			Data[itrans,1:nrow(Cur.dat),(4+n.obs.cov):(4+n.obs.cov+n.ind.cov-1)]=Cur.dat[,(4+n.obs.cov):(4+n.obs.cov+n.ind.cov-1)]
			for(icov in 1:n.ind.cov){
				rsamp=switch_sample(n=(M-nrow(Cur.dat))/n.Observers[itrans],pdf=Cov.prior.pdf[icov],cur.par=Cov.prior.parms[,icov])
				Data[itrans,(nrow(Cur.dat)+1):M,3+n.obs.cov+icov]=rep(rsamp,each=n.Observers[itrans])
			}
		}
		G.transect[itrans]=nrow(Cur.dat)/n.Observers[itrans]		#initialize abundance in each transect area to be = to total number of animals observed
		n.Records[itrans]=G.transect[itrans]*n.Observers[itrans]
		N.transect[itrans]=ifelse(grps==FALSE,G.transect[itrans],sum(Cur.dat[,which(colnames(Cur.dat)=="Group")])/n.Observers[itrans])
		#fill observer ids
		Data[itrans,(n.Records[itrans]+1):M,1]=rep(Data[itrans,1:n.Observers[itrans],1],(M-n.Records[itrans])/n.Observers[itrans])
		curcol=3
		#fill observer covariates
		if(n.obs.cov>0){
			for(icol in curcol:(curcol+n.obs.cov-1)){
				Data[itrans,(n.Records[itrans]+1):M,curcol]=rep(Data[itrans,1:n.Observers[itrans],curcol],(M-n.Records[itrans])/n.Observers[itrans])
			}
			curcol=curcol+n.obs.cov
		}
	}
	G.obs=G.transect
	
	stacked.names=colnames(Dat)[3:ncol(Dat)]
	
	stack_data<-function(Data,G.transect,n.transects,stacked.names,factor.ind){
		#convert from "sparse" 3-d data augmentation array to a rich 2-d dataframe for updating beta parameters 
		Stacked=as.data.frame(Data[1,1:G.transect[1],])
		for(itrans in 2:n.transects){
			Stacked=rbind(Stacked,Data[itrans,1:G.transect[itrans],])
		}
		colnames(Stacked)=stacked.names	#gotta reestablish variable type since 3-d array doesn't hold it
		factor.cols=which(factor.ind[stacked.names]==TRUE) 
		if(length(factor.cols)>0){
			for(icol in 1:length(factor.cols)){
				Stacked[,factor.cols[icol]]=as.factor(Stacked[,factor.cols[icol]])
			}
		}
		Stacked
	}
	Stacked=stack_data(Data,G.transect,n.transects,stacked.names,factor.ind) #a stacked form of data for updating beta parameters
	
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
	
	get_mod_matrix<-function(Cur.dat,stacked.names,factor.ind,Det.formula,Levels){
		Cur.dat=as.data.frame(Cur.dat)
		colnames(Cur.dat)=stacked.names
		factor.cols=which(factor.ind[stacked.names]==TRUE)
		if(length(factor.cols)>0){
			for(icol in 1:length(factor.cols)){
				Cur.dat[,factor.cols[icol]]=eval(parse(text=paste('factor(Cur.dat[,factor.cols[icol]],levels=Levels$',names(factor.cols)[icol],')',sep='')))
			}
		}
		DM=model.matrix(Det.formula,data=Cur.dat)
	}
	
	
	generate_inits<-function(DM.hab,DM.det,G.transect,Area.trans,Area.hab){
		Par=list(det=rnorm(ncol(DM.det),0,1),hab=rep(0,ncol(DM.hab)),cor=runif(1,0,1))
		Par$hab[1]=mean(G.transect)/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(1,0,1))
		Par
	}
	
	Par=Inits
	if(is.null(Inits)==TRUE)Par=generate_inits(DM.hab,DM.det,G.transect,Area.trans,Area.hab)	
			
	mcmc.length=Control$iter-Control$burnin
	MCMC=list(N=rep(0,mcmc.length),G=matrix(0,mcmc.length,S),Hab=data.frame(matrix(0,mcmc.length,length(Par$hab))),Det=data.frame(matrix(0,mcmc.length,length(Par$det))),cor.par=rep(0,mcmc.length))
	colnames(MCMC$Hab)=colnames(DM.hab)
	colnames(MCMC$Det)=colnames(DM.det)
	
	Y.tilde=matrix(0,n.transects,M)
	dist.pl=3+n.obs.cov
	if(i.binned==0)dist.mult=1
	if(i.binned==1)dist.mult=1/(n.bins-1)
	
	for(iiter in 1:Control$iter){
		#if((iiter%%1000)==1)cat(paste('\n ', iiter))
		cat(paste('\n ', iiter))
		Accept=list(cor.par=0,N=rep(0,n.transects))
		
		### update abundance parameters at the strata scale
		
		Lambda=DM.hab%*%Par$hab  #keep these at expected values for the moment
		
		### update abundance, distances, ind. covariates using RJMCMC
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
						Sigma[offdiag]=i.binned*Par$cor*(Cur.dat[i*n.Observers[itrans],dist.pl]-1)*dist.mult+(1-i.binned)*Cur.dat[i*n.Observers[itrans],dist.pl]
						P[i]=pmvnorm(upper=c(0,0),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:a){
						tmp.sum=tmp.sum-log(G.transect[itrans]-G.obs[itrans]+i)+log(P[i])
					}
					MH.prob=exp(a*log(Lambda[itrans]*Area.trans[itrans])+tmp.sum)
					if(runif(1)<MH.prob){
						G.transect[itrans]=G.transect[itrans]+a
						n.Records[itrans]=G.transect[itrans]*n.Observers[itrans]
						if(grps==FALSE)N.transect[itrans]=G.transect[itrans]
						else N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(stacked.names=="Group")])/n.Observers[itrans]
						Accept$N[itrans]=Accept$N[itrans]+1
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
						Sigma[offdiag]=i.binned*Par$cor*(Cur.dat[i*n.Observers[itrans],dist.pl]-1)*dist.mult+(1-i.binned)*Cur.dat[i*n.Observers[itrans],dist.pl]
						P[i]=pmvnorm(upper=rep(0,n.Observers[itrans]),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
					}
					tmp.sum=0
					for(i in 1:-a){
						tmp.sum=tmp.sum+log(G.transect[itrans]-G.obs[itrans]-i+1)-log(P[i])
					}
					MH.prob=exp(a*log(Lambda[itrans]*Area.trans[itrans])+tmp.sum)
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
				else dist.star=runif(1,cur.G)
				Cur.dat=Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),]
				cur.dist=Cur.dat[,dist.pl]
				X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
				ExpY=X%*%Par$det
				L.old=c(1:length(dist.star))
				for(i in 1:length(dist.star)){
					Sigma[offdiag]=i.binned*Par$cor*(cur.dist[i]-1)*dist.mult+(1-i.binned)*cur.dist[i]
					L.old[i]=pmvnorm(upper=rep(0,n.Observers[itrans]),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
				}	
				Cur.dat[,dist.pl]=rep(dist.star,each=n.Observers[itrans])
				X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
				ExpY=X%*%Par$det
				L.star=L.old
				for(i in 1:length(dist.star)){
					Sigma[offdiag]=i.binned*Par$cor*(dist.star[i]-1)*dist.mult+(1-i.binned)*dist.star[i]
					L.star[i]=pmvnorm(upper=rep(0,n.Observers[itrans]),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
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
							Sigma[offdiag]=i.binned*Par$cor*(cur.dist[i]-1)*dist.mult+(1-i.binned)*cur.dist[i]
							L.old[i]=pmvnorm(upper=rep(0,n.Observers[itrans]),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
						}	
						Cur.dat[,dist.pl+icov]=rep(Cov.star,each=n.Observers[itrans])
						X=get_mod_matrix(Cur.dat=Cur.dat,stacked.names,factor.ind,Det.formula,Levels)
						ExpY=X%*%Par$det
						L.star=L.old
						for(i in 1:length(Cov.star)){
							Sigma[offdiag]=i.binned*Par$cor*(cur.dist[i]-1)*dist.mult+(1-i.binned)*cur.dist[i]
							L.star[i]=pmvnorm(upper=rep(0,n.Observers[itrans]),mean=ExpY[(i*n.Observers[itrans]-n.Observers[itrans]+1):(i*n.Observers[itrans])],sigma=Sigma)
						}	
						Acc=(runif(length(L.star))<(L.star/L.old))
						Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),dist.pl+icov]=(1-rep(Acc,each=n.Observers[itrans]))*Data[itrans,(G.obs[itrans]*n.Observers[itrans]+1):(G.transect[itrans]*n.Observers[itrans]),dist.pl+icov]+rep(Acc,each=n.Observers[itrans])*rep(Cov.star,each=n.Observers[itrans])						
					}
				}
			}
			if(grps==FALSE)N.transect[itrans]=G.transect[itrans]
			else N.transect[itrans]=sum(Data[itrans,1:n.Records[itrans],which(stacked.names=="Group")])/n.Observers[itrans]
			
			
			#update Y.tilde
			EY1=Mu1+Cor*(Y.tilde[,2]-Mu2)
			Y.tilde[,1] <- rtruncnorm(N, a=ifelse(Dat[,3]==0,-Inf,0), b=ifelse(Dat[,3]==0,0,Inf), EY1, sqrt(1-Cor^2))
			EY2=Mu2+Cor*(Y.tilde[,1]-Mu1)
			Y.tilde[,2] <- rtruncnorm(N, a=ifelse(Dat[,4]==0,-Inf,0), b=ifelse(Dat[,4]==0,0,Inf), EY2, sqrt(1-Cor^2))



		}
		
		#store results if applicable
		if(iiter>Control$burnin){
			MCMC$G[iiter-Control$burnin,]=G.transect
			MCMC$N[iiter-Control$burnin]=sum(N.transect)
			MCMC$cor[iiter-Control$burnin]=Par$cor
			MCMC$Hab[iiter-Control$burnin,]=Par$hab
			MCMC$Det[iiter-Control$burnin,]=Par$det
		}
		#update betas for detection process
		
		#update correlation parameter for detection process (if applicable)


	}

#}

	