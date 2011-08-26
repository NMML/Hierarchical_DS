
#hierarchical_DS<-function(Dat,Adj,Area.hab=1,Mapping,Area.trans,Bin.length,DM.hab,Det.cov=NULL,obs.eff=TRUE,pol.eff=c(1:2),point.ind=TRUE,grps=FALSE,M,Control,adapt=TRUE)
	#Dat - matrix or data frame with the following columns
		#1-transect ID
		#2-First observer ID
		#3-Second observer ID
		#4-Observation 1 - (0/1)
		#5-Observation 2 - (0/1)
		#6-Distance (if all integers, assumed to be discrete bins
		#7-Group size 
	#Adj - adjacency matrix for habitat cells (diagonal matrix implies spatial independence)
	#Area.hab - a vector giving the area of each habitat cell (default is equal area)
	#Mapping - A vector giving the habitat cell id # for each transect
	#Area.trans - a vector giving the effective area covered by each transect as fraction of total area in the strata it is located
	#Bin.length - If distances are binned, this vector gives the relative length of each distance bin (vector must sum to one)
	#DM.hab - a data.frame object representing a design matrix for abundance intensity at strata level; column names index individual covariates
	#Det.cov - a data frame giving covariates thought to influence detection probability for each transect 
	#		(each row gives values for a separate transect). Note that distance and observer effects are handled internally
	#		(i.e., shouldn't be provided here).
	#obs.eff - include an offset for observer effect on det prob?
	#pol.eff - for continuous distance, which polynomial degrees to model (default is c(1:2); an intercept is always estimated)
	#point.ind - estimate a correlation parameter for detection probability that's an increasing function of distance?
	#grps - if FALSE, detections are assumed to all be of individual animals
	#M - maximum possible value for abundance in any one transect (in practice just set high enough that values at M and above are never sampled during MCMC)
	#			and can be fine tuned as needed
	#Control - a list object including the following slots:
		#"iter" - number of MCMC iterations
		#"burnin" - number of MCMC burnin iterations
		#"adapt" - if adapt==TRUE, this gives the number of additional MCMC iterations should be performed to adapt MCMC proposals to optimal ranges prior to final MCMC run 
		#"MH.cor" - Metropolis-hastings tuning parameter for updating the correlation parameter (if point.ind==TRUE)
		#"MH.dist" - M-H tuning parameter for updating individual distances (continuous distances only)
		#"grp.mean" - approximate mean value for group size (if all groups are detected); this value is simply used as a starting value 
		#"RJ.N" - A vector giving the maximum number of additions and deletions proposed in an iteration of the RJMCMC algorithm for each transect
	#Inits - a list object including the following initial values:
		#"Beta.hab" - a vector of initial values for parameters influencing animal distribution (should have length=number of columns of DM.hab)
		#"Beta.det" - a vector of initial values for parameters influencing detection probability (should have length=number of columns of Det.cov)
					#note that initial values for parameters describing the effects of distance, group, observers are handled internally
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
	DM.hab=data.frame(matrix(c(1:S)/S,S,1))
	colnames(DM.hab)="Cov1"
	Det.cov=NULL
	obs.eff=TRUE
	pol.eff=2 #not currently used since using distance bins
	point.ind=TRUE
	grps=TRUE
	M=200
	Control=list(iter=15000,burnin=0,MH.cor=0.05,MH.dist=0.1,grp.mean=3,RJ.N=rep(5,S),adapt=1000)
	Inits=list(Beta.hab=100,Beta.det=NULL,cor.par=0.4)    #p1 and p2 are probability of detection on the line
	adapt=FALSE
	
	
	n.observers=length(unique(c(Dat[,2],Dat[,3])))
	n.transects=length(Area.trans)
	i.binned=0
	if(sum(Dat[,6]%%1)==0){
		i.binned=1
		n.bins=unique(Dat[,6])
	}
	n.beta.hab=ncol(DM.hab)
	if(is.null(Det.cov)==1)n.det.cov=0
	else n.det.cov=ncol(DM.det)
	N.remain=Adj*0  #keep track of "remaining" abundance in each habitat area for portion of area not covered by transects
	N.total=N.remain #keep track of total abundance in area of interest
	
	#Check to make sure input values are internally consistent
	#Later...
	
	#Initialize data augmentation multi-d array, parameter vectors and matrices
	Data<-array(0,dim=c(n.transects,M,4)) #array cols are Y1,Y2,Distance,Group size
	Obs.id=matrix(0,n.transects,2)  #holds observer id's for each transect
	cur.pl=1
	N.transect=rep(0,n.transects)
	for(itrans in 1:n.transects){
		Cur.dat=Dat[which(Dat[,1]==itrans),]
		Obs.id[itrans,]=Cur.dat[1,2:3]
		Data[itrans,cur.pl:(cur.pl+nrow(Cur.dat)-1),]=Cur.dat[,4:7]
		N.transect[itrans]=nrow(Cur.dat)		#initialize abundance in each transect area to be = to total number of animals observed
		cur.pl=cur.pl+nrow(Cur.dat)
	}
	N.obs=N.transect #vector holding total number of observations by transect
	
	cor.par=Inits$cor.par
	Beta.hab=data.frame(Inits$Beta.hab)
	names(Beta.hab)=colnames(DM.hab)
	Beta.det=data.frame(1)
	names(Beta.det)="Intercept"
	if(n.det.cov>0){
		Beta.det2=data.frame(Inits$Beta.det)
		names(Beta.det2)=colnames(Det.cov)
		Beta.det=cbind(Beta.det,Beta.det2)
	}
	Beta.det2=data.frame(matrix(0,1,n.observers-1))
	cur.str="Obs1"
	if(n.observers>2){
		for(iobs in 2:(n.observers-1))cur.str=c(cur.str,paste("Obs",iobs,sep=''))
	}
	names(Beta.det2)=cur.str
	Beta.det=cbind(Beta.det,Beta.det2)
	if(i.binned==0){ #continuous distances
		Beta.det2=data.frame(matrix(0,1,length(pol.eff)))
		names(Beta.det2)=paste(rep("Dist",length(pol.eff)),pol.eff,sep='')
	}
	else{
		Beta.det2=data.frame(matrix(0,1,n.bins-1))
		names(Beta.det2)=paste(rep("Dist",n.bins-1),c(2:n.bins),sep='')
	}
	Beta.det=cbind(Beta.det,Beta.det2)
	if(grps==TRUE){
		Beta.det2=data.frame(0)
		names(Beta.det2)="Grp"
		Beta.det=cbind(Beta.det,Beta.det2)
	}
	
	mcmc.length=Control$iter-Control$burnin
	MCMC=list(N=rep(0,mcmc.length),Beta.hab=data.frame(matrix(0,mcmc.length,length(Beta.hab))),Beta.det=data.frame(matrix(0,mcmc.length,length(Beta.det))),cor.par=rep(0,mcmc.length))
	colnames(MCMC$Beta.hab)=colnames(DM.hab)
	colnames(MCMC$Beta.det)=colnames(Beta.det)
	
	#tools for quicker calculations
	Cor=matrix(c(1,cor.par,cor.par,1),2,2)
	if(i.binned==1){
		X.cov.templ=rep(0,n.observers+n.bins)
		X.cov.templ[1]=1
	}
	
	for(iiter in 1:Control$iter){
		if((iiter%%1000)==1)cat(paste('\n ', iiter))

		Sample=c(-20:20)
		Sample=Sample[-21]
		max.a=max(Sample)
		Accept=list(N=0,sigma1=0,sigma2=0,dist=0)
		
		### update abundance parameters at the strata scale
		
		Lambda=Beta.hab*DM.hab  #keep these at expected values for the moment
		
		### update abundance, distances, ind. covariates within transect (RJMCMC)
		for(itrans in 1:n.transects){
			Sample=c(-Control$RJ.N:Control$RJ.N)
			Sample=Sample[-(Control$RJ.N+1)] #don't make proposals where pop size stays the same
			X.cov.templ1=X.cov.templ		#holds partial design vector
			X.cov.templ2=X.cov.templ
			if(Obs.id[itrans,1]<n.observers)X.cov.templ1[1+Obs.id[itrans,1]]=1
			if(Obs.id[itrans,2]<n.observers)X.cov.templ2[1+Obs.id[itrans,2]]=1
			

			#update distances & group sizes for latent animals
			if(N.transect[itrans]>N.obs[itrans]){
				for(iind in (N.obs[itrans]+1):N.transect[itrans]){
					if(i.binned==1){
						dist.star=sample(c(1:n.bins),1)
						grp.star=rpois(1,grp.mean)
						X1=X.cov.templ1
						if(dist.star>1){
							X1[(n.observers+1):(n.observers+dist.star-1)]=1
						}
						X1[n.observers+n.bins]=grp.star
						X2=X1
						X2[2:n.observers]=X.cov.templ2[2:n.observers]
						LogL.star=log(pnorm(upper=c(0,0),mean=c(X1%*%Beta.det,X2%*%Beta.det),sigma=Cor)
						X1=X.cov.templ1
						if(Data[itrans,iind,3]>1){
							X1[(n.observers+1):(n.observers+Data[itrans,iind,3]-1)]=1
						}
						X1[n.observers+n.bins]=grp.star
						X2=X1
						X2[2:n.observers]=X.cov.templ2[2:n.observers]						
						LogL.old=log(pnorm(upper=c(0,0),mean=c(X1%*%Beta.det,X2%*%Beta.det),sigma=Cor)
						if(runif(1)<exp(LogL.star-LogL.old)){
							Data[itrans,iind,3]=dist.star
						}					
					}
					else{
#						dist.star=Data[iind,3]+runif(1,-Control$MH.dist,Control$MH.dist)
#						if(dist.star>0 & dist.star<1){
#							LogL.star=log(1-p1*dnorm(dist.star,0,sigma1)/norm.const1)+log(1-p2*dnorm(dist.star,0,sigma2)/norm.const2)
#							LogL.old=log(1-p1*dnorm(Dat[iind,2],0,sigma1)/norm.const1)+log(1-p2*dnorm(Dat[iind,2],0,sigma2)/norm.const2)
#							if(runif(1)<exp(LogL.star-LogL.old)){
#								Dat[iind,2]=dist.star
#								Accept$dist=Accept$dist+1
#							}
#						}
					}
				}
			}

			Dat[(N+1):M,2]=runif(M-N)
			
			
		}
		
		### update detection parameters conditional on current abundances

	}

#}

	