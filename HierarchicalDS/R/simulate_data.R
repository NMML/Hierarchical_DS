#' function to simulate distance sampling data from simple model with increasing abundance
#' intensity, no assumed spatial structure, and point independence.  If no parameters are given, uses internally defined values.
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @param Observers A (2 x S) matrix giving the observer identity for each transect
#' @param misID If TRUE (default), misidentification is assumed
#' @param X.site model.matrix for habitat covariates (defaults to a linear, quadratic effects of transect # on log scale)
#' @param n.species Number of species to simulate (current max is 2) (default is 2)
#' @param Beta.hab A (# of species X # parameters) matrix giving parameters for habitat-abundance relationship (default is linear increase for species 1, quadratic for species 2)
#' @param detect.model	A formula for the detection model.  Default is ~Observer+Distance+Group+Species (formula should consist of these key words)
#' @param Beta.det A vector giving parameters for the detection function; # of parameters must match model.matrix()!
#' @param dist.cont If TRUE, uses continuous distances on (0,1).  If FALSE, uses discrete distance classes
#' @param n.bins If dist.cont=FALSE, how many bins to use for distances.  Default is 5.
#' @param cor.par Correlation at maximum distance.  Default is 0.5
#' @param Grp.par	A vector with an entry for each species giving parameters for group size (assumed zero-truncated Poisson). Default is 3 and 1, corresponding to mean group sizes of 4 and 2 for each species
#' @param Conf.model If misID=TRUE, this gives a list vector with formulas for the confusion matrix. Formulas are provided for "getting it right" and "getting an unknown".  Default is c(~Observer+Group+Distance+Species, ~Observer). The same models are assumed for each species.
#' @param Conf.par	If misID=TRUE, ths is a list vector with parameters for each confusion model.  Default is list(c(2,.2,.4,.3,-.1,-.2,-.4,-.8,.5),c(-1,-.2,.2)) 

#' @return a distance sampling dataset
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_data<-function(S,Observers,misID=TRUE,X.site=NULL,n.species=2,Beta.hab=NULL,Beta.det=NULL,detect.model=~Observer+Distance+Group+Species,dist.cont=FALSE,n.bins=5,cor.par=0.5,Grp.par=c(3,1),Conf.model=NULL,Conf.par=NULL){
	require(mvtnorm)
	
	if(n.species>2)cat("\n Error: current max species is 2 \n")
	if(n.species==1 & misID==TRUE){
		cat("\n n.speces=1 so misID set to FALSE \n")
		misID=FALSE
	}
	if((n.bins!=5 | dist.cont==TRUE) & (is.null(Beta.det)==TRUE))cat("\n Error: if using continuous distances or a non-default distance bin #, you must input Beta.det \n")
	#process parameters
	if(is.null(X.site)==TRUE)X.site=cbind(rep(1,S),log(c(1:S)/S),(log(c(1:S)/S))^2) #covariate on abundance intensity, sp 1
	if(is.null(Beta.hab)==TRUE){
		Beta.hab=matrix(0,n.species,3)
		Beta.hab[1,]=c(log(40),1,0) 
		if(n.species==2)Beta.hab[2,]=c(log(10),-2,-1)
	}
	
	#detection parameters
	if(dist.cont==FALSE)Levels=list(Observer=sort(unique(c(Observers))),Distance=as.factor(c(1:n.bins)),Species=as.factor(1:n.species))
	else Levels=list(Observer=unique(c(Observers)),Species=as.factor(1:n.species))
	factor.ind=list(Observer=TRUE,Distance=(dist.cont==FALSE),Group=FALSE,Species=TRUE)
	
	if(is.null(Beta.det)==TRUE)Beta.det=c(1.2,-.2,-.4,-.6,-.9,-1.1,-1.3,.1,.3)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
												
	N1=round(exp(X.site%*%Beta.hab[1,]))
	N2=N1*0
	if(n.species==2)N2=round(exp(X.site%*%Beta.hab[2,]))
	cat(paste("\n True G, sp 1 = ",N1,'\n\n'))
	cat(paste("\n True G.tot, sp 1= ",sum(N1),'\n'))
	if(n.species==2){
		cat(paste("\n True G, sp 2 = ",N2,'\n\n'))
		cat(paste("\n True G.tot, sp 2= ",sum(N2),'\n'))
	}
	
	Dat=matrix(0,sum(N1)+sum(N2),8)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size
	
	#initialize confusion matrix, etc.
	if(misID==1){
		Confusion=matrix(0,2,3)
		Confusion[1,]=c(1,-1,2)
		Confusion[2,]=c(-1,1,2)
		if(is.null(Conf.model)==TRUE)Conf.model=c(~Observer+Group+Distance+Species,~Observer)
		if(is.null(Conf.par)==TRUE){
			Conf.par=list(c(2,.2,.4,.3,-.1,-.2,-.4,-.8,.5),c(-1,-.2,.2)) #parameters for getting an 'unknown'
		}
	}
	
	pl=1
	for(i in 1:S){
		cur.Observers=Observers[,i]
		if(N1[i]>0){
			for(j in 1:N1[i]){
				if(dist.cont==TRUE)dist=runif(1)
				else dist=sample(c(1:n.bins),1)
				grp.size=rpois(1,Grp.par[1])+1
				Dat1=matrix(c(cur.Observers[1],dist,grp.size,1),1,4)
				Dat2=Dat1
				Dat2[1]=cur.Observers[2]
				X1=get_mod_matrix(Cur.dat=Dat1,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				X2=get_mod_matrix(Cur.dat=Dat2,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=dist
				Dat[pl,7]=grp.size
				Dat[pl,8]=1
				mu1=X1%*%Beta.det
				mu2=X2%*%Beta.det
				if(dist.cont==FALSE)cur.cor=(dist-1)/(n.bins-1)*cor.par
				else cur.cor=dist*cor.par
				Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		if(N2[i]>0){
			for(j in 1:N2[i]){
				if(dist.cont==TRUE)dist=runif(1)
				else dist=sample(c(1:n.bins),1)
				grp.size=rpois(1,Grp.par[2])+1
				Dat1=matrix(c(cur.Observers[1],dist,grp.size,2),1,4)
				Dat2=Dat1
				Dat2[1]=cur.Observers[2]
				X1=get_mod_matrix(Cur.dat=Dat1,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				X2=get_mod_matrix(Cur.dat=Dat2,stacked.names=c("Observer","Distance","Group","Species"),factor.ind=factor.ind,Det.formula=detect.model,Levels=Levels)
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=dist
				Dat[pl,7]=grp.size
				Dat[pl,8]=2
				mu1=X1%*%Beta.det
				mu2=X2%*%Beta.det
				if(dist.cont==FALSE)cur.cor=(dist-1)/(n.bins-1)*cor.par
				else cur.cor=dist*cor.par
				Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		
	}
	cat(paste("Total N, sp 1 = ",sum(Dat[Dat[,8]==1,7])))
	if(n.species==2)cat(paste("Total N, sp 2 = ",sum(Dat[Dat[,8]==2,7])))
	Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),] #get rid of animals never observed
	
	#put things in "Jay's" format
	Dat2=matrix(0,2*nrow(Dat),7)
	ipl=1
	for(irecord in 1:nrow(Dat)){
		Dat2[ipl,1]=Dat[irecord,1]
		Dat2[ipl+1,1]=Dat[irecord,1]
		Dat2[ipl,3]=Dat[irecord,2]
		Dat2[ipl+1,3]=Dat[irecord,3]
		Dat2[ipl,4]=Dat[irecord,4]
		Dat2[ipl+1,4]=Dat[irecord,5]
		Dat2[ipl,5]=Dat[irecord,6]
		Dat2[ipl+1,5]=Dat[irecord,6]
		Dat2[ipl,6]=Dat[irecord,7]
		Dat2[ipl+1,6]=Dat[irecord,7]
		Dat2[ipl,7]=Dat[irecord,8]
		Dat2[ipl+1,7]=Dat[irecord,8]
		Dat2[ipl,2]=irecord  #match number
		Dat2[ipl+1,2]=irecord
		ipl=ipl+2
	}
	Dat2=as.data.frame(Dat2)
	colnames(Dat2)=c("Transect","Match","Observer","Obs","Distance","Group","Species")
	Dat2[,"Observer"]=as.factor(Dat2[,"Observer"])
	Dat2[,"Distance"]=as.factor(Dat2[,"Distance"])
	
	Dat=Dat2
	True.species=Dat[,"Species"]
	
	if(misID==TRUE){
		# Now, put in partial observation process
		Ind.sp1=which(Dat[,"Species"]==1)
		Ind.sp2=which(Dat[,"Species"]==2)
		Dat1=Dat[Ind.sp1,]
		Dat1[,"Species"]=0  # i think species is continuous here
		X=model.matrix(Conf.model[[1]],data=Dat1)
		MN.logit.right=exp(X%*%Conf.par[[1]])
		X=model.matrix(Conf.model[[2]],data=Dat1)
		MN.logit.unk=exp(X%*%Conf.par[[2]])
		Psi=matrix(0,nrow(Dat1),3)
		Psi[,1]=MN.logit.right/(1+MN.logit.right+MN.logit.unk)
		Psi[,3]=MN.logit.unk/(1+MN.logit.right+MN.logit.unk)
		Psi[,2]=1-Psi[,1]-Psi[,3]
		get_samp<-function(prob)sample(c(1:length(prob)),1,prob=prob)
		Cur.sp=apply(Psi,1,'get_samp')	
		Dat[Ind.sp1,"Species"]=Cur.sp
		
		Dat1=Dat[Ind.sp2,]
		Dat1[,"Species"]=1
		X=model.matrix(Conf.model[[1]],data=Dat1)
		MN.logit.right=exp(X%*%Conf.par[[1]])
		X=model.matrix(Conf.model[[2]],data=Dat1)
		MN.logit.unk=exp(X%*%Conf.par[[2]])
		Psi=matrix(0,nrow(Dat1),3)
		Psi[,2]=MN.logit.right/(1+MN.logit.right+MN.logit.unk)
		Psi[,3]=MN.logit.unk/(1+MN.logit.right+MN.logit.unk)
		Psi[,1]=1-Psi[,2]-Psi[,3]
		get_samp<-function(prob)sample(c(1:length(prob)),1,prob=prob)
		Cur.sp=apply(Psi,1,'get_samp')	
		Dat[Ind.sp2,"Species"]=Cur.sp
	}
	Dat[,"Species"]=as.integer(Dat[,"Species"])
	
	Dat=cbind(Dat[,1:4],Dat[,"Species"],Dat[,5:6])
	colnames(Dat)[5]="Species"
	G.true=N1
    if(n.species==2)G.true=cbind(G.true,N2)
	
	Out=list(Dat=Dat,G.true=G.true,True.species=True.species)
}

