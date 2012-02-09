#' function to simulate distance sampling data (with misID) and run a simple hierarchical example analysis using said dataset
#' @return see help for hierarchical_DS
#' @export
#' @keywords simulation
#' @author Paul B. Conn
run_analysis=function(){
	S=25 #this needs to be a square number (square grid assumed)
	n.real.transects=S
	n.transects=S #each transect spans two cells
	Observers=matrix(NA,2,n.transects)
	Obs.cov=array(0,dim=c(2,n.transects,1))
	Obs.cov[1,,]=1
	n.obs.cov=1  #Seat occurs in dataset so must be included here
	set.seed(93284)
	#set.seed(932866)
	for(i in 1:n.real.transects){
		Observers[,i]=sample(c(1,2,3),size=2,replace=FALSE)
		#if(i%%2==1)Observers[,((i-1)*2+1):((i-1)*2+2)]=rep(sample(c(1,2,3),size=2,replace=FALSE),2)
		#else Observers[1,((i-1)*2+1):((i-1)*2+2)]=rep(sample(c(1,2,3),size=1,replace=FALSE),2)
	}
	Levels=list(Observer=c("1","2","3"),Seat=c("0","1"),Distance=c("1","2","3","4","5"),Species=c("1","2"))
	Adj=square_adj(sqrt(S))
	Out=simulate_data_misID(S=S,Observers=Observers)
	Dat=Out$Dat
	Mapping=c(1:S)
	Area.trans=rep(1,S)
	#Dat=Dat[,-8] uncomment if not modeling species effect
	n.bins=length(unique(Dat[,7]))
	Area.hab=rep(1,S)
	Bin.length=rep(1,n.bins)
	
	Hab.cov=data.frame(cbind(log(c(1:S)/S),(log(c(1:S)/S))^2)) #covariate on abundance intensity
	colnames(Hab.cov)=c("Cov1","Cov2")
	#Hab.formula=~Cov1
	Hab.formula=c(~Cov1+Cov2,~Cov1+Cov2)
	#Det.formula=~Observer+Seat+Distance+Group
	Det.formula=~Observer+Distance+Group+Species
	misID.mat=matrix(0,2,3)  # misID matrix 
	misID.mat[1,]=c(1,-1,2)  # positive numbers specify a model, 0 denotes impossible, -1 denotes obtain by subtraction
	misID.mat[2,]=c(-1,3,4)
	misID.models=c(~Observer+Distance+Group+Species,~Observer,~Observer+Distance+Group+Species,~Observer)
	n.species=nrow(misID.mat)
	#Det.formula=~Distance+Group+Species
	#Cov.prior.pdf=c("pois1_ln","multinom")
	#Cov.prior.parms=matrix(c(1.1,1,1,1,1,1),3,2)
	#Cov.prior.parms=matrix(c(mean(Dat[,"Group"])-1,1,0,0,1,1,1,1),4,2) #prior mode near 2...
	Cov.prior.parms=array(0,dim=c(n.species,3,1))
	Cov.prior.parms[1,,1]=c(1.1,1,.1)
	Cov.prior.parms[2,,1]=c(1.1,1,.1)
	#Cov.prior.parms=matrix(c(3,0,0,.5,.3,.2),3,2)
	Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
	Cov.prior.pdf=Cov.prior.fixed
	Cov.prior.pdf[,1]=c("pois1_ln","pois1_ln")
	#Cov.prior.parms=matrix(c(1.1,1,1),3,1)
	#Cov.prior.fixed=0
	#Cov.prior.pdf="pois1"
	#Cov.prior.parms=matrix(3,1,1)
	#Cov.prior.fixed=1
	#colnames(Cov.prior.parms)="Group"
	pol.eff=2 #not currently used since using distance bins
	point.ind=TRUE
	spat.ind=TRUE #dont' include spatial dependence unless there really is spatial structure!
	fix.tau.nu=FALSE
	srr=FALSE
	srr.tol=0.2
	grps=TRUE
	M=t(Out$G.true*3)
	M[which(M<10)]=15
	#M=c(11:20)
	#M=c(21:50)
	Control=list(iter=5100,burnin=100,thin=10,MH.cor=0.2,MH.nu=c(.01,.01),MH.beta=matrix(.2,2,3),RJ.N=matrix(rep(5,S*n.species),n.species,S),adapt=1000)
	#Control=list(iter=2010,burnin=10,thin=10,MH.cor=0.2,MH.nu=.01,MH.beta=c(.2,.4),RJ.N=rep(5,S),adapt=100)
	hab=matrix(0,n.species,3) #covariates are intercept, index, index^2
	hab[1,1:3]=c(log(40),1,0)
	hab[2,1:3]=c(log(10),-2,-1)
	Inits=list(hab=hab,tau.nu=c(100,100)) #chain1
	Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.sd=100) #(1,.01) prior makes it closer to a uniform distribution near the origin
	adapt=TRUE
	
	set.seed(8327329)   #chain1
	set.seed(8327330)   #chain2
	Out=hierarchical_DS_misID(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Observers=Observers,Bin.length=Bin.length,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.formula=Hab.formula,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,pol.eff=NULL,point.ind=TRUE,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,Inits=Inits,grps=grps,M=M,Control=Control,Levels=Levels,adapt=TRUE,Prior.pars=Prior.pars,misID.mat=misID.mat,misID.models=misID.models)
	plot_obs_pred(Out)
	summary_N(Out)
	
	save(Out,file="OutCluster2b.Rdata")
	
	#plot estimated abundance
	N.mean=apply(Out$MCMC$N,2,'mean')
	Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(N.mean))))
	colnames(Abund.df)=c("y","x","Abundance")
	require(ggplot2)
	crap<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
	crap
}