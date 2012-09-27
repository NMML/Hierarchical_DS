#' script to simulate distance sampling data (with misID) and run a simple hierarchical example analysis using said dataset
#' #@return see help for hierarchical_DS.R
#' #@export
#' #@keywords simulation
#' #@author Paul B. Conn
#' 
  #1) First, simulate data
	S=16 #Number of grid cells; this needs to be a square number (square grid assumed)
  library(hierarchicalDS)
	n.transects=S #one transect per cell
	Observers=matrix(NA,2,n.transects)
	set.seed(11111) 
	for(i in 1:n.transects){
		Observers[,i]=sample(c(1,2,3),size=2,replace=FALSE)
	}
	Out=simulate_data(S=S,Observers=Observers,misID=TRUE,tau=10)
	Dat=Out$Dat
  #you can also look at the true number of groups per transect with Out$G.true; species 1 increases linearly across the x axis and is more abundant, while species 2 is less common and has a uniform distribution
	
  #2) declare inputs and call hierarchicalDS; note that this run takes about 3.5 hours.  Alternatively, load "sim_data" and proceed to plot/summarize results
	Obs.cov=array(0,dim=c(2,n.transects,1))
	Obs.cov[1,,]=1
	n.obs.cov=1 #1 observer covariate, "Seat" is included in the dataset even though it isn't modeled
	Adj=square_adj(sqrt(S))
	misID.mat=matrix(0,2,3)  # misID matrix 
	misID.mat[1,]=c(1,-1,2)  # positive numbers specify a model, 0 denotes impossible, -1 denotes obtain by subtraction
	misID.mat[2,]=c(-1,3,-1)
	misID.models=c(~1,~1,~1)
	MisID=vector("list",max(misID.mat))
	MisID[[1]]=2 #parameters for getting it right
	MisID[[2]]=1
	MisID[[3]]=3 #parameters for getting it right
	misID.symm=TRUE
	Mapping=c(1:S)
	Area.trans=rep(1,S)
	#Dat=Dat[,-8] uncomment if not modeling species effect
	n.bins=length(unique(Dat[,"Distance"]))
	Area.hab=rep(1,S)
	Bin.length=rep(1,n.bins)  #equal bin lengths	
	Hab.cov=data.frame(rep(log(c(1:sqrt(S)/sqrt(S))),each=sqrt(S))) #covariate on abundance intensity same as used to generate data
	colnames(Hab.cov)=c("Cov1")
	Hab.formula=c(~Cov1,~Cov1) #the density of species two isn't affected by the habitat covariate but we'll estimate this effect anyway
	detect=TRUE
	Det.formula=~Observer+Distance+Group
	n.species=nrow(misID.mat)
	Cov.prior.parms=array(0,dim=c(n.species,2,1))
	Cov.prior.parms[1,,1]=c(2,1)  #we'll put priors a little off from their true values; expected group sizes are 4 and 2 for each species
	Cov.prior.parms[2,,1]=c(2,1)
	Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
	Cov.prior.pdf=Cov.prior.fixed
	Cov.prior.pdf[,1]=c("pois1","pois1")  #model group size as a zero truncated poisson
	Cov.prior.n=matrix(2,2,1)
	point.ind=TRUE #include point independence
	spat.ind=FALSE #do not make spatially independent; i.e. estimate spatial autocorrelation!
  fix.tau.nu=FALSE 
	srr=FALSE
	srr.tol=0.2
	misID=TRUE
	grps=TRUE
  post.loss=TRUE   
	M=t(Out$G.true*10)
	M[which(M<30)]=50
  Control=list(iter=5500,burnin=500,thin=10,MH.cor=0.2,MH.nu=matrix(.2,2,S),MH.beta=c(.2,.4),MH.misID=matrix(0.1,3,1),RJ.N=matrix(rep(5,S*n.species),n.species,S),n.species,adapt=100,iter.fix.N=100)
  hab=matrix(0,n.species,2) #covariates are intercept, index
	hab[1,1:2]=c(log(50),0)
	hab[2,1:2]=c(log(10),-2)
	Inits=list(hab=hab,tau.nu=c(500,500),MisID=MisID) #provide some initial values to ensure MCMC doesn't start out at weird place
	misID.mu=vector("list",max(misID.mat))
	misID.sd=misID.mu
	misID.mu[[1]]=0
	misID.mu[[2]]=0
	misID.mu[[3]]=0
	misID.sd[[1]]=1.75
	misID.sd[[2]]=1.75
	misID.sd[[3]]=1.75
	Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.tau=0.01,misID.mu=misID.mu,misID.sd=misID.sd) #(1,.01) prior makes it closer to a uniform distribution near the origin
	adapt=TRUE	
	set.seed(8327329)   #chain1
	Out=hierarchical_DS(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Observers=Observers,Bin.length=Bin.length,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.formula=Hab.formula,detect=detect,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,pol.eff=NULL,point.ind=point.ind,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,misID=misID,Inits=Inits,grps=grps,M=M,Control=Control,adapt=adapt,Prior.pars=Prior.pars,misID.mat=misID.mat,misID.models=misID.models,misID.symm=misID.symm)

  ###3) plot and summarize results; note that chain would need to be run a lot longer to summarize the posterior very well!!!
  plot(Out$MCMC)
	summary_N(Out)
  post_loss(Out)
	

