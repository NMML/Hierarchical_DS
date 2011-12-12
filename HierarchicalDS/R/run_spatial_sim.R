#' function to simulate data and run a simple hierarchical example analysis using said dataset
#' @return see help for hierarchical_DS
#' @export
#' @keywords simulation
#' @author Paul B. Conn
run_spatial_sim=function(){
	S=100 #this needs to be a square number (square grid assumed)
	n.real.transects=sqrt(S)
	n.transects=n.real.transects*2 #each transect spans four cells
	Observers=matrix(NA,2,n.transects)
	Obs.cov=array(0,dim=c(2,n.transects,1))
	Obs.cov[1,,]=1
	n.obs.cov=dim(Obs.cov)[3]
	set.seed(93284)
	for(i in 1:n.real.transects){
		Observers[,((i-1)*2+1):((i-1)*2+2)]=rep(sample(c(1,2,3),size=2,replace=FALSE),2)
		#if(i%%2==1)Observers[,((i-1)*2+1):((i-1)*2+2)]=rep(sample(c(1,2,3),size=2,replace=FALSE),2)
		#else Observers[1,((i-1)*2+1):((i-1)*2+2)]=rep(sample(c(1,2,3),size=1,replace=FALSE),2)
	}
	Levels=list(Observer=c("1","2","3"),Seat=c("0","1"),Distance=c("1","2","3","4","5"),Species=c("1","2","3"))
	Adj=square_adj(n.real.transects)
	Out1=simulate_data_spatial(S=S,Observers=Observers,Adj=Adj,tau=1000) 
	Dat=Out1$Dat
	Mapping=Out1$Mapping
	Area.trans=Out1$Area.trans
	#Dat=Dat[,-8] uncomment if not modeling species effect
	n.bins=length(unique(Dat[,6]))
	Area.hab=rep(1,S)
	Bin.length=rep(1,n.bins)
	Hab.cov=data.frame(matrix(rep(log(c(1:sqrt(S))/sqrt(S)),sqrt(S)),S,1)) #covariate on abundance intensity
	colnames(Hab.cov)="Cov1"
	Hab.formula=~Cov1
	#Det.formula=~Observer+Seat+Distance+Group
	Det.formula=~Observer+Distance+Group+Species
	#Det.formula=~Distance+Group+Species
	Cov.prior.pdf=c("pois1","multinom")
	#Cov.prior.parms=matrix(c(1.1,1,1,0,1,1,1,1),4,2)
	#Cov.prior.parms=matrix(c(mean(Dat[,"Group"])-1,1,0,0,1,1,1,1),4,2) #prior mode near 2...
	Cov.prior.parms=matrix(c(3,1,0,1,1,1),3,2)
	#Cov.prior.parms=matrix(c(3,0,0,0,.1,.2,.3,.4),4,2)
	Cov.prior.fixed=c(0,0)
	#Cov.prior.pdf="pois1"
	#Cov.prior.parms=matrix(c(1.1,1,1),3,1)
	#Cov.prior.fixed=0
	#Cov.prior.pdf="pois1"
	#Cov.prior.parms=matrix(3,1,1)
	#Cov.prior.fixed=1
	colnames(Cov.prior.parms)=c("Group","Species")
	#colnames(Cov.prior.parms)="Group"
	pol.eff=2 #not currently used since using distance bins
	point.ind=TRUE
	spat.ind=TRUE #dont' include spatial dependence unless there really is spatial structure!
	grps=TRUE
	M=Out1$True.G*3
	M[which(M<10)]=15
	#M=c(11:20)
	#M=c(21:50)
	Control=list(iter=100100,burnin=100,thin=5,MH.cor=0.2,MH.nu=.01,MH.beta=c(.2,.4),RJ.N=rep(5,S),adapt=500)
	Inits=NULL
	Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.sd=c(10000,100)) #(1,.01) prior makes it closer to a uniform distribution near the origin
	adapt=TRUE
	
	set.seed(8327329)   #console 1
	Out=hierarchical_DS(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Observers=Observers,Bin.length=Bin.length,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.formula=Hab.formula,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,pol.eff=NULL,point.ind=TRUE,spat.ind=spat.ind,Inits=NULL,grps=grps,M=M,Control=Control,Levels=Levels,adapt=TRUE,Prior.pars=Prior.pars)

	#plot estimated abundance
	N.mean=apply(Out$MCMC$N,2,'mean')
	Abund.df=data.frame(cbind(rep(c(1:sqrt(S)),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(N.mean))))
	colnames(Abund.df)=c("y","x","Abundance")
	require(ggplot2)
	crap<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,750))+xlab("")+ylab("")
	crap
	
}