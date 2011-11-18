#' function to simulate data and run a simple hierarchical example analysis using said dataset
#' @return see help for hierarchical_DS
#' @export
#' @keywords simulation
#' @author Paul B. Conn
run_analysis=function(){
	S=10
	Observers=matrix(0,2,S)
	Obs.cov=array(0,dim=c(2,S,1))
	Obs.cov[1,,]=1
	n.obs.cov=dim(Obs.cov)[3]
	for(i in 1:S){
		Observers[,i]=sample(c(1,2,3),size=2,replace=FALSE)
	}
	Levels=list(Observer=c("1","2","3"),Seat=c("0","1"),Distance=c("1","2","3","4","5"),Species=c("1","2","3","4"))
	Dat=simulate_data(S=S,Observers=Observers) 
	#Dat=Dat[,-8] uncomment if not modeling species effect
	n.bins=length(unique(Dat[,6]))
	Adj=matrix(0,S,S)
	for(i in 2:S){
		Adj[i-1,i]=1
		Adj[i,i-1]=1
	}
	Area.hab=rep(1,S)
	Mapping=c(1:S)
	Area.trans=rep(0.5,S)
	Bin.length=rep(1,n.bins)
	Hab.cov=data.frame(matrix(log(c(1:S)/S),S,1))
	colnames(Hab.cov)="Cov1"
	Hab.formula=~Cov1
	#Det.formula=~Observer+Seat+Distance+Group
	Det.formula=~Observer+Seat+Distance+Group+Species
	Cov.prior.pdf=c("pois1","multinom")
	Cov.prior.parms=matrix(c(3,1,1,0,1,1,1,1),4,2)
	Cov.prior.fixed=c(0,0)
	#Cov.prior.pdf="pois1_ln"
	#Cov.prior.parms=matrix(c(1.1,1,1),3,1)
	#Cov.prior.fixed=0
	#Cov.prior.pdf="pois1"
	#Cov.prior.parms=matrix(3,1,1)
	#Cov.prior.fixed=1
	colnames(Cov.prior.parms)=c("Group","Species")
	#colnames(Cov.prior.parms)="Group"
	pol.eff=2 #not currently used since using distance bins
	point.ind=TRUE
	spat.ind=TRUE  #dont' include spatial dependence unless there really is spatial structure!
	grps=TRUE
	M=c(30,50,60,80,90,100,110,120,130,140)
	#M=c(11:20)
	#M=c(21:50)
	Control=list(iter=11100,burnin=100,thin=5,MH.cor=0.2,MH.nu=.01,MH.beta=c(.2,.4),RJ.N=rep(5,S),adapt=500)
	Inits=NULL
	Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.sd=c(10000,100)) #(1,.01) prior makes it closer to a uniform distribution near the origin
	adapt=TRUE
	
	Out=hierarchical_DS(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Observers=Observers,Bin.length=Bin.length,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.formula=Hab.formula,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,pol.eff=NULL,point.ind=TRUE,spat.ind=spat.ind,Inits=NULL,grps=grps,M=M,Control=Control,Levels=Levels,adapt=TRUE,Prior.pars=Prior.pars)
}