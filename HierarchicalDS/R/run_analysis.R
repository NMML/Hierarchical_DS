#' function to simulate data and run a simple hierarchical example analysis using said dataset
#' @param S	number of sites
#' @return see help for hierarchical_DS
#' @export
#' @keywords simulation
#' @author Paul B. Conn
run_analysis=function(S){
	Dat=simulate_data(S=S) 
	n.bins=length(unique(Dat[,6]))
	Adj=matrix(0,S,S)
	for(i in 2:S){
		Adj[i-1,i]=1
		Adj[i,i-1]=1
	}
	Area.hab=rep(1,S)
	Mapping=c(1:S)
	Area.trans=rep(1,S)
	Bin.length=rep(1,n.bins)
	Hab.cov=data.frame(matrix(log(c(1:S)/S),S,1))
	colnames(Hab.cov)="Cov1"
	n.obs.cov=1  #dummy 'seat' variable
	Hab.formula=~Cov1
	Det.formula=~Observer+Seat+Distance+Group
	Cov.prior.pdf="pois1"
	Cov.prior.parms=matrix(c(3,NA),2,1)
	colnames(Cov.prior.parms)="Group"
	pol.eff=2 #not currently used since using distance bins
	point.ind=TRUE
	spat.ind=FALSE
	grps=TRUE
	M=200
	Control=list(iter=15000,burnin=100,thin=5,MH.cor=0.05,MH.nu=.01,MH.beta=c(.2,.4),grp.mean=3,RJ.N=rep(5,S),adapt=1000)
	Inits=NULL
	Prior.pars=list(a.eta=.01,b.eta=.01,a.nu=.01,b.nu=.01,beta.sd=c(10000,100))
	adapt=TRUE
	
	Out=hierarchical_DS(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Bin.length=Bin.length,Hab.cov=Hab.cov,n.obs.cov=n.obs.cov,pol.eff=NULL,point.ind=TRUE,spat.ind=FALSE,Inits=NULL,grps=grps,M=M,Control=Control,adapt=TRUE,Prior.pars=Prior.pars)
}