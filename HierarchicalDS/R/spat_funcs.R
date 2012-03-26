#' function to sample from a specified probability density function
#' @param n number of samples desired
#' @param pdf probability density function (pois1, poisson, normal, unif.disc, unif.cont)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @param RE random effects, if present
#' @return a vector of length n samples from the desired distribution 
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_sample<-function(n,pdf,cur.par,RE){
	switch(pdf,
			pois1=rpois(n,cur.par[1])+1,
			poisson=rpois(n,cur.par[1]),
			pois1_ln=rpois(n,exp(cur.par[1]+cur.par[2]*RE))+1,
			poisson_ln=rpois(n,exp(cur.par[1]+cur.par[2]*RE)),
			normal=rnorm(n,cur.par[1],cur.par[2]),
			unif.disc=sample(cur.par[1]:cur.par[2],n,replace=TRUE),
			unif.cont=runif(n,cur.par[1],cur.par[2]),
			multinom=sample(c(1:length(cur.par)),n,replace=TRUE,prob=cur.par)
	)
}

#' function to sample from hyperpriors of a specified probability density function; note that
#' initial values for sigma of lognormal random effects are fixed to a small value (0.05) to
#' prevent numerical errors
#' @param pdf probability density function (pois1, poisson, normal, unif.disc, unif.cont)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @return a vector of length n samples from the desired distribution 
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_sample_prior<-function(pdf,cur.par){
	require(mc2d)
	switch(pdf,
			pois1=rgamma(1,cur.par[1],cur.par[2]),
			poisson=rgamma(1,cur.par[1],cur.par[2]),
			pois1_ln=c(rnorm(1,cur.par[1],cur.par[2]),0.05),
			poisson_ln=c(rnorm(1,cur.par[1],cur.par[2]),0.05),
			multinom=rdirichlet(1,cur.par)
	)
}

#' function to calculate the joint pdf for a sample of values from one of a number of pdfs
#' @param x values to be evaluated
#' @param pdf probability density function (pois1, poisson, pois1_ln, poisson_ln, normal, multinom)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @param RE random effects, if present
#' @return total log likelihood of points
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_pdf<-function(x,pdf,cur.par,RE){
	switch(pdf,
			pois1=sum(dpois(x-1,cur.par[1],log=1)),
			poisson=sum(dpois(x,cur.par[1],log=1)),
			pois1_ln=sum(dpois(x-1,exp(cur.par[1]+cur.par[2]*RE),log=1)),
			poisson_ln=sum(dpois(x,exp(cur.par[1]+cur.par[2]*RE),log=1)),
			normal=sum(dnorm(x,cur.par[1],cur.par[2],log=1)),
			multinom=sum(log(cur.par[x]))
	)
}

#' function to stack data (going from three dimensional array to a two dimensional array including only "existing" animals
#' @param Data three-d dataset
#' @param Obs.transect	current number of observations of animals in each transect (vector)
#' @param n.transects	number of transects
#' @param stacked.names	column names for new stacked dataset
#' @param factor.ind	a vector of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
#' @return a stacked dataset
#' @export
#' @keywords stack data
#' @author Paul B. Conn
stack_data<-function(Data,Obs.transect,n.transects,stacked.names,factor.ind){
	#convert from "sparse" 3-d data augmentation array to a rich 2-d dataframe for updating beta parameters 
	if(n.transects==1)Stacked=Data
	else{
		Stacked=as.data.frame(Data[1,1:2,])
		for(itrans in 1:n.transects){
			if(Obs.transect[itrans]>0)Stacked=rbind(Stacked,Data[itrans,1:Obs.transect[itrans],])
		}
		Stacked=Stacked[-c(1,2),]
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

#' function to produce a design matrix given a dataset and user-specified formula object
#' @param Cur.dat 	current dataset
#' @param stacked.names	column names for current dataset
#' @param factor.ind	a list of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
#' @param Det.formula	a formula object
#' @param Levels	A list object giving the number of levels for factor variables
#' @return a design matrix
#' @export
#' @keywords model matrix
#' @author Paul B. Conn
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
	DM
}

#' generate initial values for MCMC chain if not already specified by user
#' @param DM.hab 	design matrix for habitat model
#' @param DM.det	design matrix for detection model
#' @param G.transect a vector of the number of groups of animals in area covered by each transect		
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param point.ind	is point independence assumed (TRUE/FALSE)
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  pois1 parameter for group size
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits<-function(DM.hab,DM.det,G.transect,Area.trans,Area.hab,Mapping,point.ind,spat.ind,grp.mean){		
	Par=list(det=rnorm(ncol(DM.det),0,1),hab=rep(0,ncol(DM.hab)),cor=ifelse(point.ind,runif(1,0,.8),0),
			Nu=log(max(G.transect)/mean(Area.trans)*exp(rnorm(length(Area.hab)))),Eta=rnorm(length(Area.hab)),
			tau.eta=runif(1,0.5,2),tau.nu=runif(1,0.5,2))
	Par$hab[1]=mean(G.transect)/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(1,0,1))
	Par$G=exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu)))
	Par$N=Par$G+rpois(length(Par$G),grp.mean*Par$G)
	if(spat.ind==1)Par$Eta=0*Par$Eta
	Par
}

#' generate initial values for misID model if not already specified by user
#' @param DM.hab 	a list of design matrices for the habitat model (elements are named sp1,sp2, etc.)
#' @param DM.det	design matrix for detection model
#' @param N.hab.par  vector giving number of parameters in the habitat model for each species
#' @param G.transect a matrix of the number of groups of animals in area covered by each transect; each row gives a separate species		
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param point.ind	is point independence assumed (TRUE/FALSE)
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  a vector giving the pois1 parameter for group size (one entry for each species)
#' @param misID    if TRUE, indicates that misidentification is incorporated into modeling
#' @param misID.mat a matrix specifying which elements of the misID matrix are linked to model equations
#' @param N.par.misID a vector giving the number of parameters for each misID model (in multinomial logit space)
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits_misID<-function(DM.hab,DM.det,N.hab.par,G.transect,Area.trans,Area.hab,Mapping,point.ind,spat.ind,grp.mean,misID,misID.mat,N.par.misID){		
	n.species=nrow(G.transect)
	n.cells=length(Area.hab)
	MisID=NULL
	if(misID){
		n.misID.eq=max(misID.mat)
		ncol.misID=max(N.par.misID)
		MisID=matrix(runif(n.misID.eq*ncol.misID,-.5,.5),n.misID.eq,ncol.misID)
		diag.mods=diag(misID.mat)
		diag.mods=diag.mods[which(diag.mods>0)]
		if(length(diag.mods)>0)MisID[diag.mods,1]=MisID[diag.mods,1]+2 #ensure that the highest probability is for a non-misID
	}
	hab=matrix(0,n.species,max(N.hab.par))
	Nu=matrix(0,n.species,n.cells)
	for(isp in 1:n.species){
		Nu[isp,]=log(max(G.transect[isp,])/mean(Area.trans)*exp(rnorm(length(Area.hab),0,0.1)))
	}
	Par=list(det=rnorm(ncol(DM.det),0,1),hab=hab,cor=ifelse(point.ind,runif(1,0,.8),0),
			Nu=Nu,Eta=matrix(rnorm(n.species*n.cells),n.species,n.cells),
			tau.eta=runif(n.species,0.5,2),tau.nu=runif(n.species,0.5,2),MisID=MisID)
	Par$hab[,1]=log(apply(G.transect,1,'mean')/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(n.species,0,1)))
	Par$G=exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu)))
	for(isp in 1:n.species)Par$N[isp,]=Par$G[isp,]+rpois(n.cells,grp.mean[isp]*Par$G[isp,])
	if(spat.ind==1)Par$Eta=0*Par$Eta
	Par
}


#' compute the first derivative of log_lambda likelihood component for Langevin-Hastings
#' @param Mu	expected value for all cells
#' @param Nu	current observed valus (all cells)
#' @param Sampled Vector giving the cell identities for all sampled cells
#' @param Area Proportional area of each sampled cell that is covered by one or more transects
#' @param N	    number of groups in each transect
#' @param var.nu	variance of the overdispersion process
#' @return a gradient value
#' @export
#' @keywords gradient, Langevin-Hastings
#' @author Paul B. Conn
log_lambda_gradient<-function(Mu,Nu,Sampled,Area,N,var.nu){
	Grad=(Mu[Sampled]-Nu[Sampled])/var.nu+N-Area*exp(Mu[Sampled])
	Grad 
}

#' compute the likelihood for nu parameters
#' @param Log.lambda 	Log of poisson intensities for total areas sampled in each sampled strata
#' @param DM	the design matrix
#' @param Beta	linear predictor parameters for the log of abundance intensity	
#' @param Eta	a vector of spatial random effects
#' @param SD	standard deviation of the overdispersion process
#' @param N		a vector giving the current iteration's number of groups in the area 
#' @param Sampled Index for which cells were actually sampled
#' @param Area	Total area sampled in each sampled cell
#' @return the log likelihood associated with the data and the current set of parameters 
#' @export
#' @keywords log likelihood
#' @author Paul B. Conn
log_lambda_log_likelihood<-function(Log.lambda,DM,Beta,Eta=0,SD,N,Sampled,Area){
	Pred.log.lam=(DM%*%Beta+Eta)[Sampled]	
	logL=sum(dnorm(Log.lambda,Pred.log.lam,SD,log=1)) #normal component
	logL=logL+sum(N*(Log.lambda+log(Area))-Area*exp(Log.lambda))
	return(logL)
} 	

#' SIMULATE AN ICAR PROCESS 
#' @param Q Precision matrix for the ICAR process
#' @return Spatial random effects
#' @export 
#' @keywords ICAR, simulation
#' @author Devin Johnson
rrw <- function(Q){
	v <- eigen(Q, TRUE)
	val.inv <- sqrt(ifelse(v$values>sqrt(.Machine$double.eps), 1/v$values, 0))
	P <- v$vectors
	sim <- P%*%diag(val.inv)%*%rnorm(dim(Q)[1], 0, 1)
	X <- rep(1,length(sim))
	if(sum(val.inv==0)==2) X <- cbind(X, 1:length(sim))
	sim <- sim-X%*%solve(crossprod(X), crossprod(X,sim))
	return(sim)
}

#' Produce an adjacency matrix for a vector
#' @param x length of vector
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn
linear_adj <- function(x){
	Adj1=matrix(0,x,x)
	Adj2=matrix(0,x,x)
	diag.min.1=diag(x-1)
	Adj1[2:x,1:(x-1)]=diag.min.1
	Adj2[1:(x-1),2:x]=diag.min.1
	Adj=Adj1+Adj2
	Adj
}	
	
#' Produce an adjacency matrix for a square grid
#' @param x number of cells on side of grid
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn
square_adj <- function(x){
	Ind=matrix(c(1:x^2),x,x)
	Adj=matrix(0,x^2,x^2)
	for(i in 1:x){
		for(j in 1:x){
			if(i==1 & j==1){
				Adj[Ind[i,j],Ind[i,j]+1]=1
				Adj[Ind[i,j],Ind[i,j]+x]=1
				Adj[Ind[i,j],Ind[i,j]+x+1]=1
			}
			if(i==1 & j>1 & j<x){
				Adj[Ind[i,j],Ind[i,j]+1]=1
				Adj[Ind[i,j],Ind[i,j]+x]=1
				Adj[Ind[i,j],Ind[i,j]-x]=1
				Adj[Ind[i,j],Ind[i,j]+x+1]=1
				Adj[Ind[i,j],Ind[i,j]-x+1]=1
			}
			if(i==1 & j==x){
				Adj[Ind[i,j],Ind[i,j]+1]=1
				Adj[Ind[i,j],Ind[i,j]-x]=1	
				Adj[Ind[i,j],Ind[i,j]-x+1]=1
			}
			if(i>1 & i<x & j==1){
				Adj[Ind[i,j],Ind[i,j]+1]=1
				Adj[Ind[i,j],Ind[i,j]+x]=1
				Adj[Ind[i,j],Ind[i,j]-1]=1
				Adj[Ind[i,j],Ind[i,j]+x-1]=1
				Adj[Ind[i,j],Ind[i,j]+x+1]=1
			}
			if(i>1 & i<x & j>1 & j<x){
				cur.nums=c(Ind[i,j]-x-1,Ind[i,j]-x,Ind[i,j]-x+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+x-1,Ind[i,j]+x,Ind[i,j]+x+1)
				Adj[Ind[i,j],cur.nums]=1
			}
			if(i>1 & i<x & j==x){
				Adj[Ind[i,j],Ind[i,j]+1]=1
				Adj[Ind[i,j],Ind[i,j]-x]=1
				Adj[Ind[i,j],Ind[i,j]-1]=1	
				Adj[Ind[i,j],Ind[i,j]-x-1]=1
				Adj[Ind[i,j],Ind[i,j]-x+1]=1
				
			}
			if(i==x & j==1){
				Adj[Ind[i,j],Ind[i,j]+x]=1
				Adj[Ind[i,j],Ind[i,j]-1]=1	
				Adj[Ind[i,j],Ind[i,j]+x-1]=1
			}
			if(i==x & j>1 & j<x){
				Adj[Ind[i,j],Ind[i,j]+x]=1
				Adj[Ind[i,j],Ind[i,j]-1]=1								
				Adj[Ind[i,j],Ind[i,j]-x]=1
				Adj[Ind[i,j],Ind[i,j]+x-1]=1
				Adj[Ind[i,j],Ind[i,j]-x-1]=1
			}
			if(i==x & j==x){
				Adj[Ind[i,j],Ind[i,j]-1]=1								
				Adj[Ind[i,j],Ind[i,j]-x]=1
				Adj[Ind[i,j],Ind[i,j]-x-1]=1
			}				
		}
	}
	return(Adj)
}

#' estimate optimal 'a' parameter for linex loss function
#' @param Pred.G  Predicted group abundance
#' @param Obs.G	Observed group abundance
#' @param min.a Minimum value for linex 'a' parameter
#' @param max.a Maximum value for linex 'a' parameter
#' @return The optimal tuning parameter for linex loss function as determined by minimum sum of squares 
#' @export
#' @keywords linex
#' @author Paul B. Conn
calc_linex_a<-function(Pred.G,Obs.G,min.a=0.00001,max.a=1.0){
	Y=apply(Obs.G,2,mean)
	linex_ssq<-function(a,X,Y){
		Theta=exp(-a*X)
		Theta=-1/a*log(apply(Theta,2,'mean'))
		return(sum((Y-Theta)^2))
	}
	a=optimize(f=linex_ssq,interval=c(min.a,max.a),X=Pred.G,Y=Y)
	a
} 	

#' plot 'observed' versus predicted values for abundance of each species at each transect
#' @param Out  Output list from "mcmc_ds.R" 	
#' @return NULL 
#' @export
#' @keywords diagnostics, plot
#' @author Paul B. Conn
plot_obs_pred<-function(Out){
	n.species=dim(Out$Pred.N)[1]
	par(mfrow=c(n.species,1))
	for(isp in 1:n.species){
		a.linex=calc_linex_a(Out$Pred.N[isp,,],Out$Obs.N[isp,,])$minimum
		max.x=max(c(apply(Out$Obs.N[isp,,],2,'mean'),apply(Out$Pred.N[isp,,],2,'mean')))
		plot(apply(Out$Obs.N[isp,,],2,'mean'),apply(Out$Pred.N[isp,,],2,'mean'),pch=1,xlim=c(0,max.x),ylim=c(0,max.x),xlab="Observed",ylab="Predicted")
		points(apply(Out$Obs.N[isp,,],2,'mean'),apply(Out$Pred.N[isp,,],2,'median'),pch=2)
		Theta=exp(-a.linex*Out$Pred.N[isp,,])
		Theta=-1/a.linex*log(apply(Theta,2,'mean'))
		points(apply(Out$Obs.N[isp,,],2,'mean'),Theta,pch=3)
		abline(a=0,b=1)
		legend(max.x*.1,max.x*.8,c("Mean","Median","Linex"),pch=c(1,2,3))
	}
}

#' calculate parameter estimates and confidence intervals for various loss functions
#' @param Out  Output list from "mcmc_ds.R" 	
#' @return summary.N  list vector, with the first list index indicating species
#' @export
#' @keywords summary
#' @author Paul B. Conn
summary_N<-function(Out){
	n.species=dim(Out$Pred.N)[1]
	summary.N=vector('list',n.species)
	for(isp in 1:n.species){
	  a.linex=calc_linex_a(Out$Pred.N[isp,,],Out$Obs.N[isp,,])$minimum
	  Theta=exp(-a.linex*Out$Post$N[isp,,])
	  Theta=-1/a.linex*log(apply(Theta,2,'mean'))
	  summary.N[[isp]]=list(mean=sum(apply(Out$Post$N[isp,,],2,'mean')),median=sum(apply(Out$Post$N[isp,,],2,'median')),linex=sum(Theta))
  	}  
	summary.N
}

#' Mrds probit detection and related functions
#'
#' For independent observers, probit.fct computes observer-specific detection functions,
#' conditional detection functions, delta dependence function, duplicate detection function (seen by both),
#' and pooled detection function (seen by at least one).
#'
#' The vectors of covariate values can be of different lengths because expand.grid is used to create a
#' dataframe of all unique combinations of the distances and covariate values and the detection and related
#' values are computed for each combination.  The covariate vector observer=1:2 is automatically included.
#'
#' @param x vector of perpendicular distances
#' @param formula linear probit formula for detection using distance and other covariates
#' @param beta parameter values
#' @param rho maximum correlation at largest distance
#' @param ... any number of named vectors of covariates used in the formula
#' @return dat dataframe with distance, observer, any covariates specified in ... and detection probability p,
#' conditional detection probability pc, dupiicate detection dup, pooled detection pool and
#' dependence pc/p=delta.
#' @export
#' @author Jeff Laake
#' @examples
#' test=probit.fct(0:10,~distance,c(1,-.15),.8,size=1:3)
#' par(mfrow=c(1,2))
#' with(test[test$observer==1,],
#' {plot(distance,p,ylim=c(0,1),xlab="Distance",ylab="Detection probability")
#' points(distance,pc,pch=2)
#' points(distance,dup,pch=3)
#' points(distance,pool,pch=4)
#' legend(1,.2,legend=c("Detection","Conditional detection","Duplicate detection","Pooled detection"),pch=1:4,bty="n")
#' plot(distance,delta,xlab="Distance",ylab="Dependence")
#' })

probit.fct=function(x,formula,beta,rho,...)
{
	require(mvtnorm)
#  Create dataframe and apply formula to get design matrix
	dat=expand.grid(distance=x,observer=1:2,...)
	xmat=model.matrix(formula,dat)
#  Make sure length of beta matches number of columns of design matrix
	if(ncol(xmat)!=length(beta))stop("Mismatch between beta and formula")
#  Compute XB and partition for 2 observers
	xbeta=xmat%*%beta
	xbeta1=xbeta[dat$observer==1]
	xbeta2=xbeta[dat$observer==2]
#  Compute rho values
	distance=dat$distance[dat$observer==1]
	rhox=rho*distance/max(distance)
#  Compute detection observer-specific p1,p2 and duplicate p3
	p1=pnorm(xbeta1,0,1)
	p2=pnorm(xbeta2,0,1)
	p3=apply(cbind(xbeta1,xbeta2,rhox),1,function(x)
				pmvnorm(lower=c(-x[1],-x[2]),corr=matrix(c(1,x[3],x[3],1),ncol=2,nrow=2)))
#  Compute conditional detection prob
	p1c2=p3/p2
	p2c1=p3/p1
#  Store values in dataframe
	dat$p[dat$observer==1]=p1
	dat$p[dat$observer==2]=p2
	dat$pc[dat$observer==1]=p1c2
	dat$pc[dat$observer==2]=p2c1
	dat$dup[dat$observer==1]=p3
	dat$dup[dat$observer==2]=p3
	dat$pool[dat$observer==1]=p1+p2-p3
	dat$pool[dat$observer==2]=p1+p2-p3
	dat$delta=dat$pc/dat$p
	return(dat)
}

#' function to convert HierarchicalDS MCMC list vector (used in estimation) into an mcmc object (cf. coda package) 
#' @param MCMC list vector providing MCMC samples for each parameter type 
#' @param N.hab.par see help for mcmc_ds.R
#' @param Cov.par.n see help for mcmc_ds.R
#' @param Hab.names see help for mcmc_ds.R
#' @param Cov.names see help for mcmc_ds.R
#' @param Det.names see help for mcmc_ds.R
#' @param MisID.names see help for mcmc_ds.R
#' @param N.par.misID see help for mcmc_ds.R
#' @param misID.mat see help for mcmc_ds.R
#' @param misID see help for mcmc_ds.R
#' @param fix.tau.nu see help for mcmc_ds.R
#' @param spat.ind see help for mcmc_ds.R
#' @param point.ind see help for mcmc_ds.R
#' @export
#' @keywords MCMC, coda
#' @author Paul B. Conn
convert.HDS.to.mcmc<-function(MCMC,N.hab.par,Cov.par.n,Hab.names,Det.names,Cov.names,MisID.names,N.par.misID=NULL,misID.mat=NULL,fix.tau.nu=FALSE,misID=TRUE,spat.ind=TRUE,point.ind=TRUE){
	require(coda)
	if(misID==TRUE & (is.null(N.par.misID)|is.null(misID.mat)))cat("\n Error: must provide N.par.misID and misID.mat whenever misID=TRUE \n")
	n.species=nrow(MCMC$Hab)
	n.iter=length(MCMC$Hab[1,,1])
	n.col=n.species*2+sum(N.hab.par)+ncol(MCMC$Det)+point.ind+(1-spat.ind)*n.species+(1-fix.tau.nu)*n.species+sum(Cov.par.n)*n.species+misID*sum(N.par.misID)
	Mat=matrix(0,n.iter,n.col)
	Mat[,1:n.species]=t(MCMC$N.tot)
	counter=n.species
	col.names=paste("Abund.sp",c(1:n.species),sep='')
	for(isp in 1:n.species){
	  Mat[,counter+isp]=rowSums(MCMC$G[isp,,]) #total abundance of groups
	  col.names=c(col.names,paste("Groups.sp",isp,sep=''))
    }
	counter=counter+n.species
	for(isp in 1:n.species){  #habitat parameters
	  Mat[,(counter+1):(counter+N.hab.par[isp])]=MCMC$Hab[isp,,1:N.hab.par[isp]]
	  col.names=c(col.names,paste("Hab.sp",isp,Hab.names[[isp]],sep=''))
	  counter=counter+sum(N.hab.par[isp])
  	}
	Mat[,(counter+1):(counter+ncol(MCMC$Det))]=as.matrix(MCMC$Det)
	col.names=c(col.names,paste("Det.",Det.names,sep=''))
	counter=counter+ncol(MCMC$Det)
	if(point.ind==TRUE){
		Mat[,counter+1]=MCMC$cor
		col.names=c(col.names,"rho")
		counter=counter+1
	}
	if(spat.ind==FALSE){
		Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta)
		col.names=c(col.names,paste("tau.eta.sp",c(1:n.species),sep=''))
		counter=counter+n.species
	}
	if(fix.tau.nu==FALSE){
		Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.nu)
		col.names=c(col.names,paste("tau.nu.sp",c(1:n.species),sep=''))
		counter=counter+n.species
	}
	if(is.null(Cov.par.n)==FALSE){
		max.par=max(Cov.par.n)
		for(isp in 1:n.species){
			for(ipar in 1:length(Cov.par.n)){
				Mat[,(counter+1):(counter+Cov.par.n[ipar])]=MCMC$Cov.par[isp,,((ipar-1)*max.par+1):((ipar-1)*max.par+Cov.par.n[ipar])]
				counter=counter+Cov.par.n[ipar]
				col.names=c(col.names,paste("Cov.sp",isp,".",Cov.names[[ipar]],sep=''))
			}
		}
	}
	if(misID==TRUE){
		for(imod in 1:max(misID.mat)){
			Mat[,(counter+1):(counter+N.par.misID[imod])]=MCMC$MisID[,imod,1:N.par.misID[imod]]
			counter=counter+N.par.misID[imod]
			col.names=c(col.names,paste("misID.mod",imod,".",MisID.names[[imod]],sep=''))
		}
	}	
	colnames(Mat)=col.names
	Mat=mcmc(Mat)
  	Mat
}

#' function to plot kernel density estimates for mcmc objects
#' @param MCMC an mcmc object with columns referencing different parameter types (column names are used for plotting labels)
#' @param xaxt if TRUE (default), prints x-axis values
#' @export
#' @keywords MCMC, kernel density, ggplot2
#' @author Paul B. Conn
dens.plot<-function(MCMC,xaxt=TRUE){
	require(ggplot2)
	Plot.df=data.frame(matrix(0,nrow(MCMC)*ncol(MCMC),2))
	Plot.df[,1]=rep(colnames(MCMC),each=nrow(MCMC))
	Plot.df[,2]=as.vector(MCMC)
	colnames(Plot.df)=c("type","Value")
	if(xaxt==FALSE)qplot(Value,data=Plot.df,geom="density")+scale_y_continuous("Density")+opts(axis.text.x=theme_blank(),axis.text.y=theme_blank())+facet_wrap(~type,scale="free")
	else qplot(Value,data=Plot.df,geom="density")+scale_y_continuous("Density")+opts(axis.text.y=theme_blank())+facet_wrap(~type,scale="free")
}


#' function to export posterior summaries from an mcmc object to a table
#' @aliases table.mcmc
#' @S3method table mcmc
#' @method table mcmc
#' @param MCMC An mcmc object with columns referencing different parameter types (column names are used for plotting labels)
#' @param file A file name to ouput to (including path); if null (default), outputs to screen
#' @param type What type of table to produce (either "csv" or "tex")
#' @param a Value to use for credible intervals.  For example, alpha=0.05 results in 95\% credible intervals
#' @export
#' @keywords MCMC, table
#' @author Paul B. Conn
table.mcmc<-function(MCMC,file=NULL,type="csv",a=0.05){
	require(xtable)
	Out.tab=data.frame(matrix(0,ncol(MCMC),5))
	colnames(Out.tab)=c("Parameter","Mean","Median","Lower","Upper")
	MCMC=as.matrix(MCMC)
	Out.tab[,1]=colnames(MCMC)
	Out.tab[,2]=colMeans(MCMC)
	Out.tab[,3]=apply(MCMC,2,'median')
	Out.tab[,4]=apply(MCMC,2,'quantile',a/2)
	Out.tab[,5]=apply(MCMC,2,'quantile',1-a/2)
	if(is.null(file))print(Out.tab)
	else{
		if(type=="csv")write.csv(Out.tab,file=file)
		if(type=="tex"){
			Out.tab=xtable(Out.tab)
			print(Out.tab,file=file)
		}
		if(type!="csv" & type!="tex")cat("\n Error: unknown table type.  No table was printed to file.")
	}
}


#' MCMC output from running example in Hierarchical DS 
#' 
#' @name simdata 
#' @docType data 
#' @author Paul Conn \email{paul.conn@@noaa.gov} 
#' @keywords data 
NULL 