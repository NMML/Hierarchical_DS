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

#' function to sample from hyperpriors of a specified probability density function
#' @param pdf probability density function (pois1, poisson, normal, unif.disc, unif.cont)
#' @param cur.par a vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @return a vector of length n samples from the desired distribution 
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_sample_prior<-function(pdf,cur.par){
	switch(pdf,
			pois1=rgamma(1,cur.par[1],cur.par[2]),
			poisson=rgamma(1,cur.par[1],cur.par[2]),
			pois1_ln=c(rnorm(1,cur.par[1],cur.par[2]),runif(1,0,cur.par[3])),
			poisson_ln=c(rnorm(1,cur.par[1],cur.par[2]),runif(1,0,cur.par[3])),
			multinom=rdirichlet(1,cur.par)
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
	Stacked=as.data.frame(Data[1,1:2,])
	for(itrans in 1:n.transects){
		if(Obs.transect[itrans]>0)Stacked=rbind(Stacked,Data[itrans,1:Obs.transect[itrans],])
	}
	Stacked=Stacked[-c(1,2),]
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
#' @param factor.ind	a vector of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
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
generate_inits<-function(DM.hab,G.transect,Area.trans,Area.hab,Mapping,spat.ind){		
	Par=list(hab=rep(0,ncol(DM.hab)),
			Nu=log(max(G.transect)/mean(Area.trans)*exp(rnorm(length(Area.hab)))),Eta=rnorm(length(Area.hab)),
			tau.eta=runif(1,0.5,2),tau.nu=runif(1,0.5,2))
	Par$hab[1]=mean(G.transect)/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(1,0,1))
	Par$G=exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu)))
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
#' @param Par  Posterior samples from MCMC 	
#' @param Mapping	Index of which areas belong to which transects
#' @param Area	Propotional area of a transect wrt addressed grid cell
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

#' plot 'observed' versus predicted values for 
#' @param Out  Output list from "mcmc_ds.R" 	
#' @return NULL 
#' @export
#' @keywords diagnostics, plot
#' @author Paul B. Conn
plot_obs_pred<-function(Out){
	a.linex=calc_linex_a(Out$Pred.G,Out$Obs.G)$minimum
	max.x=max(c(apply(Out$Obs.G,2,'mean'),apply(Out$Pred.G,2,'mean')))
	plot(apply(Out$Obs.G,2,'mean'),apply(Out$Pred.G,2,'mean'),pch=1,xlim=c(0,max.x),ylim=c(0,max.x),xlab="Observed",ylab="Predicted")
	points(apply(Out$Obs.G,2,'mean'),apply(Out$Pred.G,2,'median'),pch=2)
	Theta=exp(-a.linex*Out$Pred.G)
	Theta=-1/a.linex*log(apply(Theta,2,'mean'))
	points(apply(Out$Obs.G,2,'mean'),Theta,pch=3)
	abline(a=0,b=1)
	legend(max.x*.1,max.x*.8,c("Mean","Median","Linex"),pch=c(1,2,3))
}

#' calculate parameter estimates and confidence intervals for various loss functions
#' @param Out  Output list from "mcmc_ds.R" 	
#' @return summary.N 
#' @export
#' @keywords summary
#' @author Paul B. Conn
summary_N<-function(Out){
	a.linex=calc_linex_a(Out$Pred.G,Out$Obs.G)$minimum
	Theta=exp(-a.linex*Out$MCMC$G)
	Theta=-1/a.linex*log(apply(Theta,2,'mean'))
	summary.N=list(mean=sum(apply(Out$MCMC$G,2,'mean')),median=sum(apply(Out$MCMC$G,2,'median')),linex=sum(Theta))
	summary.N
}

