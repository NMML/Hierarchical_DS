#' function to sample from a specified probability density function
#' @param n number of samples desired
#' @param pdf probability density function (pois1, poisson, normal, unif.disc, unif.cont)
#' @param cur.par a length two vector giving parameters for the specified distribution; only the first is used for single parameter distributions
#' @return a vector of length n samples from the desired distribution 
#' @export
#' @keywords probability density
#' @author Paul B. Conn
switch_sample<-function(n,pdf,cur.par){
	switch(pdf,
			pois1=rpois(n,cur.par[1])+1,
			poisson=rpois(n,cur.par[1]),
			normal=rnorm(n,cur.par[1],cur.par[2]),
			unif.disc=sample(cur.par[1]:cur.par[2],n,replace=TRUE),
			unif.cont=runif(n,cur.par[1],cur.par[2])
	)
}

#' function to stack data (going from three dimensional array to a two dimensional array including only "existing" animals
#' @param Data three-d dataset
#' @param G.transect	current number of groups of animals in each transect (vector)
#' @param n.transects	number of transects
#' @param stacked.names	column names for new stacked dataset
#' @param factor.ind	a vector of indicator variables (1 = factor/categorical variable, 0 = continuous variable)
#' @return a stacked dataset
#' @export
#' @keywords stack data
#' @author Paul B. Conn
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

#' function to produce a design matrix given
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
generate_inits<-function(DM.hab,DM.det,G.transect,Area.trans,Area.hab,Mapping,point.ind,spat.ind,grp.mean){		
	Par=list(det=rnorm(ncol(DM.det),0,1),hab=rep(0,ncol(DM.hab)),cor=ifelse(point.ind,runif(1,0,1),0),
			Nu=log(max(G.transect)/mean(Area.trans)*exp(rnorm(length(Area.hab)))),Eta=rnorm(length(Area.hab)),
			tau.eta=runif(1,0.5,2),tau.nu=runif(1,0.5,2))
	Par$hab[1]=mean(G.transect)/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(1,0,1))
	Par$G=exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu)))
	Par$N=Par$G+rpois(length(Par$G),grp.mean*Par$G)
	if(spat.ind==1)Par$Eta=0*Par$Eta
	Par
}

#' compute the first derivative of log_lambda likelihood component for Langevin-Hastings
#' @param index 	which element to compute the derivative for (for use with apply functions)
#' @param Mu	expected value 
#' @param Nu	current observed value	
#' @param N	    number of groups in each transect
#' @param var.nu	variance of the overdispersion process
#' @return a gradient value
#' @export
#' @keywords gradient, Langevin-Hastings
#' @author Paul B. Conn
log_lambda_gradient<-function(index,Mu,Nu,N,var.nu){
	return(-(Nu[index]-Mu[index])/var.nu+N[index]-exp(Nu[index]))
}

#' compute the likelihood for nu parameters
#' @param Log.lambda 	which element to compute the derivative for (for use with apply functions)
#' @param DM	the design matrix
#' @param Beta	linear predictor parameters for the log of abundance intensity	
#' @param Eta	a vector of spatial random effects
#' @param SD	standard deviation of the overdispersion process
#' @param N		a vector giving the current iteration's number of groups in the area covered by each transect
#' @return the log likelihood associated with the data and the current set of parameters 
#' @export
#' @keywords log likelihood
#' @author Paul B. Conn
log_lambda_log_likelihood<-function(Log.lambda,DM,Beta,Eta=0,SD,N){
	Pred.log.lam=DM%*%Beta+Eta	
	logL=sum(dnorm(Log.lambda,Pred.log.lam,SD,log=1)) #year 1
	logL=logL+sum(N*Log.lambda-exp(Log.lambda))
	return(logL)
} 	

