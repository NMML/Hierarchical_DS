#' function to simulate distance sampling data from simple model with increasing abundance
#' intensity, no assumed spatial structure, and point independence
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @param Observers A (2 x S) matrix giving the observer identity for each transect
#' @param Adj a matrix describing adjacency of neighboring cells
#' @param tau precision parameter of the ICAR process
#' @return a list object with the following slots: 
#' 		"Dat" A distance sampling dataset, 
#' 		"Mapping" a vector specifying which cells are associated with each transect, 
#' 		"Area.trans" a vector specifying the proportion of area of each cell covered by a given transect
#' 		"True.G"	true number of groups per transect
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_data_spatial<-function(S,Observers,Adj,tau){
	require(mvtnorm)
	require(Matrix)
	set.seed(2074278) #for 15 by 15?
	set.seed(2074332)
	
	if(sqrt(S)%%1 >0)cat("\nError: S should be a square number\n")
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q=Matrix(tau*Q)
	#simulate icar process
	Eta=rrw(Q)
	#Eta=rep(0,S)
	
	SP=matrix(Eta,sqrt(S),sqrt(S))
	#SP=0
	#S=10 #number of sites
	#set.seed(207434)
	
	#process parameters
	lambda.grp=3
	Species.prop=c(.5,.3,.2)  #3 species!
	X.site=cbind(rep(1,S),rep(log(c(1:sqrt(S))/sqrt(S)),sqrt(S))) #covariate on abundance intensity
	Beta.site=c(log(200),1) 
	
	#detection parameters
	n.bins=5 #n.bins=5 hardwired elsewhere
	#Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.2,0,0,0,0)
	#Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.1,0,.2,-.4,-.2)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
	Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.1,0,.2,-.4)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
												#in this version, distance pars are additive (i.e., bin 3 gets bin 2 and bin 3 effect).
	cor.par=0.5 #correlation in max age bin (linear from zero)
		
	#sample transects - assume each covers 1/4 of a cell, length=2 cells
	#set.seed(12234
	n.transects=sqrt(S)
	x.base=c(1:n.transects)
	y.base=round(runif(n.transects,0.5,sqrt(S)-.5))
	#y.base=c(1,4,2,1,4)
	#y.base=c(1,1,2,3,4,5,6,7,8,9)
	ids=factor(c(1:n.transects))
	df<-data.frame(id=ids,x=x.base,y=y.base)
	
	#plot expected abundance
	Abund=exp(X.site%*%Beta.site+as.vector(SP)) 
	Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
	colnames(Abund.df)=c("y","x","Abundance")
	require(ggplot2)
	plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,300))+xlab("")+ylab("")
	plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/8,xmax=x+1/8,ymin=y-.5,ymax=y+1.5),fill="maroon")
	plot1
	
	#Determine mapping, fraction of cell occupied by each transect
	Mapping=rep(0,n.transects*2)
	Area.trans=Mapping
	for(itrans in 1:n.transects){
		Mapping[((itrans-1)*2+1):((itrans-1)*2+2)]=(itrans-1)*sqrt(S)+((n.transects-y.base[itrans]):(n.transects-y.base[itrans]+1))
		Area.trans[((itrans-1)*2+1):((itrans-1)*2+2)]=0.25
	}
	
	Exp.grp.abund=exp(X.site%*%Beta.site+as.vector(SP)) 
	cat(paste("\n Expected abundance = ",sum(Exp.grp.abund)))
	#N=rpois(n.transects*2,Area.trans*Exp.grp.abund[Mapping])
	N=rpois(n.transects*2,lambda=Area.trans*Exp.grp.abund[Mapping])
	
	Out=list(Dat=N,Mapping=Mapping,Area.trans=Area.trans,True.G=N,Tot.abund=sum(Exp.grp.abund)) 
	return(Out)
}

