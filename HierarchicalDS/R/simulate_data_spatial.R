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
	set.seed(2074278)
	
	if(sqrt(S)%%1 >0)cat("\nError: S should be a square number\n")
	Q=-Adj
	diag(Q)=apply(Adj,2,'sum')
	Q=Matrix(tau*Q)
	#simulate icar process
	Eta=rrw(Q)
	#Eta=rep(0,S)
	
	SP=matrix(Eta,sqrt(S),sqrt(S))
	#S=10 #number of sites
	
	#process parameters
	lambda.grp=3
	Species.prop=c(.1,.2,.3,.4)  #4 species!
	X.site=cbind(rep(1,S),rep(log(c(1:sqrt(S))/sqrt(S)),sqrt(S))) #covariate on abundance intensity
	Beta.site=c(log(100),1) 
	
	#detection parameters
	n.bins=5 #n.bins=5 hardwired elsewhere
	#Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.2,0,0,0,0)
	Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.1,0,.2,-.4,-.2)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
												#in this version, distance pars are additive (i.e., bin 3 gets bin 2 and bin 3 effect).
	cor.par=0.5 #correlation in max age bin (linear from zero)
		
	#sample transects - assume each covers 1/4 of a cell, length=2 cells
	set.seed(111124)
	n.transects=sqrt(S)
	x.base=c(1:n.transects)
	y.base=runif(n.transects,0.5,sqrt(S)-2.5)
	ids=factor(c(1:n.transects))
	df<-data.frame(id=ids,x=x.base,y=y.base)
	
	#plot expected abundance
	Abund=4*exp(X.site%*%Beta.site+as.vector(SP))  #sum=11971; expected total abundance =11971*4=47885
	Abund.df=data.frame(cbind(rep(c(1:sqrt(S)),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
	colnames(Abund.df)=c("y","x","Abundance")
	require(ggplot2)
	plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,750))+xlab("")+ylab("")
	plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/8,xmax=x+1/8,ymin=y,ymax=y+3),fill="maroon")
	plot1
	
	#Determine mapping, fraction of cell occupied by each transect
	Mapping=rep(0,n.transects*4)
	Area.trans=Mapping
	for(itrans in 1:n.transects){
		Mapping[((itrans-1)*4+1):((itrans-1)*4+4)]=(itrans-1)*15+c((n.transects-floor(y.base[itrans]+.5)-2):(n.transects-floor(y.base[itrans]+.5)+1))
		Area.trans[((itrans-1)*4+1):((itrans-1)*4+4)]=.25*c(ceiling(y.base[itrans]+0.5)-y.base[itrans]-0.5,1,1,y.base[itrans]+0.5-floor(y.base[itrans]+0.5))
	}
	
	Exp.grp.abund=exp(X.site%*%Beta.site+as.vector(SP))
	N=rpois(n.transects*4,Area.trans*Exp.grp.abund[Mapping])
	
	Dat=matrix(0,sum(N),8)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size,species
	X=rep(0,length(Beta.det))
	pl=1
	for(i in 1:(n.transects*4)){
		cur.Observers=Observers[,i]
		if(N[i]>0){
			for(j in 1:N[i]){
				if(is.na(cur.Observers[2])){
					X1=X
					X1[cur.Observers[1]]=1
					Dat[pl,1]=i
					Dat[pl,2]=cur.Observers[1]
					Dat[pl,3]=cur.Observers[2]
					Dat[pl,6]=sample(c(1:n.bins),1)
					Dat[pl,7]=rpois(1,lambda.grp)+1
					cur.sp=sample(c(1,2,3,4),1,prob=Species.prop)
					Dat[pl,8]=cur.sp
					if(Dat[pl,6]>1)X1[4:(2+Dat[pl,6])]=1
					X1[8]=Dat[pl,7]
					temp=c(0,0,0,0)
					temp[cur.sp]=1
					X1[9:12]=temp
					mu1=X1%*%Beta.det
					Dat[pl,4]=rnorm(1,mu1,1)
					Dat[pl,4]=(Dat[pl,4]>0)*1.0
					Dat[pl,5]=NA
					pl=pl+1					
				}
				else{
					X1=X
					X2=X
					X1[cur.Observers[1]]=1
					X2[cur.Observers[2]]=1
					Dat[pl,1]=i
					Dat[pl,2]=cur.Observers[1]
					Dat[pl,3]=cur.Observers[2]
					Dat[pl,6]=sample(c(1:n.bins),1)
					Dat[pl,7]=rpois(1,lambda.grp)+1
					cur.sp=sample(c(1,2,3,4),1,prob=Species.prop)
					Dat[pl,8]=cur.sp
					if(Dat[pl,6]>1){
						X1[4:(2+Dat[pl,6])]=1
						X2[4:(2+Dat[pl,6])]=1
					}
					X1[8]=Dat[pl,7]
					X2[8]=Dat[pl,7]
					temp=c(0,0,0,0)
					temp[cur.sp]=1
					X1[9:12]=temp
					X2[9:12]=temp
					mu1=X1%*%Beta.det
					mu2=X2%*%Beta.det
					cur.cor=(Dat[pl,6]-1)/(n.bins-1)*cor.par
					Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
					Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
					pl=pl+1
				}
			}
		}
	}
	
	Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),]
	
	#put things in "Jay's" format
	Dat2=rbind(Dat,Dat)
	ipl=1
	for(irecord in 1:nrow(Dat)){
		if(is.na(Dat[irecord,3])){ #one observer
			Dat2[ipl,1]=Dat[irecord,1]
			Dat2[ipl,3]=Dat[irecord,2]
			Dat2[ipl,4]=Dat[irecord,4]
			Dat2[ipl,5]=1  #observer covariate that has no effect
			Dat2[ipl,6]=Dat[irecord,6]
			Dat2[ipl,7]=Dat[irecord,7]
			Dat2[ipl,8]=Dat[irecord,8]
			Dat2[ipl,2]=irecord  #match number
			ipl=ipl+1
		}
		else{ #two observers
			Dat2[ipl,1]=Dat[irecord,1]
			Dat2[ipl+1,1]=Dat[irecord,1]
			Dat2[ipl,3]=Dat[irecord,2]
			Dat2[ipl+1,3]=Dat[irecord,3]
			Dat2[ipl,4]=Dat[irecord,4]
			Dat2[ipl+1,4]=Dat[irecord,5]
			Dat2[ipl,5]=1  #observer covariate that has no effect
			Dat2[ipl+1,5]=0
			Dat2[ipl,6]=Dat[irecord,6]
			Dat2[ipl+1,6]=Dat[irecord,6]
			Dat2[ipl,7]=Dat[irecord,7]
			Dat2[ipl+1,7]=Dat[irecord,7]
			Dat2[ipl,8]=Dat[irecord,8]
			Dat2[ipl+1,8]=Dat[irecord,8]
			Dat2[ipl,2]=irecord  #match number
			Dat2[ipl+1,2]=irecord
			ipl=ipl+2
		}
	}
	Dat2=as.data.frame(Dat2[1:(ipl-1),])
	colnames(Dat2)=c("Transect","Match","Observer","Obs","Seat","Distance","Group","Species")
	Dat2[,"Observer"]=as.factor(Dat2[,"Observer"])
	Dat2[,"Distance"]=as.factor(Dat2[,"Distance"])
	Dat2[,"Seat"]=as.factor(Dat2[,"Seat"])
	Dat2[,"Species"]=as.factor(Dat2[,"Species"])
	Out=list(Dat=Dat2,Mapping=Mapping,Area.trans=Area.trans,True.G=N,Tot.abund=Abund) 
	Out
}