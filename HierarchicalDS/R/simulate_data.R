#' function to simulate distance sampling data from simple model with increasing abundance
#' intensity, no assumed spatial structure, and point independence
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @return a distance sampling dataset
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_data<-function(S){
	require(mvtnorm)
	#S=10 #number of sites
	
	#process parameters
	lambda.grp=3
	Species.prop=c(.1,.2,.3,.4)  #4 species!
	X.site=cbind(rep(1,S),log(c(1:S)/S)) #covariate on abundance intensity
	Beta.site=c(log(100),1) 
	
	#detection parameters
	n.bins=5 #n.bins=5 hardwired elsewhere
	Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.2,0,0,0,0)
	#Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.2,0,.2,-.4,-.2)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
												#in this version, distance pars are additive (i.e., bin 3 gets bin 2 and bin 3 effect).
	cor.par=0.5 #correlation in max age bin (linear from zero)
		
	N=rpois(S,0.5*exp(X.site%*%Beta.site))
	#N=exp(X.site%*%Beta.site)
	
	Dat=matrix(0,sum(N),8)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size,species
	X=rep(0,length(Beta.det))
	pl=1
	for(i in 1:S){
		Observers=sample(c(1:3),2)
		for(j in 1:N[i]){
			X1=X
			X2=X
			X1[Observers[1]]=1
			X2[Observers[2]]=1
			Dat[pl,1]=i
			Dat[pl,2]=Observers[1]
			Dat[pl,3]=Observers[2]
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
	
	Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),]
	
	#put things in "Jay's" format
	Dat2=rbind(Dat,Dat)
	ipl=1
	for(irecord in 1:nrow(Dat)){
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
		Dat2[ipl,2]=irecord  #match number
		Dat2[ipl+1,2]=irecord
		ipl=ipl+2
	}
	Dat2=as.data.frame(Dat2)
	colnames(Dat2)=c("Transect","Match","Observer","Obs","Seat","Distance","Group","Species")
	Dat2[,"Observer"]=as.factor(Dat2[,"Observer"])
	Dat2[,"Seat"]=as.factor(Dat2[,"Seat"])
	Dat2 
}