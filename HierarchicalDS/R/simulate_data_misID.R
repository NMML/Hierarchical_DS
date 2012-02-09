#' function to simulate distance sampling data from simple model with increasing abundance
#' intensity, no assumed spatial structure, and point independence
#' @param S number of spatial strata (a single transect is placed in each strata and assumed to cover the whole strata)
#' @param Observers A (2 x S) matrix giving the observer identity for each transect
#' @return a distance sampling dataset
#' @export
#' @keywords distance sampling, simulation
#' @author Paul B. Conn
simulate_data_misID<-function(S,Observers){
	require(mvtnorm)
	set.seed(122412)
	
	#process parameters
	lambda.grp1=3
	lambda.grp2=1
	X.site1=cbind(rep(1,S),log(c(1:S)/S),(log(c(1:S)/S))^2) #covariate on abundance intensity, sp 1
	X.site2=X.site1 #covariate on abundance intensity, sp 1
	Beta.site1=c(log(40),1,0) 
	Beta.site2=c(log(10),-2,-1)
	plot(exp(X.site1%*%Beta.site1))
	points(exp(X.site2%*%Beta.site2))
	
	#detection parameters
	n.bins=5 #n.bins=5 hardwired elsewhere
	#Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.2,0,0,0,0)
	Beta.det=c(1.2,1.0,0.8,-.6,-.3,-.2,-.2,.1,0,.3)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
												#in this version, distance pars are additive (i.e., bin 3 gets bin 2 and bin 3 effect).
	cor.par=0.5 #correlation in max age bin (linear from zero)
		
	#N=rpois(S,0.5*exp(X.site%*%Beta.site))
	N1=round(exp(X.site1%*%Beta.site1))
	N2=round(exp(X.site2%*%Beta.site2))
	cat(paste("\n True G, sp 1 = ",N1,'\n\n'))
	cat(paste("\n True G.tot, sp 1= ",sum(N1),'\n'))
	cat(paste("\n True G, sp 2 = ",N2,'\n\n'))
	cat(paste("\n True G.tot, sp 2= ",sum(N2),'\n'))
	
	Dat=matrix(0,sum(N1)+sum(N2),8)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size
	Confusion=matrix(0,2,3)
	Confusion[1,]=c(1,-1,2)
	Confusion[2,]=c(-1,1,2)
	Conf.model=c(~Observer+Group+Distance+Species,~Observer)
	Conf.param1=c(10,.2,.4,.3,-.1,-.2,-.4,-.8,.5) #parameters for getting it right
	#basically turning misID off with -10 intercept below
	Conf.param2=c(-10,-.2,.2) #parameters for getting an 'unknown'
	
	X=rep(0,length(Beta.det))
	pl=1
	for(i in 1:S){
		cur.Observers=Observers[,i]
		if(N1[i]>0){
			for(j in 1:N1[i]){
				X1=X
				X2=X
				X1[cur.Observers[1]]=1
				X2[cur.Observers[2]]=1
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=sample(c(1:n.bins),1)
				Dat[pl,7]=rpois(1,lambda.grp1)+1
				cur.sp=1
				Dat[pl,8]=cur.sp
				if(Dat[pl,6]>1){
					X1[4:(2+Dat[pl,6])]=1
					X2[4:(2+Dat[pl,6])]=1
				}
				X1[8]=Dat[pl,7]
				X2[8]=Dat[pl,7]
				temp=c(0,0)
				temp[cur.sp]=1
				X1[9:10]=temp
				X2[9:10]=temp
				mu1=X1%*%Beta.det
				mu2=X2%*%Beta.det
				cur.cor=(Dat[pl,6]-1)/(n.bins-1)*cor.par
				Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		if(N2[i]>0){
			for(j in 1:N2[i]){
				X1=X
				X2=X
				X1[cur.Observers[1]]=1
				X2[cur.Observers[2]]=1
				Dat[pl,1]=i
				Dat[pl,2]=cur.Observers[1]
				Dat[pl,3]=cur.Observers[2]
				Dat[pl,6]=sample(c(1:n.bins),1)
				Dat[pl,7]=rpois(1,lambda.grp2)+1
				cur.sp=2
				Dat[pl,8]=cur.sp
				if(Dat[pl,6]>1){
					X1[4:(2+Dat[pl,6])]=1
					X2[4:(2+Dat[pl,6])]=1
				}
				X1[8]=Dat[pl,7]
				X2[8]=Dat[pl,7]
				temp=c(0,0)
				temp[cur.sp]=1
				X1[9:10]=temp
				X2[9:10]=temp
				mu1=X1%*%Beta.det
				mu2=X2%*%Beta.det
				cur.cor=(Dat[pl,6]-1)/(n.bins-1)*cor.par
				Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
				Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
				pl=pl+1
			}
		}
		
	}
	#Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),]
	
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
		Dat2[ipl,8]=Dat[irecord,8]
		Dat2[ipl+1,8]=Dat[irecord,8]
		Dat2[ipl,2]=irecord  #match number
		Dat2[ipl+1,2]=irecord
		ipl=ipl+2
	}
	Dat2=as.data.frame(Dat2)
	colnames(Dat2)=c("Transect","Match","Observer","Obs","Seat","Distance","Group","Species")
	Dat2[,"Observer"]=as.factor(Dat2[,"Observer"])
	Dat2[,"Distance"]=as.factor(Dat2[,"Distance"])
	Dat2[,"Seat"]=as.factor(Dat2[,"Seat"])
	
	Dat=Dat2
	# Now, put in partial observation process
	Ind.sp1=which(Dat[,"Species"]==1)
	Ind.sp2=which(Dat[,"Species"]==2)
	Dat1=Dat[Ind.sp1,]
	X=model.matrix(Conf.model[[1]],data=Dat1)
	MN.logit.right=exp(X%*%Conf.param1)
	X=model.matrix(Conf.model[[2]],data=Dat1)
	MN.logit.unk=exp(X%*%Conf.param2)
	Psi=matrix(0,nrow(Dat1),3)
	Psi[,1]=MN.logit.right/(1+MN.logit.right+MN.logit.unk)
	Psi[,3]=MN.logit.unk/(1+MN.logit.right+MN.logit.unk)
	Psi[,2]=1-Psi[,1]-Psi[,3]
	get_samp<-function(prob)sample(c(1:length(prob)),1,prob=prob)
	Cur.sp=apply(Psi,1,'get_samp')	
	Dat[Ind.sp1,"Species"]=Cur.sp
	
	Dat1=Dat[Ind.sp2,]
	X=model.matrix(Conf.model[[1]],data=Dat1)
	MN.logit.right=exp(X%*%Conf.param1)
	X=model.matrix(Conf.model[[2]],data=Dat1)
	MN.logit.unk=exp(X%*%Conf.param2)
	Psi=matrix(0,nrow(Dat1),3)
	Psi[,2]=MN.logit.right/(1+MN.logit.right+MN.logit.unk)
	Psi[,3]=MN.logit.unk/(1+MN.logit.right+MN.logit.unk)
	Psi[,1]=1-Psi[,2]-Psi[,3]
	get_samp<-function(prob)sample(c(1:length(prob)),1,prob=prob)
	Cur.sp=apply(Psi,1,'get_samp')	
	Dat[Ind.sp2,"Species"]=Cur.sp
	
	Dat[,"Species"]=as.integer(Dat[,"Species"])
	
	Dat=cbind(Dat[,1:4],Dat[,"Species"],Dat[,5:7])
	colnames(Dat)[5]="Species"
    G.true=cbind(N1,N2)
	
	Out=list(Dat=Dat,G.true=G.true)
}

