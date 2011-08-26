# simulate data for hierarchical abundance estimation
simulate_data<-function(S){
	S=10 #number of sites
	
	#process parameters
	lambda.grp=3
	X.site=c(1:S)/S #covariate on abundance intensity
	beta.site=100 
	
	#detection parameters
	n.bins=5 #n.bins=5 hardwired elsewhere
	Beta.det=c(1.2,1.0,0.8,-.8,-.6,-.4,-.2,.2)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size
	cor.par=0.5 #correlation in max age bin (linear from zero)
	
	
	N=round(X.site*beta.site)
	
	Dat=matrix(0,sum(N),7)  #rows are site, observer 1 ID, obs 2 ID,  Y_1, Y_2, Distance, Group size
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
			if(Dat[pl,6]>1){
				X1[4:(2+Dat[pl,6])]=1
				X2[4:(2+Dat[pl,6])]=1
			}
			X1[8]=Dat[pl,7]
			X2[8]=Dat[pl,7]
			mu1=X1%*%Beta.det
			mu2=X2%*%Beta.det
			cur.cor=(Dat[pl,6]-1)/(n.bins-1)*cor.par
			Dat[pl,4:5]=rmvnorm(1,c(mu1,mu2),matrix(c(1,cur.cor,cur.cor,1),2,2))
			Dat[pl,4:5]=(Dat[pl,4:5]>0)*1.0
			pl=pl+1
		}
	}
	
	Dat=Dat[which(Dat[,4]>0 | Dat[,5]>0),]
	Dat 
}