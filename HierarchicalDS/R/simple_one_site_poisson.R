

# simulate & analyze data for one site (one observer); compare typical dat aug (andy) (Bin approx to Poisson) w/ 
# exact poisson algorithm (addition, deletion, translation)

n.sims=1
Est=rep(0,n.sims)
lambda.true=1000
log.lambda.true=log(lambda.true)
#N.true=rpois(n.sims,lambda.true)
N.true=rep(lambda.true,n.sims)
M.over.N=2 #use different values of M in data augmentation
M.sim=M.over.N*N.true
N.est1=matrix(0,n.sims,2)  #keep track of posterior mean and mode of Binomial data aug
N.est2=matrix(0,n.sims,2)  #keep track of " " for Poisson rev jump algorithm

sigma.det=.4 	#assume exp detection function in form of exp(-distance/alpha)^1/beta;  



bin_data_aug<-function(Dat,Control){ #assumes uniform priors
	M=nrow(Dat)
	N=rep(0,Control$iter-Control$burnin)
	N[1]=sum(Data[,1])
	psi=N[1]/M
	#update z, d

}

sample_distance<-function(a,sigma,norm.const){ 
	#draw distance(s) from [D | y=0,theta]
	D=rep(0,a)
	for(i in 1:a){
		D[i]=-1
		while(D[i]<0){
			tmp=runif(1)
			if(runif(1)<(1-dnorm(tmp,0,sigma)/norm.const))D[i]=tmp
		}
			
	}
	D
}

pois_data_aug<-function(Dat,Control,Inits){ #assumes uniform priors
	M=nrow(Dat)
	MCMC=list(N=rep(0,Control$iter-Control$burnin),sigma=rep(0,Control$iter-Control$burnin))
	N=sum(Dat[,1])
	n.obs=sum(Dat[,2])
	sigma=Inits$sigma
	norm.const=dnorm(0,0,sigma)
	Sample=c(-20:20)
	Sample=Sample[-21]
	max.a=max(Sample)
	Accept=list(N=0,sigma=0)
	
	#begin MCMC
	for(iiter in 1:Control$iter){
		if((iiter%%1000)==1)cat(paste('\n ', iiter))
		
		#update distances
    Dat[(n.obs+1):N,3]=sample_distance(a=(N-n.obs),sigma=sigma,norm.const=norm.const)
    Dat[(N+1):M,3]=runif(M-N)
    		
		#additions/deletions
		a=sample(Sample,1)
		if(a>0){  #current a results in an addition
			if((N+a)>M)cat('\n Proposal for N exceeds M!')
			else{
				#D.prop=sample_distance(a,sigma,norm.const)
				tmp.sum=0
				for(i in 1:a){
					tmp.sum=tmp.sum-log(N-n.obs+i)+log(1-dnorm(Dat[N+i,3],0,sigma)/norm.const)
				}
				MH.prob=exp(a*log.lambda.true+tmp.sum)
				if(runif(1)<MH.prob){
				  Dat[(N+1):(N+a),1]=1
				  #Dat[(N+1):(N+a),3]=D.prop
				  N=N+a
				  Accept$N=Accept$N+1
				}
			}			
		}
		else{  #current a results in a subtraction (a<0)
			if((N+a)>=n.obs){  #can't delete actual observations
				tmp.sum=0
				for(i in 0:(a+1)){
					tmp.sum=tmp.sum+log(N-n.obs+i)-log(1-dnorm(Dat[N+i,3],0,sigma)/norm.const)
				}
				MH.prob=exp(a*log.lambda.true+tmp.sum)
			  if(runif(1)<MH.prob){
					Dat[(N:(N+a+1)),1]=0
					N=N+a
					Accept$N=Accept$N+1
				}
			}
		}
		
		#translations (probably don't need right now)
		
		
		#update sigma (half normal det. param); assume jeffrey's prior of 1/sig
	  sigma.star=sigma+runif(1,-Control$MH.sig.sd,Control$MH.sig.sd)
	  if(sigma.star>0){
	    norm.const.star=dnorm(0,0,sigma.star)
      p.star=dnorm(Dat[1:N,3],0,sigma.star)/norm.const.star
      p.old=dnorm(Dat[1:N,3],0,sigma)/norm.const 
      logL.old=sum(log(p.old[1:n.obs]))-log(sigma)
      if(N>n.obs)logL.old=logL.old+sum(log(1-p.old[(n.obs+1):N]))
      logL.star=sum(log(p.star[1:n.obs]))-log(sigma.star)
      if(N>n.obs)logL.star=logL.star+sum(log(1-p.star[(n.obs+1):N]))
      if(runif(1)<exp(logL.star-logL.old)){
        Accept$sigma=Accept$sigma+1
        sigma=sigma.star
        norm.const=norm.const.star        
      }
    }
	
    if(iiter>Control$burnin){
      MCMC$N[iiter-Control$burnin]=N
      MCMC$sigma[iiter-Control$burnin]=sigma
    }
	}
	Out=list(MCMC=MCMC,Accept=Accept)
}

for(isim in 1:n.sims){
	#1) generate data
	Distances=runif(N.true[isim])
	Det.p=dnorm(Distances,0,sigma.det)/dnorm(0,0,sigma.det)
	Y=rbinom(N.true[isim],1,Det.p)
	Obs.dist=Distances[which(Y==1)]
	Not.obs.dist=Distances[which(Y==0)]
	
	#2) Assemble data 
	Dat=matrix(0,M.sim[isim],3)  #first column - z, second column - y, third column - d
	Dat[1:N.true[isim],1]=1  #start at true values
	Dat[1:sum(Y),2]=1
	Dat[1:sum(Y),3]=Obs.dist
	Dat[(sum(Y)+1):N.true[isim],3]=Not.obs.dist
	Dat[(N.true[isim]+1):M.sim[isim],3]=runif(M.sim[isim]-N.true[isim])

	
	#3) run pois data 
	Control=list(iter=15000,burnin=0,MH.sig.sd=0.04)
	Inits=list(sigma=sigma.det)
	Out=pois_data_aug(Dat,Control,Inits)
  plot(Out$MCMC$N,type="l")
	N.est2[isim,1]=mean(Out$MCMC$N)
	N.est2[isim,2]=median(Out$MCMC$N)
	
}
