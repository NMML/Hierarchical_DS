

# simulate & analyze data for one site (two observers); 

n.sims=50
Est=rep(0,n.sims)
lambda.true=100
log.lambda.true=log(lambda.true)
#N.true=rpois(n.sims,lambda.true)
N.true=rep(lambda.true,n.sims)
M.over.N=2 #use different values of M in data augmentation
M.sim=M.over.N*N.true
N.est1=matrix(0,n.sims,2)  #keep track of posterior mean and mode of Binomial data aug
N.est2=matrix(0,n.sims,2)  #keep track of " " for Poisson rev jump algorithm
p1=0.9  #probably of detection on trackline for obs 1
P1.est=N.est1
P2.est=N.est1
p2=0.8  # "" obs 2
sigma.det1=.4 	#assume exp detection function in form of exp(-distance/alpha)^1/beta;  
sigma.det2=.3 	#assume exp detection function in form of exp(-distance/alpha)^1/beta;  


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
	MCMC=list(N=rep(0,Control$iter-Control$burnin),sigma1=rep(0,Control$iter-Control$burnin),sigma2=rep(0,Control$iter-Control$burnin),p1=rep(0,Control$iter-Control$burnin),p2=rep(0,Control$iter-Control$burnin))
	N=sum(Dat[,1])
	n.obs=sum((Dat[,3]+Dat[,4])>0)
	sigma1=Inits$sigma1
	sigma2=Inits$sigma2
	p1=Inits$p1
	p2=Inits$p2
	norm.const1=dnorm(0,0,sigma1)
	norm.const2=dnorm(0,0,sigma2)
	Sample=c(-20:20)
	Sample=Sample[-21]
	max.a=max(Sample)
	Accept=list(N=0,sigma1=0,sigma2=0,dist=0)
	
	#begin MCMC
	for(iiter in 1:Control$iter){
		if((iiter%%1000)==1)cat(paste('\n ', iiter))
		
		#update distances
		if(N>n.obs){
		  for(iind in (n.obs+1):N){
		    dist.star=Dat[iind,2]+runif(1,-Control$MH.dist.sd,Control$MH.dist.sd)
		    if(dist.star>0 & dist.star<1){
		      LogL.star=log(1-p1*dnorm(dist.star,0,sigma1)/norm.const1)+log(1-p2*dnorm(dist.star,0,sigma2)/norm.const2)
		      LogL.old=log(1-p1*dnorm(Dat[iind,2],0,sigma1)/norm.const1)+log(1-p2*dnorm(Dat[iind,2],0,sigma2)/norm.const2)
              if(runif(1)<exp(LogL.star-LogL.old)){
                Dat[iind,2]=dist.star
                Accept$dist=Accept$dist+1
			  }
		    }
          }
        }
      
    	Dat[(N+1):M,2]=runif(M-N)
    		
		#additions/deletions
		a=sample(Sample,1)
		if(a>0){  #current a results in an addition
			if((N+a)>M)cat('\n Proposal for N exceeds M!')
			else{
				tmp.sum=0
				for(i in 1:a){
					tmp.sum=tmp.sum-log(N-n.obs+i)+log(1-p1*dnorm(Dat[N+i,2],0,sigma1)/norm.const1)+log(1-p2*dnorm(Dat[N+i,2],0,sigma2)/norm.const2)
				}
				MH.prob=exp(a*log.lambda.true+tmp.sum)
				if(runif(1)<MH.prob){
				  Dat[(N+1):(N+a),1]=1
				  N=N+a
				  Accept$N=Accept$N+1
				}
			}			
		}
		else{  #current a results in a subtraction (a<0)
			if((N+a)>=n.obs){  #can't delete actual observations
				tmp.sum=0
				for(i in 0:(a+1)){
					tmp.sum=tmp.sum+log(N-n.obs+i)-log(1-p1*dnorm(Dat[N+i,2],0,sigma1)/norm.const1)-log(1-p2*dnorm(Dat[N+i,2],0,sigma2)/norm.const2)
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
	  	sigma1.star=sigma1+runif(1,-Control$MH.sig.sd,Control$MH.sig.sd)
	  	if(sigma1.star>0){
	    	norm.const.star=dnorm(0,0,sigma1.star)
     		p.star=p1*dnorm(Dat[1:N,2],0,sigma1.star)/norm.const.star
      		p.old=p1*dnorm(Dat[1:N,2],0,sigma1)/norm.const1 
      		logL.old=sum(Dat[1:N,3]*log(p.old))+sum((1-Dat[1:N,3])*log(1-p.old))-log(sigma1)
			logL.star=sum(Dat[1:N,3]*log(p.star ))+sum((1-Dat[1:N,3])*log(1-p.star))-log(sigma1.star)
     		if(runif(1)<exp(logL.star-logL.old)){
        		Accept$sigma1=Accept$sigma1+1
        		sigma1=sigma1.star
        		norm.const1=norm.const.star        
      		}
    	}

		sigma2.star=sigma2+runif(1,-Control$MH.sig.sd,Control$MH.sig.sd)
		if(sigma2.star>0){
			norm.const.star=dnorm(0,0,sigma2.star)
			p.star=p2*dnorm(Dat[1:N,2],0,sigma2.star)/norm.const.star
			p.old=p2*dnorm(Dat[1:N,2],0,sigma2)/norm.const2 
			logL.old=sum(Dat[1:N,4]*log(p.old))+sum((1-Dat[1:N,4])*log(1-p.old))-log(sigma2)
			logL.star=sum(Dat[1:N,4]*log(p.star))+sum((1-Dat[1:N,4])*log(1-p.star))-log(sigma2.star)
			if(runif(1)<exp(logL.star-logL.old)){
				Accept$sigma2=Accept$sigma2+1
				sigma2=sigma2.star
				norm.const2=norm.const.star        
			}
		}
		
		#update p(det) on trackline
		p1.star=p1+runif(1,-Control$MH.p.sd,Control$MH.p.sd)
		if(p1.star>0 & p1.star<1){
			p.star=p1.star*dnorm(Dat[1:N,2],0,sigma1)/norm.const1
			p.old=p1*dnorm(Dat[1:N,2],0,sigma1)/norm.const1
			logL.old=sum(Dat[1:N,3]*log(p.old))+sum((1-Dat[1:N,3])*log(1-p.old))
			logL.star=sum(Dat[1:N,3]*log(p.star))+sum((1-Dat[1:N,3])*log(1-p.star))
			if(runif(1)<exp(logL.star-logL.old)){
				Accept$p1=Accept$p1+1
				p1=p1.star
			}
		}	
		
		p2.star=p2+runif(1,-Control$MH.p.sd,Control$MH.p.sd)
		if(p2.star>0 & p2.star<1){
			p.star=p2.star*dnorm(Dat[1:N,2],0,sigma2)/norm.const2
			p.old=p2*dnorm(Dat[1:N,2],0,sigma2)/norm.const2
			logL.old=sum(Dat[1:N,4]*log(p.old))+sum((1-Dat[1:N,4])*log(1-p.old))
			logL.star=sum(Dat[1:N,4]*log(p.star))+sum((1-Dat[1:N,4])*log(1-p.star))
			if(runif(1)<exp(logL.star-logL.old)){
				Accept$p2=Accept$p2+1
				p2=p2.star
			}
		}		
		
		#store results if applicable
   		if(iiter>Control$burnin){
      		MCMC$N[iiter-Control$burnin]=N
     		MCMC$sigma1[iiter-Control$burnin]=sigma1
			MCMC$sigma2[iiter-Control$burnin]=sigma2
			MCMC$p1[iiter-Control$burnin]=p1
			MCMC$p2[iiter-Control$burnin]=p2
    	}
	}
	Out=list(MCMC=MCMC,Accept=Accept)
}

for(isim in 1:n.sims){
	#1) generate data
	Distances=runif(N.true[isim])
	Det.p1=p1*dnorm(Distances,0,sigma.det1)/dnorm(0,0,sigma.det1)
	Det.p2=p2*dnorm(Distances,0,sigma.det2)/dnorm(0,0,sigma.det2)
	Y1=rbinom(N.true[isim],1,Det.p1)
	Y2=rbinom(N.true[isim],1,Det.p2)
	Y=Y1+Y2
	Obs.dist=Distances[which(Y>0)]
	Y1=Y1[which(Y>0)]
	Y2=Y2[which(Y>0)]
	Not.obs.dist=Distances[which(Y==0)]
	
	#2) Assemble data 
	Dat=matrix(0,M.sim[isim],4)  #first column - z, second column - x, third column - Y1, fourth - Y2
	Dat[1:N.true[isim],1]=1  #start at true values
	n.obs=sum(Y>0)
	Dat[1:n.obs,3]=Y1
	Dat[1:n.obs,4]=Y2
	Dat[1:n.obs,2]=Obs.dist
	Dat[(n.obs+1):N.true[isim],2]=Not.obs.dist
	Dat[(N.true[isim]+1):M.sim[isim],2]=runif(M.sim[isim]-N.true[isim])

	
	#3) run pois data 
	Control=list(iter=15000,burnin=0,MH.sig.sd=0.07,MH.dist.sd=0.2,MH.p.sd=0.03)
	Inits=list(sigma1=sigma.det1,sigma2=sigma.det2,p1=1,p2=1)    #p1 and p2 are probability of detection on the line
	Out=pois_data_aug(Dat,Control,Inits)
    plot(Out$MCMC$N,type="l")
	N.est2[isim,1]=mean(Out$MCMC$N)
	N.est2[isim,2]=median(Out$MCMC$N)
	P1.est[isim,1]=mean(Out$MCMC$p1)
	P1.est[isim,2]=median(Out$MCMC$p1)
	P2.est[isim,1]=mean(Out$MCMC$p2)
	P2.est[isim,2]=median(Out$MCMC$p2)
	Out$Accept$dist=Out$Accept$dist/(Control$iter*(mean(Out$MCMC$N)-n.obs)) #approximate acceptance rate for distances
}
