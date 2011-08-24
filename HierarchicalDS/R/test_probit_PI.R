# test probit model for multiple observers (pop size assumed known)
require(mvtnorm)
require(truncnorm)

N=1000
Dat=matrix(0,N,4) #Z, Distance, Y1, Y2
Dat[,2]=runif(N)
Dat[,1]=rep(1,N)

rho.true=0.5
beta.obs1=1.2	#det on trans line = 0.885
beta.obs2=0.8   #det on trans line = 0.788
beta.lin=-.05
beta.quad=-2
#Beta=c(beta.obs1,beta.obs2,beta.lin,beta.quad)
Beta=c(beta.obs1,beta.obs2,beta.quad)
n.beta=length(Beta)

X1=matrix(0,N,n.beta) 
X1[,1]=1
#X1[,3]=Dat[,2]
#X1[,4]=Dat[,2]^2
X1[,3]=Dat[,2]^2
X2=X1
X2[,1]=0
X2[,2]=1

Mu1=X1%*%Beta
Mu2=X2%*%Beta

Y.tilde=matrix(0,N,2)

for(i in 1:N){
	Y.tilde[i,]=rmvnorm(1,mean=c(Mu1[i],Mu2[i]),sigma=matrix(c(1,rho.true,rho.true,1),2,2))
}
Dat[,3]=(Y.tilde[,1]>0)
Dat[,4]=(Y.tilde[,2]>0)


#conduct inference
n.iter=10000
MCMC.beta=matrix(0,n.iter,length(Beta))
MCMC.cor=rep(0,n.iter)
MCMC.beta[1,]=Beta
MCMC.cor[1]=0.5
MH.sig=0.05

full_cond_cor<-function(Y.tilde,Mu,rho){ #on log scale!
	Delta=Y.tilde-Mu
	-.5*(nrow(Y.tilde)*log(1-rho^2)+sum(Delta[,1]^2-2*rho*Delta[,1]*Delta[,2]+Delta[,2]^2)/(1-rho^2))
}

st <- Sys.time()

for(iiter in 2:n.iter){
	print(iiter)	
	MCMC.cor[iiter]=MCMC.cor[iiter-1]
	#update beta  #just using basic formulae from chapter 14 of Gelman et al.
	X=rbind(X1-MCMC.cor[iiter]*X2,X2-MCMC.cor[iiter]*X1)
	V.inv <- crossprod(X,X) 	
	Y=c(Y.tilde[,1]-MCMC.cor[iiter]*Y.tilde[,2],Y.tilde[,2]-MCMC.cor[iiter]*Y.tilde[,1])
	M.z <- solve(V.inv, crossprod(X,Y))
	MCMC.beta[iiter,] <- M.z + solve(chol(V.inv), rnorm(n.beta,0,sqrt(1-MCMC.cor[iiter]^2)))
	
	#update cor  (metropolis-hastings)
	Mu1=X1%*%MCMC.beta[iiter,]
	Mu2=X2%*%MCMC.beta[iiter,]
	sig.star=MCMC.cor[iiter]+runif(1,-MH.sig,MH.sig)
	if(sig.star>0 & sig.star<1){		
		logP.old=full_cond_cor(Y.tilde,cbind(Mu1,Mu2),MCMC.cor[iiter])
		Sig.new=matrix(c(1,sig.star,sig.star,1),2,2)
		logP.new=full_cond_cor(Y.tilde,cbind(Mu1,Mu2),sig.star)
		if(runif(1)<exp(logP.new-logP.old)){
			MCMC.cor[iiter]=sig.star
		}
	}
	
	#update Y.tilde parms
	EY1=Mu1+MCMC.cor[iiter]*(Y.tilde[,2]-Mu2)
	Y.tilde[,1] <- rtruncnorm(N, a=ifelse(Dat[,3]==0,-Inf,0), b=ifelse(Dat[,3]==0,0,Inf), EY1, sqrt(1-MCMC.cor[iiter]^2))
	EY2=Mu2+MCMC.cor[iiter]*(Y.tilde[,1]-Mu1)
	Y.tilde[,2] <- rtruncnorm(N, a=ifelse(Dat[,4]==0,-Inf,0), b=ifelse(Dat[,4]==0,0,Inf), EY2, sqrt(1-MCMC.cor[iiter]^2))

	
}
tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))
cat("\n time elapsed (seconds): ", tpi, "\n")

