
require(Matrix)
require(ggplot2)
S=100
set.seed(1093284)
tau=10
tau.e=10

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

Adj=matrix(0,S,S)
for(i in 2:S){
	Adj[i-1,i]=1
	Adj[i,i-1]=1
}
Q=-Adj
diag(Q)=apply(Adj,2,'sum')
Q=Matrix(tau*Q)

Eta=rrw(Q)
#Eta=rep(0,S)
	
SP=matrix(Eta,sqrt(S),sqrt(S))

	
X.site=cbind(rep(1,S),rep(log(c(1:sqrt(S))/sqrt(S)),sqrt(S))) #covariate on abundance intensity
Beta.site=c(log(100),1) 

Dat=matrix(X.site%*%Beta.site+Eta+rnorm(S,0,sqrt(1/tau.e)),sqrt(S),sqrt(S))
	
Abund.df=data.frame(cbind(rep(c(1:sqrt(S)),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Dat))))
colnames(Abund.df)=c("y","x","Abundance")
require(ggplot2)
plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,6))+xlab("")+ylab("")
plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/8,xmax=x+1/8,ymin=y-.5,ymax=y+1.5),fill="maroon")
plot1

### first, try to estimate beta using traditional ICAR approach
n.iter=1000
MCMC=list(Beta=matrix(0,2,n.iter),tau=rep(0,n.iter),tau.e=rep(0,n.iter))
MCMC$Beta[1,1]=4
MCMC$Beta[2,1]=0.5

MCMC$tau[1]=5
MCMC$tau.e[1]=.2

X=X.site
Dat=as.vector(Dat)

V.inv <- crossprod(X,X) 
Eta=rnorm(S,0,sqrt(1/tau))
n.obs.tau=S-1
a.eta=.01
b.eta=.01

for(iiter in 2:n.iter){
	#update beta
	M.z <- solve(V.inv, crossprod(X,Dat-Eta))
	MCMC$Beta[,iiter] <- as.vector(M.z + solve(chol(V.inv), rnorm(2,0,1)))
	
	#update tau
	MCMC$tau[iiter] <- rgamma(1, n.obs.tau/2 + a.eta, as.numeric(crossprod(Eta, Q %*% Eta)/2) + b.eta)
	
	#Update Eta
	Dat.minus.Exp=Dat-X%*%MCMC$Beta[,iiter]
	V.eta.inv <- diag(S) + tau*Q
	M.eta <- solve(V.eta.inv, Dat.minus.Exp)
	Eta <- M.eta + solve(chol(V.eta.inv), rnorm(S,0,1))
	#center using eq 2.30 of Rue and Held		
	Eta=Eta-mean(Eta)	
}


#Now, the spatially restricted version
P.c=diag(S)-X%*%solve(crossprod(X,X),t(X))
Eigen.P.c=eigen(P.c)
Ind=which(Eigen.P.c$values>0.01)
L.t=Eigen.P.c$vectors[,Ind]
L=t(L.t)
Qt=L%*%Q%*%L.t
cross.L=L%*%L.t

n.theta=nrow(Qt)
Theta=rnorm(n.theta,0,sqrt(1/tau))

for(iiter in 2:n.iter){
	#update beta
	Eta=L.t%*%Theta
	M.z <- solve(V.inv, crossprod(X,Dat-Eta))
	MCMC$Beta[,iiter] <- as.vector(M.z + solve(chol(V.inv), rnorm(2,0,1)))
	
	#update tau
	MCMC$tau[iiter] <- rgamma(1, n.obs.tau/2 + a.eta, as.numeric(crossprod(Theta, Qt %*% Theta)/2) + b.eta)
	
	#Update Theta
	Dat.minus.Exp=Dat-X%*%MCMC$Beta[,iiter]
	V.eta.inv <- cross.L + tau*Qt
	M.eta <- solve(V.eta.inv, L%*%Dat.minus.Exp)
	Theta <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(n.theta,0,1))
	#Theta=Theta-mean(Theta)  should be superfluous!

}


crap<-read.table('Dormann_data.txt',header=TRUE)



	