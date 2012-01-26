# Analyze simulation data for PLoS One paper
#
# Paul Conn 1/25/12
rm(list=ls())
setwd("c:/users/paul.conn/git/Hierarchical_DS/Output/sims")
load("Out1a.Rdata")
Out1=Out
load("Out1b.Rdata")

N.tot=c(Out1$MCMC$N.tot,Out$MCMC$N.tot)
N=rbind(Out1$MCMC$N,Out$MCMC$N)

#plot estimated abundance
S=ncol(Out$MCMC$N)
N.mean=apply(N,2,'mean')
Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(N.mean))))
colnames(Abund.df)=c("y","x","Abundance")
require(ggplot2)
crap<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
crap

###posteriors vs. truth
#par(mfrow=c(5,5))
Hab=rbind(Out$MCMC$Hab,Out1$MCMC$Hab)
Det=rbind(Out$MCMC$Det,Out1$MCMC$Det)
cor=c(Out$MCMC$cor,Out1$MCMC$cor)
tau=c(Out$MCMC$tau.eta,Out1$MCMC$tau.eta)
Cov.par=rbind(Out$MCMC$Cov.par,Out1$MCMC$Cov.par)
Posts=cbind(N.tot,Hab,Det,cor,tau,Cov.par[,c(1,4,5)])
colnames(Posts)=c("N","Hab.0","Hab.lat","Det.0","Det.Obs2","Det.Obs3","Det.Dist2","Det.Dist3","Det.Dist4","Det.Dist5","Det.Group","Det.Sp2","Det.Sp3","cor","tau.eta","Grp.size","Sp.prop1","Sp.prop2")

n=length(N.tot)
Plot.mat=data.frame(matrix(0,n,2))
Plot.mat[,1]=as.character(Plot.mat[,1])
colnames(Plot.mat)=c("Name","Posterior")
Plot.mat[,1]=colnames(Posts)[1]
Plot.mat[,2]=Posts[,1]
for(i in 2:ncol(Posts)){
	Plot.mat[((i-1)*n+1):(i*n),1]=colnames(Posts)[i]
	Plot.mat[((i-1)*n+1):(i*n),2]=Posts[,i]
}

V.lines=data.frame(matrix(0,length(colnames(Posts)),2))
V.lines[,1]=colnames(Posts)
V.lines[,2]=c(95978,log(200),1,1.2,-.2,-.4,-.8,-1.4,-1.8,-2,.1,.2,-.4,0.5,0.15,3,0.5,0.3)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
colnames(V.lines)=c("Name","truth")

qplot(Posterior,data=Plot.mat,geom="density",xaxt='n')+geom_vline(aes(xintercept=truth),data=V.lines,colour="red")+scale_y_continuous("Density")+opts(axis.text.x=theme_blank(),axis.text.y=theme_blank())+facet_wrap(~Name,scale="free")




load("Outcluster1a.Rdata")
Out1=Out
load("Outcluster1b.Rdata")

N.tot=c(Out1$MCMC$N.tot[501:2500],Out$MCMC$N.tot[501:2500])
N=rbind(Out1$MCMC$N[501:2500,],Out$MCMC$N[501:2500,])
Hab=c(as.vector(Out1$MCMC$Hab[501:2500,]),as.vector(Out$MCMC$Hab[501:2500,]))

#plot estimated abundance
S=ncol(Out$MCMC$N)
N.mean=apply(N,2,'mean')
Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(N.mean))))
colnames(Abund.df)=c("y","x","Abundance")
require(ggplot2)
crap<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
crap

#plot posteriors vs. truth
Det=rbind(Out$MCMC$Det[501:2500,],Out1$MCMC$Det[501:2500,])
cor=c(Out$MCMC$cor[501:2500],Out1$MCMC$cor[501:2500])
Cov.par=rbind(Out$MCMC$Cov.par[501:2500,],Out1$MCMC$Cov.par[501:2500,])
Posts=cbind(N.tot,Det,cor,Cov.par[,c(1,4,5)])

colnames(Posts)=c("N","Det.0","Det.Obs2","Det.Obs3","Det.Dist2","Det.Dist3","Det.Dist4","Det.Dist5","Det.Group","Det.Sp2","Det.Sp3","cor","Grp.size","Sp.prop1","Sp.prop2")

n=length(N.tot)
Plot.mat=data.frame(matrix(0,n,2))
Plot.mat[,1]=as.character(Plot.mat[,1])
colnames(Plot.mat)=c("Name","Posterior")
Plot.mat[,1]=colnames(Posts)[1]
Plot.mat[,2]=Posts[,1]
for(i in 2:ncol(Posts)){
	Plot.mat[((i-1)*n+1):(i*n),1]=colnames(Posts)[i]
	Plot.mat[((i-1)*n+1):(i*n),2]=Posts[,i]
}

V.lines=data.frame(matrix(0,length(colnames(Posts)),2))
V.lines[,1]=colnames(Posts)
V.lines[,2]=c(48649,1.5,-.2,-.5,-.4,-.7,-.9,-1,.1,.2,-.4,0.5,3,0.5,0.3)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
colnames(V.lines)=c("Name","truth")

qplot(Posterior,data=Plot.mat,geom="density",xaxt='n')+geom_vline(aes(xintercept=truth),data=V.lines,colour="red")+scale_y_continuous("Density")+opts(axis.text.x=theme_blank(),axis.text.y=theme_blank())+facet_wrap(~Name,scale="free")



