# Analyze simulation data for PLoS One paper
#
# Paul Conn 1/25/12
rm(list=ls())
require(ggplot2)
setwd("c:/users/paul.conn/git/Hierarchical_DS/Output/sims")
load("Out_two_sp1.Rdata")
Out1=Out
load("Out_two_sp2.Rdata")
par(mfrow=c(2,2))
plot(Out$MCMC$N.tot[1,],type="l",col=1)
lines(Out1$MCMC$N.tot[1,],type="l",col=2)
plot(Out$MCMC$N.tot[2,],type="l",col=1)
lines(Out1$MCMC$N.tot[2,],type="l",col=2)
plot(Out$MCMC$Hab[1,,1],type="l",col=1)
lines(Out1$MCMC$Hab[1,,1],type="l",col=2)
plot(Out$MCMC$Det[,1],type="l",col=1)
lines(Out1$MCMC$Det[,1],type="l",col=2)


#no evidence of a burnin really needed; however, remove first 10,000 iterations from each chain
burnin=20000
thin=20
i.first=burnin/thin+1
i.last=length(Out$MCMC$N.tot[1,])
N.tot.1=c(Out$MCMC$N.tot[1,i.first:i.last],Out1$MCMC$N.tot[1,i.first:i.last])
N.tot.2=c(Out$MCMC$N.tot[2,i.first:i.last],Out1$MCMC$N.tot[2,i.first:i.last])
N.tot=N.tot.1+N.tot.2
Hab.sp1=rbind(Out$MCMC$Hab[1,i.first:i.last,],Out1$MCMC$Hab[1,i.first:i.last,])
Hab.sp2=rbind(Out$MCMC$Hab[2,i.first:i.last,],Out1$MCMC$Hab[2,i.first:i.last,])
G.sp1=rbind(Out$MCMC$G[1,i.first:i.last,],Out1$MCMC$G[1,i.first:i.last,])
G.sp2=rbind(Out$MCMC$G[2,i.first:i.last,],Out1$MCMC$G[2,i.first:i.last,])
G.tot.sp1=rowSums(rbind(Out$MCMC$G[1,i.first:i.last,],Out1$MCMC$G[1,i.first:i.last,]))
G.tot.sp2=rowSums(rbind(Out$MCMC$G[2,i.first:i.last,],Out1$MCMC$G[2,i.first:i.last,]))
G.tot=G.tot.sp1+G.tot.sp2

#true abundance
S=ncol(Out$MCMC$G[1,,])
X.site1=cbind(rep(1,S),log(c(1:S)/S),(log(c(1:S)/S))^2) #covariate on abundance intensity, sp 1
X.site2=X.site1 #covariate on abundance intensity, sp 1
Beta.site1=c(log(40),1,0) 
Beta.site2=c(log(10),-2,-1)
N1=round(exp(X.site1%*%Beta.site1))
N2=round(exp(X.site2%*%Beta.site2))

#plot estimated and true abundance
S=ncol(G.sp1)
library(gplots)
layout(matrix(c(1,2)),heights=c(lcm(8),lcm(8)),widths=c(lcm(10),1))
par(mar=c(3,4,1,2))
plotCI(gap=0,x=c(1:S),y=colMeans(G.sp1),li=apply(G.sp1,2,'quantile',0.05),xlab='',ui=apply(G.sp1,2,'quantile',0.95),xaxt='n',ylab="Abundance (# of groups)",cex.lab=1.2,cex.axis=1.2,,xlim=c(-.01,25))
lines(c(1:S),N1,lwd=2,col=2)
text(2,41,"A.",cex=1.4)
par(mar=c(4,4,0,2))
plotCI(gap=0,x=c(1:S),y=colMeans(G.sp2),li=apply(G.sp2,2,'quantile',0.05),ui=apply(G.sp2,2,'quantile',0.95),xlab='Transect',ylab="Abundance (# of groups)",cex.lab=1.2,cex.axis=1.2,xlim=c(-.01,25))
lines(c(1:S),N2,lwd=2,col=2)
text(2,29,"B.",cex=1.4)


###posteriors vs. truth
#par(mfrow=c(5,5))
Hab=cbind(Hab.sp1[,1:2],Hab.sp2)
Det=rbind(Out$MCMC$Det[i.first:i.last,],Out1$MCMC$Det[i.first:i.last,])
cor=c(Out$MCMC$cor[i.first:i.last],Out1$MCMC$cor[i.first:i.last])
tau=c(Out$MCMC$tau.nu[i.first:i.last],Out1$MCMC$tau.nu[i.first:i.last])
Cov.par1=c(Out$MCMC$Cov.par[1,i.first:i.last,1],Out1$MCMC$Cov.par[1,i.first:i.last,1])
Cov.par2=c(Out$MCMC$Cov.par[2,i.first:i.last,1],Out1$MCMC$Cov.par[2,i.first:i.last,1])

#do some cleaning up!
#rm(Out,Out1,N1,N2,Beta.site1,Beta.site2,G.sp1,G.sp2,G.tot,G.tot.sp1,N.tot,X.site1,X.site2)

Posts=cbind(N.tot.1,N.tot.2,Hab,Det,cor,Cov.par1,Cov.par2)
colnames(Posts)=c("N.sp1","N.sp2","Hab.Interc.sp1","Hab.Lin.sp1","Hab.Interc.sp2","Hab.Lin.sp2","Hab.Quad.sp2",
		"Det.Interc","Det.Obs2","Det.Obs3","Det.Dist2","Det.Dist3","Det.Dist4","Det.Dist5","Det.Group","Det.Sp2","cor",
		"Cov.mu.size.sp1","Cov.mu.size.sp2")

n=length(N.tot.1)
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
V.lines[,2]=c(2062,873,log(40),1,log(10),-2,-1,1.2,-.2,-.4,-.6,-.9,-1.1,-1.3,.1,.3,0.5,3,1)  #obs 1 (bin 1), obs 2, obs 3, offset for bin 2, ..., offset for bin n.bins, grp size,species
colnames(V.lines)=c("Name","truth")

qplot(Posterior,data=Plot.mat,geom="density",xaxt='n')+geom_vline(aes(xintercept=truth),data=V.lines,colour="red")+scale_y_continuous("Density")+opts(axis.text.x=theme_blank(),axis.text.y=theme_blank())+facet_wrap(~Name,scale="free")


#Plot realized versus true detection probability curves
Cov=matrix(0,5,9)
Cov[1,]=c(1,1,0,0,0,0,0,4,0)
for(icov in 2:5){
	Cov[icov,]=Cov[1,]
	Cov[icov,icov+2]=1
}
Cov2=matrix(0,5,9)
Cov2[1,]=c(1,1,0,0,0,0,0,2,1)
for(icov in 2:5){
	Cov2[icov,]=Cov2[1,]
	Cov2[icov,icov+2]=1
}

Beta.post=colMeans(Posts[,8:16])
Beta.true=c(1.2,-.2,-.4,-.6,-.9,-1.1,-1.3,.1,.3)
Preds.post=pnorm(Cov%*%Beta.post,0,1)
Preds.true=pnorm(Cov%*%Beta.true,0,1)
plot(c(1:5),Preds.true,type="l",lwd=2,ylim=c(0,1),cex.axis=1.3,cex.lab=1.3,xlab="Distance bin",ylab="Detection probability")
lines(c(1:5),Preds.post,lty=2,lwd=2)
Preds.post2=pnorm(Cov2%*%Beta.post,0,1)
Preds.true2=pnorm(Cov2%*%Beta.true,0,1)
lines(c(1:5),Preds.true2,lty=3,lwd=2)
lines(c(1:5),Preds.post2,lty=4,lwd=2)

legend(1.5,.4,c("True, species 1","Posterior, species 1","True, species 2","Posterior, species 2"),lty=c(1:4),lwd=rep(2,4),cex=1.3)


########## spatial sims only!

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



