#' script to process golf tee data and run MCMC estimation 
#' @return see help for hierarchicalDS
#' @export
#' @keywords golf tee data
#' @author Paul B. Conn

Tee.data<-read.table("c:/users/paul.conn/git/hierarchicalDS/allteesdata.txt")
plot.data=TRUE
if(plot.data){
	setwd("/Users/paul.conn/git/hierarchicalDS/")
	gex.data=read.csv("gex.data.csv")
	attach(gex.data)
	I.obs=(rowSums(gex.data[,16:23])>0)
	plot(east,north,cex=(s+3)/5,col=gray(1-(c+.5)/max(c+.5)),pch=17-2*v-15*(1-I.obs),xaxt='n',yaxt='n',xlab='',ylab='') #s is group size, v is visibility (0/1), c is color (0/1=green/yellow, referred to as "sex" in the book)
	detach(gex.data)
	glines=read.csv("glines.csv")
	attach(glines)
	for(i in 1:(length(glines[,1])/2)) lines(c(x[2*i-1],x[2*i]),c(y[2*i-1],y[2*i]),col='black',lwd=2)
	for(i in 1:(length(glines[,1])/2)) lines(c(x[2*i-1]-4,x[2*i]-4),c(y[2*i-1],y[2*i]),col='black',lty=3,lwd=2)
	for(i in 1:(length(glines[,1])/2)) lines(c(x[2*i-1]+4,x[2*i]+4),c(y[2*i-1],y[2*i]),col='black',lty=3,lwd=2)
	gbnd=read.csv("gbnd.csv")
	detach(glines)
	attach(gbnd)
	for(i in 1:(length(gbnd[,1])-1)){
		if(i!=15)lines(c(x[i],x[i+1]),c(y[i],y[i+1]),col="black",lty=3,lwd=2)
		#if(i !=15)lines(c(x[i],x[i+1]),c(y[i],y[i+1]),col=type[i]*Stratum[i]+2,lwd=2)
		else lines(c(x[i],x[i+1]),c(y[i],y[i+1]),col="red",lwd=2)
	}
	text(40,25,"Strata 1",cex=1.4)
	text(40,-17,"Strata 2",cex=1.4)
}


set.seed(12345) #chain 1
set.seed(12346)

#1) format data for hierarchicalDS; remove unobserved tees
#note observations o10 correspond to observer group 1 and o11 corresponds to observer group 2
n.transects=11
n.species=2
Observed=Tee.data[,"o10"]+Tee.data[,"o11"]
Observed=(Observed>0)
Dat<-as.matrix(Tee.data[Observed,])
n.obs=sum(Observed)
Dat2<-matrix(0,n.obs*2,8)
for(i in 1:n.obs){
	Dat2[i*2-1,]=c(Dat[i,c("Strip","Animal")],1,Dat[i,"o10"],Dat[i,"sex"],abs(Dat[i,"SignPerp"]),Dat[i,c("size","exposure")])
	Dat2[i*2,]=c(Dat[i,c("Strip","Animal")],2,Dat[i,"o11"],Dat[i,"sex"],abs(Dat[i,"SignPerp"]),Dat[i,c("size","exposure")])
}
colnames(Dat2)=c("Transect","Match","Observer","Observation","Species","Distance","Group","Exposed")
Dat<-data.frame(Dat2)
Dat[,"Distance"]=Dat[,"Distance"]/4  #standardize so in (0,1) range
Dat[,"Species"]=as.factor(Dat[,"Species"]+1)  #by default, all parameters with a dirichlet prior are drawn at integer values starting at 1 (not 0!)
Dat[,"Exposed"]=as.factor(Dat[,"Exposed"]+1)


#2) Setup some other inputs
Observers=matrix(rep(c(1,2),n.transects),2,n.transects)
Det.formula=~Observer+Distance+Group+Species*Exposed+Group:Species
Cov.prior.pdf=c("pois1_ln","multinom","multinom")

Cov.prior.fixed=rep(FALSE,3)
Cov.prior.parms=array(0,dim=c(n.species,3,2))
Cov.prior.parms[1,,]=matrix(c(1.1,1,1,1,1,0),3,2)
Cov.prior.parms[2,,]=matrix(c(1.1,1,1,1,1,0),3,2)
#Cov.prior.parms=matrix(c(3,0,0,.5,.3,.2),3,2)
Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
Cov.prior.pdf=Cov.prior.fixed
Cov.prior.pdf[,1]=c("pois1_ln","pois1_ln")
Cov.prior.pdf[,2]=rep("multinom",2)
Cov.prior.n=c(3,2)

Control=list(iter=100100,burnin=100,thin=20,MH.cor=0.2,MH.nu=c(.01,.01),RJ.N=matrix(rep(3,n.transects*n.species),n.species,n.transects),adapt=1000)

#Out=hierarchicalDS(Dat=Dat,Adj=1,Area.hab=1,Mapping=1,Area.trans=1,Observers=Observers,n.obs.cov=0,Hab.formula=~1,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,pol.eff=c(1:3),point.ind=TRUE,spat.ind=TRUE,grps=TRUE,M=500,Control=Control,adapt=TRUE,Prior.pars=Prior.pars)


#for debugging, define everything beforehand
Adj=diag(n.transects)
Mapping=c(1:n.transects)
Area.trans=rep(1,n.transects)
Area.hab=c(10,30,30,27,21,12,23,23,15,12,7)
Area.hab=Area.hab/mean(Area.hab)
Bin.length=NULL
Hab.cov=data.frame(Strata=c(0,0,0,0,0,0,1,1,1,1,1)) #strata
Obs.cov=NULL
Hab.formula=c(~Strata,~Strata)
misID=NULL
misID.mat=NULL
misID.models=NULL
n.obs.cov=0
pol.eff=c(1:3)
point.ind=TRUE
spat.ind=TRUE
fix.tau.nu=FALSE
srr=FALSE
srr.tol=0.5
misID=FALSE
Inits=NULL
Levels=NULL
grps=TRUE
M=matrix(50,n.species,n.transects)
adapt=TRUE
Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.sd=c(10000,100)) #(1,.01) prior makes it closer to a uniform distribution near the origin
Inits=list(tau.nu=c(100,100),hab=matrix(c(2,2.4,0,0),2,2)) #chain1
adapt=TRUE
Out=hierarchicalDS(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,Observers=Observers,Bin.length=Bin.length,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.formula=Hab.formula,Det.formula=Det.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,pol.eff=NULL,point.ind=TRUE,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,misID=misID,Inits=Inits,grps=grps,M=M,Control=Control,Levels=Levels,adapt=TRUE,Prior.pars=Prior.pars,misID.mat=misID.mat,misID.models=misID.models)
#save(Out,file="c:/users/paul.conn/git/hierarchicalDS/Output/golf_tee2.Rdat")

#pull in data and plot
load("c:/users/paul.conn/git/hierarchicalDS/Output/golf_tee1.Rdat")
Out1=Out
load("c:/users/paul.conn/git/hierarchicalDS/Output/golf_tee2.Rdat")
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
burnin=10000
thin=10
i.first=burnin/thin+1
i.last=length(Out$MCMC$N.tot[1,])
N.tot.1=c(Out$MCMC$N.tot[1,i.first:i.last])
N.tot.2=c(Out$MCMC$N.tot[2,i.first:i.last])
N.tot=N.tot.1+N.tot.2
Hab.strata.sp1=c(Out$MCMC$Hab[1,i.first:i.last,2],Out1$MCMC$Hab[1,i.first:i.last,2])
Hab.strata.sp2=c(Out$MCMC$Hab[2,i.first:i.last,2],Out1$MCMC$Hab[2,i.first:i.last,2])
G.tot.col0=rowSums(rbind(Out$MCMC$G[1,i.first:i.last,],Out1$MCMC$G[1,i.first:i.last,]))
G.tot.col1=rowSums(rbind(Out$MCMC$G[2,i.first:i.last,],Out1$MCMC$G[2,i.first:i.last,]))
G.tot=G.tot.col0+G.tot.col1
Prop.exp0.sp1=c(Out$MCMC$Cov.par[1,i.first:i.last,4],Out1$MCMC$Cov.par[1,i.first:i.last,4]) #estimated proportion exposed
Prop.exp0.sp2=c(Out$MCMC$Cov.par[2,i.first:i.last,4],Out1$MCMC$Cov.par[2,i.first:i.last,4]) #estimated proportion exposed

G.exp0.sp1=G.tot.col0*Prop.exp0.sp1
G.exp1.sp1=G.tot.col0*(1-Prop.exp0.sp1)
G.exp0.sp2=G.tot.col1*Prop.exp0.sp2
G.exp1.sp2=G.tot.col1*(1-Prop.exp0.sp2)

#plot posteriors vs. truth vs. mrds point estimates from Laake et al. '04
G.mat=cbind(G.exp0.sp1,G.exp1.sp1,G.exp0.sp2,G.exp1.sp2)

colnames(G.mat)=c("A. Green, not exposed","B. Green, exposed","C. Yellow, not exposed","D. Yellow, exposed")

n=nrow(G.mat)
Plot.mat=data.frame(matrix(0,n,2))
Plot.mat[,1]=as.character(Plot.mat[,1])
colnames(Plot.mat)=c("Name","Posterior")
Plot.mat[,1]=colnames(G.mat)[1]
Plot.mat[,2]=G.mat[,1]
for(i in 2:ncol(G.mat)){
	Plot.mat[((i-1)*n+1):(i*n),1]=colnames(G.mat)[i]
	Plot.mat[((i-1)*n+1):(i*n),2]=G.mat[,i]
}

V.lines1=data.frame(matrix(0,length(colnames(G.mat)),2))
V.lines1[,1]=colnames(G.mat)
V.lines1[,2]=c(61,47,76,66)  
colnames(V.lines1)=c("Name","truth")

V.lines2=data.frame(matrix(0,length(colnames(G.mat)),2))
V.lines2[,1]=colnames(G.mat)
V.lines2[,2]=c(42,74,58,78)  
colnames(V.lines2)=c("Name","mrds")

qplot(Posterior,data=Plot.mat,geom="density",xaxt='n')+geom_vline(aes(xintercept=truth),data=V.lines1,colour="red")+geom_vline(aes(xintercept=mrds),data=V.lines2,colour="blue")+scale_y_continuous("Density")+opts(axis.text.x=theme_blank(),axis.text.y=theme_blank())+facet_wrap(~Name,scale="free")


### plot predicted, empirical group sizes

#hist(Tee.data[,"size"])
Cov.par.sp1=rbind(Out$MCMC$Cov.par[1,i.first:i.last,],Out1$MCMC$Cov.par[1,i.first:i.last,])
Cov.par.means.sp1=apply(Cov.par.sp1,2,'mean')
Cov.par.sp2=rbind(Out$MCMC$Cov.par[2,i.first:i.last,],Out1$MCMC$Cov.par[2,i.first:i.last,])
Cov.par.means.sp2=apply(Cov.par.sp2,2,'mean')

#Dist.sp1=1+rpois(10000,exp(Cov.par.means.sp1[1]+Cov.par.means.sp1[2]*rnorm(10000,0,1)))
#hist(Dist.sp1,freq=0)

Size.df=data.frame(matrix(0,20000+nrow(Tee.data),2))
Size.df[,1]=as.character(Size.df[,1])
Size.df[1:108,1]="Empirical: green"
Size.df[109:10108,1]="Predicted: green"
Size.df[10109:10250,1]="Empirical: yellow"
Size.df[10251:20250,1]="Predicted: yellow"
colnames(Size.df)=c("type","value")
Size.df[1:108,2]=Tee.data[which(Tee.data[,"sex"]==0),"size"]
Size.df[109:10108,2]=1+rpois(10000,exp(Cov.par.means.sp1[1]+Cov.par.means.sp1[2]*rnorm(10000,0,1)))
Size.df[10109:10250,2]=Tee.data[which(Tee.data[,"sex"]==1),"size"]
Size.df[10251:20250,2]=1+rpois(10000,exp(Cov.par.means.sp2[1]+Cov.par.means.sp2[2]*rnorm(10000,0,1)))
Size.df=Size.df[Size.df[,2]<=12,]

qplot(value,data=Size.df,geom="bar",binwidth=1)+aes(y = ..density..)+
		scale_y_continuous("Density")+scale_x_continuous("Group size")+
		#opts(axis.text.y=theme_blank())+
		facet_grid(type~.)

set.seed(12345)
mass.emp.green=tabulate(Tee.data[which(Tee.data[,"sex"]==0),"size"])
mass.emp.green=mass.emp.green/sum(mass.emp.green)
mass.emp.yel=tabulate(Tee.data[which(Tee.data[,"sex"]==1),"size"])
mass.emp.yel=mass.emp.yel/sum(mass.emp.yel)
mass.post.green=tabulate(1+rpois(10000,exp(Cov.par.means.sp1[1]+Cov.par.means.sp1[2]*rnorm(10000,0,1))))
mass.post.green=mass.post.green[1:12]
mass.post.green=mass.post.green/sum(mass.post.green)
mass.post.yel=tabulate(1+rpois(10000,exp(Cov.par.means.sp2[1]+Cov.par.means.sp2[2]*rnorm(10000,0,1))))
mass.post.yel=mass.post.yel/sum(mass.post.yel)
par(mfrow=c(2,1),mar=c(3,4,2,2),cex.axis=1.2,cex.lab=1.2)
mass.green=matrix(0,2,12)
mass.yellow=mass.green
mass.green[1,1:8]=mass.emp.green
mass.green[2,]=mass.post.green
mass.yellow[1,1:8]=mass.emp.yel
mass.yellow[2,]=mass.post.yel
barplot(mass.green,beside=T,col=c("darkgreen","lightgreen"),ylab="Mass",legend.text=c("Green empirical","Green posterior"),args.legend=list(cex=1.2))
par(mar=c(5,4,0,2))
barplot(mass.yellow,beside=T,col=c("yellow3","lightyellow"),ylab="Mass",xlab="Group size",names.arg=c(1:12),legend.text=c("Yellow empirical","Yellow posterior"),args.legend=list(cex=1.2))

#plot  observer-specific detection functions,
#conditional detection functions, delta dependence function, duplicate detection function (seen by both),
#and pooled detection function (seen by at least one).
#' @param x vector of perpendicular distances
#' @param formula linear probit formula for detection using distance and other covariates
#' @param beta parameter values
#' @param rho maximum correlation at largest distance
Det=rbind(Out$MCMC$Det[i.first:i.last,],Out1$MCMC$Det[i.first:i.last,])
Rho=rbind(Out$MCMC$cor[i.first:i.last],Out1$MCMC$cor[i.first:i.last])
test=probit.fct(x=c(0:100)/100,formula=~observer+distance+Group+Color*Exposed+Group:Color,beta=colMeans(Det),rho=mean(Rho),Group=c(1,4,7),Color=c(1:2),Exposed=c(1:2))
layout(matrix(c(1,2)),heights=c(lcm(8),lcm(8)),widths=c(lcm(10),1))
with(test[test$observer==2 & test$Group==1 & test$Color==1 & test$Exposed==1,],
{
par(mar=c(3,4,0,2))
plot(distance,p,ylim=c(0,1),xlab='',xaxt='n',ylab="Detection probability",type="l",lwd=2)
lines(distance,pc,lty=2,lwd=2)
lines(distance,dup,lty=4,lwd=2)
lines(distance,pool,lty=3,lwd=2)
legend(.01,.4,legend=c("Individual","Conditional","Duplicate","Pooled"),lty=c(1,2,4,3),lwd=rep(2,4))
par(mar=c(4,4,0,2))
plot(distance,delta,xlab="Distance",ylab="Dependence",type="l",lwd=2)
})


