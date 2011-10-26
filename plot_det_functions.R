Beta.det=c(1.2,1.0,0.8,-.2,-.05,.1,0,.2,-.2,-.2,.05,0,-.05,-.1)
X1=c(0,1,0,1,1,4,1,0,0,0,1,0,0,0)
X2=c(0,1,0,2,4,4,1,0,0,0,4,0,0,0)
X3=c(0,1,0,3,9,4,1,0,0,0,9,0,0,0)
X4=c(0,1,0,4,16,4,1,0,0,0,16,0,0,0)
X5=c(0,1,0,5,25,4,1,0,0,0,25,0,0,0)
plot(c(1:5),pnorm(rep(0,5),mean=c(X1%*%Beta.det,X2%*%Beta.det,X3%*%Beta.det,X4%*%Beta.det,X5%*%Beta.det),lower.tail=FALSE),type="l",ylim=c(0,1))
pl=8
pl2=12
for(isp in 2:4){
	X1[pl]=1
	X1[pl-1]=0
	X2[pl]=1
	X2[pl-1]=0
	X3[pl]=1
	X3[pl-1]=0
	X4[pl]=1
	X4[pl-1]=0
	X5[pl]=1
	X5[pl-1]=0
	X1[pl2]=1
	X1[pl2-1]=0
	X2[pl2]=4
	X2[pl2-1]=0
	X3[pl2]=9
	X3[pl2-1]=0
	X4[pl2]=16
	X4[pl2-1]=0
	X5[pl2]=25
	X5[pl2-1]=0
	lines(c(1:5),pnorm(rep(0,5),mean=c(X1%*%Beta.det,X2%*%Beta.det,X3%*%Beta.det,X4%*%Beta.det,X5%*%Beta.det),lower.tail=FALSE),lty=pl-6)
	pl=pl+1
	pl2=pl2+1
}
legend(4,.9,c('Sp1','Sp2','Sp3','Sp4'),lty=c(1:4))


