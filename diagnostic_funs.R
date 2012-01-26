#diagnostic functions

plot_det_distance<-function(Cur.X,Det){
	Cur.Y.tilde=rep(0,5)
	Cur.X[4:7]=0
	Cur.Y.tilde[1]=Cur.X%*%Det
	Cur.X[4]=1
	Cur.Y.tilde[2]=Cur.X%*%Det	
	for(i in 3:5){
		Cur.X[i+1]=0
		Cur.X[i+2]=1
		Cur.Y.tilde[i]=Cur.X%*%Det	
	}
	plot(c(1:5),pnorm(rep(0,5),Cur.Y.tilde,rep(1,5),lower.tail=FALSE),type="l",ylim=c(0,1),lwd=3)
}


