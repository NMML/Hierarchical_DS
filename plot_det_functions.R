Beta.det=c(1,-.4,-.5,-1.2,.05,-.2,.2)
X1=c(1,0,0,0,4,0,0)
X2=c(1,0,.2,.04,4,0,0)
X3=c(1,0,.4,.16,4,0,0)
X4=c(1,0,.6,.36,4,0,0)
X5=c(1,0,.8,.64,4,0,0)  #estimated detection probability
X6=c(1,0,1,1,4,0,0)  #estimated detection probability

plot(c(0,.2,.4,.6,.8,1),pnorm(rep(0,6),mean=c(X1%*%Beta.det,X2%*%Beta.det,X3%*%Beta.det,X4%*%Beta.det,X5%*%Beta.det,X6%*%Beta.det),lower.tail=FALSE),type="l",ylim=c(0,1))

Beta.det=c(1.2,1.0,0.8,-.6,-.6,-.9,-1,0,.3)
X1=c(0,1,0,0,0,0,0,4,.5)
X2=c(0,1,0,1,0,0,0,4,.5)
X3=c(0,1,0,0,1,0,0,4,.5)
X4=c(0,1,0,0,0,1,0,4,.5)
X5=c(0,1,0,0,0,0,1,4,.5)  #estimated detection probability
X1=c(0,1,0,0,0,0,0,4,.5)
X2=c(0,1,0,1,0,0,0,4,.5)
X3=c(0,1,0,0,1,0,0,4,.5)
X4=c(0,1,0,0,0,1,0,4,.5)
X5=c(0,1,0,0,0,0,1,4,.5)  #estimated detection probability

lines(c(1:5),pnorm(rep(0,5),mean=c(X1%*%Beta.det,X2%*%Beta.det,X3%*%Beta.det,X4%*%Beta.det,X5%*%Beta.det),lower.tail=FALSE),lty=2)


legend(4,.9,c('True','Estimated'),lty=c(1:2))


