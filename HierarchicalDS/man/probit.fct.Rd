\name{probit.fct}
\alias{probit.fct}
\title{Mrds probit detection and related functions}
\usage{
  probit.fct(x, formula, beta, rho, ...)
}
\arguments{
  \item{x}{vector of perpendicular distances}

  \item{formula}{linear probit formula for detection using
  distance and other covariates}

  \item{beta}{parameter values}

  \item{rho}{maximum correlation at largest distance}

  \item{...}{any number of named vectors of covariates used
  in the formula}
}
\value{
  dat dataframe with distance, observer, any covariates
  specified in ... and detection probability p, conditional
  detection probability pc, dupiicate detection dup, pooled
  detection pool and dependence pc/p=delta.
}
\description{
  For independent observers, probit.fct computes
  observer-specific detection functions, conditional
  detection functions, delta dependence function, duplicate
  detection function (seen by both), and pooled detection
  function (seen by at least one).
}
\details{
  The vectors of covariate values can be of different
  lengths because expand.grid is used to create a dataframe
  of all unique combinations of the distances and covariate
  values and the detection and related values are computed
  for each combination.  The covariate vector observer=1:2
  is automatically included.
}
\examples{
test=probit.fct(0:10,~distance,c(1,-.15),.8,size=1:3)
par(mfrow=c(1,2))
with(test[test$observer==1,],
{plot(distance,p,ylim=c(0,1),xlab="Distance",ylab="Detection probability")
points(distance,pc,pch=2)
points(distance,dup,pch=3)
points(distance,pool,pch=4)
legend(1,.2,legend=c("Detection","Conditional detection","Duplicate detection","Pooled detection"),pch=1:4,bty="n")
plot(distance,delta,xlab="Distance",ylab="Dependence")
})
}
\author{
  Jeff Laake
}

