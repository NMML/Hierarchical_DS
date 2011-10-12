\name{log_lambda_log_likelihood}
\alias{log_lambda_log_likelihood}
\title{compute the likelihood for nu parameters...}
\usage{log_lambda_log_likelihood(Log.lambda, DM, Beta, Eta=0, SD, N)}
\description{compute the likelihood for nu parameters}
\value{the log likelihood associated with the data and the current set of parameters}
\keyword{log}
\keyword{likelihood}
\author{Paul B. Conn}
\arguments{\item{Log.lambda}{which element to compute the derivative for (for use with apply functions)}
\item{DM}{the design matrix}
\item{Beta}{linear predictor parameters for the log of abundance intensity}
\item{Eta}{a vector of spatial random effects}
\item{SD}{standard deviation of the overdispersion process}
\item{N}{a vector giving the current iteration's number of groups in the area covered by each transect}}
