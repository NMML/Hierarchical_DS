\name{log_lambda_gradient}
\alias{log_lambda_gradient}
\title{compute the first derivative of log_lambda likelihood component for Langevin-Hastings...}
\usage{log_lambda_gradient(index, Mu, Nu, N, var.nu)}
\description{compute the first derivative of log_lambda likelihood component for Langevin-Hastings}
\value{a gradient value}
\keyword{gradient,}
\keyword{Langevin-Hastings}
\author{Paul B. Conn}
\arguments{\item{index}{which element to compute the derivative for (for use with apply functions)}
\item{Mu}{expected value}
\item{Nu}{current observed value}
\item{N}{number of groups in each transect}
\item{var.nu}{variance of the overdispersion process}}
