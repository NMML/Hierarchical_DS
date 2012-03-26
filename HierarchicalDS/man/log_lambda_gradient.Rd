\name{log_lambda_gradient}
\alias{log_lambda_gradient}
\title{compute the first derivative of log_lambda likelihood component for Langevin-Hastings}
\usage{
  log_lambda_gradient(Mu, Nu, Sampled, Area, N, var.nu)
}
\arguments{
  \item{Mu}{expected value for all cells}

  \item{Nu}{current observed valus (all cells)}

  \item{Sampled}{Vector giving the cell identities for all
  sampled cells}

  \item{Area}{Proportional area of each sampled cell that
  is covered by one or more transects}

  \item{N}{number of groups in each transect}

  \item{var.nu}{variance of the overdispersion process}
}
\value{
  a gradient value
}
\description{
  compute the first derivative of log_lambda likelihood
  component for Langevin-Hastings
}
\author{
  Paul B. Conn
}
\keyword{gradient,}
\keyword{Langevin-Hastings}

