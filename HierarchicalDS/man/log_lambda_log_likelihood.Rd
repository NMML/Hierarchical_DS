\name{log_lambda_log_likelihood}
\alias{log_lambda_log_likelihood}
\title{compute the likelihood for nu parameters}
\usage{
  log_lambda_log_likelihood(Log.lambda, DM, Beta, Eta = 0,
    SD, N, Sampled, Area)
}
\arguments{
  \item{Log.lambda}{Log of poisson intensities for total
  areas sampled in each sampled strata}

  \item{DM}{the design matrix}

  \item{Beta}{linear predictor parameters for the log of
  abundance intensity}

  \item{Eta}{a vector of spatial random effects}

  \item{SD}{standard deviation of the overdispersion
  process}

  \item{N}{a vector giving the current iteration's number
  of groups in the area}

  \item{Sampled}{Index for which cells were actually
  sampled}

  \item{Area}{Total area sampled in each sampled cell}
}
\value{
  the log likelihood associated with the data and the
  current set of parameters
}
\description{
  compute the likelihood for nu parameters
}
\author{
  Paul B. Conn
}
\keyword{likelihood}
\keyword{log}

