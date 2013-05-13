\name{convert.HDS.to.mcmc}
\alias{convert.HDS.to.mcmc}
\title{function to convert HierarchicalDS MCMC list vector (used in estimation) into an mcmc object (cf. coda package)}
\usage{
  convert.HDS.to.mcmc(MCMC, N.hab.pois.par, N.hab.bern.par,
    Cov.par.n, Hab.pois.names, Hab.bern.names, Det.names,
    Cov.names, MisID.names, N.par.misID = NULL,
    misID.mat = NULL, fix.tau.nu = FALSE, misID = TRUE,
    spat.ind = TRUE, point.ind = TRUE)
}
\arguments{
  \item{MCMC}{list vector providing MCMC samples for each
  parameter type}

  \item{N.hab.pois.par}{see help for mcmc_ds.R}

  \item{N.hab.bern.par}{see help for mcmc_ds.R}

  \item{Cov.par.n}{see help for mcmc_ds.R}

  \item{Hab.pois.names}{see help for mcmc_ds.R}

  \item{Hab.bern.names}{see help for mcmc_ds.R}

  \item{Cov.names}{see help for mcmc_ds.R}

  \item{Det.names}{see help for mcmc_ds.R}

  \item{MisID.names}{see help for mcmc_ds.R}

  \item{N.par.misID}{see help for mcmc_ds.R}

  \item{misID.mat}{see help for mcmc_ds.R}

  \item{misID}{see help for mcmc_ds.R}

  \item{fix.tau.nu}{see help for mcmc_ds.R}

  \item{spat.ind}{see help for mcmc_ds.R}

  \item{point.ind}{see help for mcmc_ds.R}
}
\description{
  function to convert HierarchicalDS MCMC list vector (used
  in estimation) into an mcmc object (cf. coda package)
}
\author{
  Paul B. Conn
}
\keyword{coda}
\keyword{MCMC,}

