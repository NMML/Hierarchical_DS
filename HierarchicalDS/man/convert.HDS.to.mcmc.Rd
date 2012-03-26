\name{convert.HDS.to.mcmc}
\alias{convert.HDS.to.mcmc}
\title{function to convert HierarchicalDS MCMC list vector (used in estimation) into an mcmc object (cf. coda package)}
\usage{
  convert.HDS.to.mcmc(MCMC, N.hab.par, Cov.par.n,
    Hab.names, Det.names, Cov.names, MisID.names,
    N.par.misID = NULL, misID.mat = NULL,
    fix.tau.nu = FALSE, misID = TRUE, spat.ind = TRUE,
    point.ind = TRUE)
}
\arguments{
  \item{MCMC}{list vector providing MCMC samples for each
  parameter type}

  \item{N.hab.par}{see help for mcmc_ds.R}

  \item{Cov.par.n}{see help for mcmc_ds.R}

  \item{Hab.names}{see help for mcmc_ds.R}

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

