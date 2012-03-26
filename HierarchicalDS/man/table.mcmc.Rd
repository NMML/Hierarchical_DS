\name{table.mcmc}
\alias{table.mcmc}
\title{function to export posterior summaries from an mcmc object to a table}
\usage{
  \method{table}{mcmc} (MCMC, file = NULL, type = "csv",
    a = 0.05)
}
\arguments{
  \item{MCMC}{An mcmc object with columns referencing
  different parameter types (column names are used for
  plotting labels)}

  \item{file}{A file name to ouput to (including path); if
  null (default), outputs to screen}

  \item{type}{What type of table to produce (either "csv"
  or "tex")}

  \item{a}{Value to use for credible intervals.  For
  example, alpha=0.05 results in 95\% credible intervals}
}
\description{
  function to export posterior summaries from an mcmc
  object to a table
}
\author{
  Paul B. Conn
}
\keyword{MCMC,}
\keyword{table}

