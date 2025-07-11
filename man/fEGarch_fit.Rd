% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{fEGarch_fit}
\alias{fEGarch_fit}
\title{Generic for Fitting EGARCH Family Models}
\usage{
fEGarch_fit(
  spec,
  rt,
  drange = c(0, 1),
  meanspec = mean_spec(),
  Drange = c(0, 1),
  n_test = 0,
  start_pars = NULL,
  LB = NULL,
  UB = NULL,
  control = list(),
  parallel = TRUE,
  ncores = max(1, future::availableCores() - 1),
  trunc = floor(0.4 * length(rt)),
  presample = 50,
  Prange = c(1, 5)
)
}
\arguments{
\item{spec}{the generic is currently without use.}

\item{rt}{the generic is currently without use.}

\item{drange}{the generic is currently without use.}

\item{meanspec}{the generic is currently without use.}

\item{Drange}{the generic is currently without use.}

\item{n_test}{the generic is currently without use.}

\item{start_pars}{the generic is currently without use.}

\item{LB}{the generic is currently without use.}

\item{UB}{the generic is currently without use.}

\item{control}{the generic is currently without use.}

\item{parallel}{the generic is currently without use.}

\item{ncores}{the generic is currently without use.}

\item{trunc}{the generic is currently without use.}

\item{presample}{the generic is currently without use.}

\item{Prange}{the generic is currently without use.}
}
\value{
The generic is currently without use. Nothing is returned.
}
\description{
The generic is currently without use.
}
