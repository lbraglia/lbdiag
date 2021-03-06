% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutoff_da.R
\name{cutoff_da}
\alias{cutoff_da}
\title{Calculate diagnostic accuracy measures for several cutoffs of a quantitative
marker.}
\usage{
cutoff_da(
  cutoffs = NULL,
  direction = c(">", "<"),
  test = NULL,
  refstd = NULL,
  digits = 4,
  ...
)
}
\arguments{
\item{cutoffs}{Cutoffs considered}

\item{direction}{> or <: if > a test value greater than the cutoff
will be interpreted as a positive test, in < a test value below
the threshold}

\item{test}{Test}

\item{refstd}{Reference standard}

\item{digits}{Rounding digits}

\item{...}{parameters passed to da}
}
\value{
A data.frame for diagnostic accuracy studies.
}
\description{
Calculate diagnostic accuracy measures for several cutoffs of a quantitative
marker.
}
