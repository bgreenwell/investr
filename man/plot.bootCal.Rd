% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invest.R
\name{plot.bootCal}
\alias{plot.bootCal}
\title{Plots method for bootstrap calibration}
\usage{
\method{plot}{bootCal}(x, ...)
}
\arguments{
\item{x}{An object that inherits from class \code{"bootCal"}.}

\item{...}{Additional optional arguments. At present, no optional arguments 
are used.}
}
\value{
\code{x} is returned invisibly.
}
\description{
The \code{\link{plot}} method for \code{"bootCal"} objects. In
particular, this method takes a \code{"bootCal"} object and produces plots 
for the bootstrap replicates of the inverse estimate.
}
