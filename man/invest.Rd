\name{invest}
\alias{invest}
\alias{invest.lm}
\alias{invest.nls}
\title{Calibration for linear and nonlinear regression.}
\usage{
  invest(x, ...)

  \method{invest}{lm} (object, y0,
    interval = c("inversion", "Wald"), level = 0.95,
    mean.response = FALSE, lower, upper,
    tol = .Machine$double.eps^0.25, maxiter = 1000,
    adjust = c("none", "bonferroni"), k, ...)

  \method{invest}{nls} (object, y0,
    interval = c("inversion", "Wald"), level = 0.95,
    mean.response = FALSE, lower, upper,
    tol = .Machine$double.eps^0.25, maxiter = 1000,
    adjust = c("none", "Bonferroni", "Scheffe"), k, ...)
}
\arguments{
  \item{object}{An appropriate object of class \code{lm} or
  \code{nls}.}

  \item{y0}{The value of the observed response(s) or
  specified value of the mean response.}

  \item{interval}{The type of interval required.}

  \item{level}{A numeric scalar between 0 and 1 giving the
  confidence level for the interval to be calculated.}

  \item{mean.response}{Logical indicating whether
  confidence intervals should correspond to an individual
  response (\code{FALSE}) or a mean response
  (\code{TRUE}).}

  \item{lower}{The lower endpoint of the interval to be
  searched.}

  \item{upper}{The upper endpoint of the interval to be
  searched.}

  \item{tol}{The desired accuracy passed on to
  \code{uniroot}.}

  \item{maxiter}{The maximum number of iterations passed on
  to \code{uniroot}. (\code{TRUE}).}

  \item{adjust}{A logical value indicating if an adjustment
  should be made to the critical value used in calculating
  the confidence interval.  This is useful for when the
  calibration curve is to be used multiple, say k, times.}

  \item{k}{The number times the calibration curve is to be
  used for computing a confidence interval. Only needed
  when \code{adjust = TRUE}.}

  \item{...}{Additional optional arguments. At present no
  optional arguments are used.}
}
\value{
  An object of class \code{calibrate} containing the
  following components: \item{estimate}{The estimate of
  x0.} \item{lwr}{The lower confidence bound on x0.}
  \item{upr}{The upper confidence bound on x0.}
  \item{se}{An estimate of the standard error (Wald
  interval only).} \item{interval}{The method used for
  calculating \code{lower} and \code{upper}.}
}
\description{
  A general function for calibration with linear and
  nonlinear regression models in R.
}
\examples{
data(Puromycin, package = "datasets")
Puromycin2 <- Puromycin[Puromycin$state == "treated", ]
Puro2.nls <- nls(rate ~ (theta1 * conc) / (theta2 + conc),
                 data = Puromycin2, start = list(theta1 = 200, theta2 = 1))
plotFit(Puro2.nls, interval = "both")
invest(Puro2.nls, y0 = 100, interval = "inversion")
invest(Puro2.nls, y0 = 100, interval = "inversion", mean.response = TRUE)
}
\references{
  Huet et. al. (20??). ????
}
\seealso{
  The function \code{\link{calibrate}} for simple linear
  regression models and the function \code{\link{dose.p}}
  in the \code{MASS} package for calibrating bioassays.
}

