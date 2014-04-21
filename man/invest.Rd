\name{invest}
\alias{invest}
\alias{invest.lm}
\alias{invest.nls}
\title{Calibration for Linear and Nonlinear Regression Models.}
\usage{
invest(object, ...)

\method{invest}{lm}(object, y0, interval = c("inversion", "Wald"),
  level = 0.95, mean.response = FALSE, lower, upper,
  tol = .Machine$double.eps^0.25, maxiter = 1000, adjust = c("none",
  "Bonferroni"), k, ...)

\method{invest}{nls}(object, y0, interval = c("inversion", "Wald"),
  level = 0.95, mean.response = FALSE, lower, upper,
  tol = .Machine$double.eps^0.25, maxiter = 1000, adjust = c("none",
  "Bonferroni"), k, ...)
}
\arguments{
  \item{object}{An object that inherits from class
  \code{lm} or \code{nls}.}

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
  \code{uniroot}. Recommend a minimum of 1e-10.}

  \item{maxiter}{The maximum number of iterations passed on
  to \code{uniroot}. (\code{TRUE}).}

  \item{adjust}{A logical value indicating if an adjustment
  should be made to the critical value used in calculating
  the confidence interval.  This is useful for when the
  calibration curve is to be used multiple, say k, times.}

  \item{k}{The number times the calibration curve is to be
  used for computing a confidence interval. Only needed
  when \code{adjust = "Bonferroni"}.}

  \item{...}{Additional optional arguments. At present, no
  optional arguments are used.}
}
\value{
An object of class \code{calibrate} containing the
following components: \describe{ \item{\code{estimate}}{The
estimate of x0.} \item{\code{lwr}}{The lower confidence
limit for x0.} \item{\code{upr}}{The upper confidence limit
for x0.} \item{\code{se}}{An estimate of the standard error
(Wald interval only).} \item{\code{interval}}{The method
used for calculating \code{lower} and \code{upper} (only
used by \code{print} method).} }
}
\description{
The function \code{invest} computes the inverse estimate
and a condfidence interval for the unknown predictor value
that corresponds to an observed value of the response (or
vector thereof) or specified value of the mean response.
See the references listed below for more details.
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
Graybill, F. A., and Iyer, H. K. Regression analysis:
Concepts and Applications. Belmont, Calif: Duxbury Press,
1994.

Huet, S., Bouvier, A., Poursat, M-A., and Jolivet, E.
Statistical Tools for Nonlinear Regression: A Practical
Guide with S-PLUS and R Examples. New York: Springer, 2004.

Seber, G. A. F., and Wild, C. J.. Nonlinear regression. New
York: Wiley, 1989.
}

