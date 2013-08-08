\name{calibrate}
\alias{calibrate}
\alias{calibrate.data.frame}
\alias{calibrate.default}
\alias{calibrate.formula}
\alias{calibrate.list}
\alias{calibrate.lm}
\alias{calibrate.matrix}
\alias{print.calibrate}
\title{Calibration for simple linear regression.}
\usage{
  calibrate(x, ...)

  \method{calibrate}{default} (x, y, y0,
    interval = c("inversion", "Wald"), level = 0.95,
    mean.response = FALSE,
    adjust = c("none", "bonferroni", "scheffe"), k, ...)

  \method{calibrate}{matrix} (z, ...)

  \method{calibrate}{data.frame} (z, ...)

  \method{calibrate}{list} (z, ...)

  \method{calibrate}{formula} (formula, data = NULL, ...,
    subset, na.action = na.fail)

  \method{calibrate}{lm} (object, ...)

  \method{print}{calibrate} (x, digits = 4, ...)
}
\arguments{
  \item{x}{A vector of predictor values.}

  \item{y}{A vector of response values.}

  \item{z}{A matrix , data frame, or list containing the
  data.}

  \item{formula}{A formula of the form \code{y ~ x}.}

  \item{object}{An R object of class \code{lm}.}

  \item{data}{an optional data frame, list or environment
  (or object coercible by \code{as.data.frame} to a data
  frame) containing the variables in the model. If not
  found in data, the variables are taken from
  \code{environment(formula)}, typically the environment
  from which \code{lm} is called.}

  \item{subset}{An optional vector specifying a subset of
  observations to be used in the fitting process.}

  \item{na.action}{a function which indicates what should
  happen when the data contain \code{NA}s. The default is
  set by the \code{na.action} setting of options, and is
  \code{na.fail} if that is unset. The ‘factory-fresh’
  default is \code{na.omit}. Another possible value is
  \code{NULL}, no action. Value \code{na.exclude} can be
  useful.}

  \item{y0}{The value of the observed response(s) or
  specified value of the mean response.}

  \item{interval}{The method to use for forming a
  confidence interval.}

  \item{level}{A numeric scalar between 0 and 1 giving the
  confidence level for the interval to be calculated.}

  \item{mean.response}{Logicial indicating whether
  confidence intervals should correspond to an observed
  response(s) (\code{FALSE}) or a specified value of the
  mean response (\code{TRUE}). Default is \code{FALSE}.}

  \item{adjust}{A logical value indicating if an adjustment
  should be made to the critical value used in calculating
  the confidence interval.  This useful for when the
  calibration curve is to be used multiple, say k, times.}

  \item{k}{The number times the calibration curve is to be
  used for computing a confidence interval. Only needed
  when \code{adjust = TRUE}.}

  \item{...}{Additional optional arguments. At present no
  optional arguments are used.}
}
\value{
  An object of class \code{calibrate} containing the
  following components: \item{estimate}{The maximum
  likelihood estimate of x0.} \item{lwr}{The lower
  confidence bound on x0.} \item{upr}{The upper confidence
  bound on x0.} \item{se}{An estimate of the standard error
  (Wald interval only).} \item{interval}{The method used
  for calculating \code{lower} and \code{upper}.}
}
\description{
  Calibration for simple linear regression models in R.
}
\examples{
\donttest{
## Inverting a prediction interval for an individual response
arsenic.lm <- lm(measured ~ actual, data = arsenic)
plotFit(arsenic.lm, interval = "prediction", shade = TRUE,
        col.pred = "lightblue")
calibrate(arsenic.lm, y0 = 3, interval = "inversion")

## Inverting a confidence interval for the mean response
crystal.lm <- lm(weight ~ time, data = crystal)
plotFit(crystal.lm, interval = "confidence", shade = TRUE,
        col.conf = "lightblue")
calibrate(crystal.lm, y0 = 8, interval = "inversion", mean.response = TRUE)

## Wald interval and approximate standard error based on the delta method
calibrate(crystal.lm, y0 = 8, interval = "Wald", mean.response = TRUE)

## Alterntively, we can use the car package to compute the standard error (this
## is trickier though when mean.respone = FALSE, hence, it is better to use the
## calibrate function).
library(car)
deltaMethod(crystal.lm, g = "(8 - b0) / b1", parameterNames = c("b0", "b1"))
}
}
\references{
  Graybill, F. A., and Iyer, H. K. Regression analysis:
  Concepts and Applications. Belmont, Calif: Duxbury Press,
  1994.
}

