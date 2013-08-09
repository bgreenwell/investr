\name{plotFit}
\alias{plotFit}
\alias{plotFit.lm}
\alias{plotFit.nls}
\title{Plotting Confidence/Prediction Bands}
\usage{
  plotFit(object, ...)

  \method{plotFit}{lm} (object,
    interval = c("none", "both", "confidence", "prediction"),
    level = 0.95,
    adjust = c("none", "Bonferroni", "Scheffe", "W-H"), k,
    ..., shade = FALSE, extend.range = FALSE, hide = TRUE,
    col.conf = if (shade) grey(0.7) else "black",
    col.pred = if (shade) grey(0.9) else "black",
    border.conf = col.conf, border.pred = col.pred,
    col.fit = "black", lty.conf = if (shade) 1 else 2,
    lty.pred = if (shade) 1 else 3, lty.fit = 1,
    lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, xlab,
    ylab, xlim, ylim = NULL)

  \method{plotFit}{nls} (object,
    interval = c("none", "both", "confidence", "prediction"),
    level = 0.95,
    adjust = c("none", "Bonferroni", "Scheffe", "W-H"), k,
    ..., shade = FALSE, extend.range = FALSE, hide = TRUE,
    col.conf = if (shade) grey(0.7) else "black",
    col.pred = if (shade) grey(0.9) else "black",
    border.conf = col.conf, border.pred = col.pred,
    col.fit = "black", lty.conf = if (shade) 1 else 2,
    lty.pred = if (shade) 1 else 3, lty.fit = 1,
    lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, xlab,
    ylab, xlim, ylim = NULL)
}
\arguments{
  \item{object}{An object that inherits from class
  \code{lm} or \code{nls}.}

  \item{interval}{A character string indicating if a
  prediction band, confidence band, both, or none should be
  plotted.}

  \item{level}{The desired confidence level.}

  \item{adjust}{A character string indicating the type of
  adjustment (if any) to make to the confidence/prediction
  bands.}

  \item{k}{An integer to be used in computing the critical
  value for the confidence/prediction bands. Only needed
  when \code{adjust = "Bonferroni"} or \code{adjust =
  "Scheffe"}.}

  \item{shade}{A logical value indicating if the band
  should be shaded.}

  \item{extend.range}{A logical value indicating if the
  fitted regression line and bands (if any) should extend
  to the edges of the plot. Default is \code{FALSE}.}

  \item{col.conf}{Shade color for confidence band.}

  \item{col.pred}{Shade color for prediction band.}

  \item{col.fit}{The color to use for the fitted line.}

  \item{border.conf}{The color to use for the confidence
  band border.}

  \item{border.pred}{The color to use for the prediction
  band border.}

  \item{lty.conf}{Line type to use for confidence band
  border.}

  \item{lty.pred}{Line type to use for prediction band
  border.}

  \item{lty.fit}{Line type to use for the fitted regression
  line.}

  \item{lwd.conf}{Line width to use for confidence band
  border.}

  \item{lwd.pred}{Line width to use for prediction band
  border.}

  \item{lwd.fit}{Line width to use for the fitted
  regression line.}

  \item{n}{The number of predictor values at which to
  evaluate the fitted model (larger implies a smoother
  plot).}

  \item{xlab}{A title for the x axis.}

  \item{ylab}{A title for the y axis.}

  \item{xlim}{The x limits (x1, x2) of the plot.}

  \item{ylim}{The y limits (y1, y2) of the plot.}

  \item{hide}{A logical value indicating if the fitted
  model should be plotted on top of the points (FALSE) or
  behind them (TRUE). Default is TRUE.}

  \item{...}{Additional optional arguments passed on to
  \code{plot}.}
}
\description{
  Plots fitted model for an object of class \code{lm} or
  \code{nls}vwith the option of adding a confidence and/or
  prediction band.
}
\note{
  Confidence/prediction bands for nonlinear regression
  (i.e., objects of class \code{nls}) are based on a linear
  approximation as described in Bates & Watts (2007). The
  source code was inpired by the \code{plotfit} function in
  the \code{nlstools} package.
}
\examples{
\donttest{
## A linear regression example
data(cars, package = "datasets")
library(splines)
cars.lm1 <- lm(dist ~ speed, data = cars)
cars.lm2 <- lm(dist ~ speed + I(speed^2), data = cars)
cars.lm3 <- lm(dist ~ speed + I(speed^2) + I(speed^3), data = cars)
cars.lm4 <- lm(dist ~ ns(speed, df = 3), data = cars)
par(mfrow = c(2, 2))
plotFit(cars.lm1, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150),
        main = "linear")
plotFit(cars.lm2, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150),
        main = "quadratic")
plotFit(cars.lm3, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150),
        main = "cubic")
plotFit(cars.lm4, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150),
        main = "cubic spline")

## A nonlinear regression example
par(mfrow = c(1, 1))
library(RColorBrewer) # requires that RColorBrewer be installed
blues <- brewer.pal(9, "Blues")
data(Puromycin, package = "datasets")
Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
                start = c(Vm = 200, K = 0.05))
plotFit(Puro.nls, interval = "both", pch = 19, shade = TRUE,
        col.conf = blues[4], col.pred = blues[2])
}
}
\references{
  Bates, D. M., and Watts, D. G. Nonlinear Regression
  Analysis and its Applications. New York: Wiley, 2007.

  F. Baty and M. L. Delignette-Muller (2012), nlstools:
  Tools for Nonlinear Regression Diagnostics.
}

