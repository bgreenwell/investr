\name{plotFit}
\alias{plotFit}
\alias{plotFit.lm}
\alias{plotFit.nls}
\title{Plotting Confidence/Prediction Bands}
\usage{
  plotFit(x, ...)

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
  \code{lm}.}

  \item{newdata}{A named list or data frame in which to
  look for variables with which to plot the inference band
  at. If \code{newdata} is missing the inference band is
  plotted over the original data points. At present this
  argument is ignored.}

  \item{interval}{A character string indicating if a
  prediction brand, confidence band, or both should be
  plotted.}

  \item{level}{A numeric scalar between 0 and 1 giving the
  confidence level for the bands to be plotted.}

  \item{shade}{A logical value indicating if the band
  should be shaded.}

  \item{col.conf}{Shade color for confidence band.}

  \item{col.pred}{Shade color for prediction band.}

  \item{col.fit}{The color to use for the fitted line. The
  default, \code{NULL}, means to use
  \code{\link{par}("fg")}}

  \item{border.conf}{The color to use for the confidence
  band border. The default, \code{NULL}, means to use
  \code{\link{par}("fg")}. Use \code{border = NA} to omit
  borders.}

  \item{border.pred}{The color to use for the prediction
  band border. The default, \code{NULL}, means to use
  \code{\link{par}("fg")}. Use \code{border = NA} to omit
  borders.}

  \item{lty.conf}{Line type to use for confidence band
  border.}

  \item{lty.pred}{Line type to use for prediction band
  border.}

  \item{lwd.conf}{Line width to use for confidence band
  border.}

  \item{lwd.pred}{Line width to use for prediction band
  border.}

  \item{n}{The number of predictor values at which to
  compute the inference band(s).}

  \item{xlab}{A title for the x axis: see \link{title}.}

  \item{ylab}{A title for the y axis: see \link{title}.}

  \item{xlim}{The x limits (x1, x2) of the plot. Note that
  x1 > x2 is allowed and leads to a ‘reversed axis’.}

  \item{hide}{A logical value indicating if the inference
  band should be plotted on top of (FALSE) or behind (TRUE)
  the points. Default is TRUE.}

  \item{...}{Additional optional arguments.}
}
\description{
  Plots fitted model (with confidence and/or prediction
  bands) for object of class \code{lm} or \code{nls}.
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
library(RColorBrewer)
blues <- brewer.pal(9, "Blues")
data(Puromycin, package = "datasets")
Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
                start = c(Vm = 200, K = 0.05))
plotFit(Puro.nls, interval = "both", pch = 19, col.conf = blues[4],
        col.pred = blues[2])
}
}

