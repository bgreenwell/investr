# Plotting fitted models

Generic function for plotting predictions from various types of fitted
models. `plotFit` currently supports objects of class
[`stats::lm()`](https://rdrr.io/r/stats/lm.html),
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), and
[`stats::nls()`](https://rdrr.io/r/stats/nls.html). A default method
also exists which may be used for plotting the fitted mean response from
other model fits (e.g.,
[`MASS::lqs()`](https://rdrr.io/pkg/MASS/man/lqs.html) and
[`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) from the **MASS**
package).

## Usage

``` r
plotFit(object, ...)

# Default S3 method
plotFit(
  object,
  type = c("response", "link"),
  interval = c("none", "both", "confidence", "prediction"),
  level = 0.95,
  data,
  adjust = c("none", "Bonferroni", "Scheffe"),
  k,
  ...,
  shade = FALSE,
  extend.range = FALSE,
  hide = TRUE,
  col.conf = if (shade) grDevices::grey(0.7) else "black",
  col.pred = if (shade) grDevices::grey(0.9) else "black",
  border.conf = col.conf,
  border.pred = col.pred,
  col.fit = "black",
  lty.conf = if (shade) 1 else 2,
  lty.pred = if (shade) 1 else 3,
  lty.fit = 1,
  lwd.conf = 1,
  lwd.pred = 1,
  lwd.fit = 1,
  n = 500,
  xlab,
  ylab,
  xlim,
  ylim
)
```

## Arguments

- object:

  A fitted model object. Typically, an object that inherits from class
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html),
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html), or
  [`stats::nls()`](https://rdrr.io/r/stats/nls.html) (but others may
  work too).

- ...:

  Additional optional arguments passed on to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- type:

  The type of prediction required. The default is on the scale of the
  response variable; the alternative `"link"` is on the scale of the
  linear predictor. This option is only used when plotting
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) objects.

- interval:

  A character string indicating if a prediction band, confidence band,
  both, or none should be plotted.

- level:

  The desired confidence level.

- data:

  An optional data frame containing the variables in the model.

- adjust:

  A character string indicating the type of adjustment (if any) to make
  to the confidence/prediction bands.

- k:

  An integer to be used in computing the critical value for the
  confidence/prediction bands. Only needed when `adjust = "Bonferroni"`,
  or when `adjust = "Scheffe"` and `interval = "prediction"`.

- shade:

  A logical value indicating if the band should be shaded.

- extend.range:

  A logical value indicating if the fitted regression line and bands (if
  any) should extend to the edges of the plot. Default is `FALSE`.

- hide:

  A logical value indicating if the fitted model should be plotted on
  top of the points (`FALSE`) or behind them (`TRUE`). Default is
  `TRUE`.

- col.conf:

  Shade color for confidence band.

- col.pred:

  Shade color for prediction band.

- border.conf:

  The color to use for the confidence band border.

- border.pred:

  The color to use for the prediction band border.

- col.fit:

  The color to use for the fitted line.

- lty.conf:

  Line type to use for confidence band border.

- lty.pred:

  Line type to use for prediction band border.

- lty.fit:

  Line type to use for the fitted regression line.

- lwd.conf:

  Line width to use for confidence band border.

- lwd.pred:

  Line width to use for prediction band border.

- lwd.fit:

  Line width to use for the fitted regression line.

- n:

  The number of predictor values at which to evaluate the fitted model
  (larger gives a smoother plot).

- xlab:

  A title for the x axis.

- ylab:

  A title for the y axis.

- xlim:

  The x limits (x1, x2) of the plot.

- ylim:

  The y limits (y1, y2) of the plot.

## Value

No return value (called for side effects).

## Note

By default, the plotted intervals are unadjusted (i.e., pointwise)
intervals. For simultaneous intervals, use `adjust = "Bonferroni"` or
`adjust = "Scheffe"`. For the Bonferroni adjustment, you must specify a
value for `k`, the number of intervals for which the coverage is to hold
simultaneously. For the Scheffe adjustment, specifying a value for `k`
is only required when `interval = "prediction"`; if
`interval = "confidence"`, `k` is set equal to \\p\\, the number of
regression parameters. For example, if `object` is a simple linear
regression model, then calling `plotFit` with `interval = "confidence"`
and `adjust = "Scheffe"` will plot the [Working-Hotelling
band](https://en.wikipedia.org/wiki/Working-Hotelling_procedure).

Confidence/prediction bands for nonlinear regression (i.e., objects of
class [`stats::nls()`](https://rdrr.io/r/stats/nls.html)) are based on
the linear approximation described in Bates & Watts (2007).

## References

Bates, D. M., and Watts, D. G. (2007) *Nonlinear Regression Analysis and
its Applications*. Wiley.

Florent Baty, Christian Ritz, Sandrine Charles, Martin Brutsche,
Jean-Pierre Flandrois, Marie-Laure Delignette-Muller (2015). A Toolbox
for Nonlinear Regression in R: The Package nlstools. *Journal of
Statistical Software*, **66**(5), 1-21.

## See also

The `plotfit` function in the nlstools package.

## Examples

``` r
# A nonlinear least squares example (see ?datasets::Puromycin and 
# ?investr::predFit)
data(Puromycin, package = "datasets")
Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
                start = c(Vm = 200, K = 0.05))
plotFit(Puro.nls, interval = "both", pch = 19, shade = TRUE, 
        col.conf = "skyblue4", col.pred = "lightskyblue2")  
#> Error in eval(stats::getCall(object)$data): object 'Puromycin2' not found
```
