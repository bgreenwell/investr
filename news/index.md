# Changelog

## investr 1.5.0

- New
  [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  methods for `"survfit"` and `"Surv"` objects (from the survival
  package) estimate the survival time at which a Kaplan-Meier curve
  reaches a given survival probability (e.g., median survival time),
  with confidence limits obtained by inverting the curve’s pointwise
  confidence band
  ([\#30](https://github.com/bgreenwell/investr/issues/30)).

- [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  results with `interval = "percentile"` now also inherit from class
  `"boot"`, so they work directly with functions from the boot package
  such as [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html)
  ([\#32](https://github.com/bgreenwell/investr/issues/32)).

- `nlme` moved from Imports to Suggests; `"lme"`-specific code paths now
  check for it with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)
  ([\#47](https://github.com/bgreenwell/investr/issues/47)).

- Fixed
  [`calibrate.lm()`](https://bgreenwell.github.io/investr/reference/calibrate.md)
  silently mishandling models with
  [`poly()`](https://rdrr.io/r/stats/poly.html)-transformed predictors;
  it now reads the predictor from the model frame and errors
  informatively instead
  ([\#48](https://github.com/bgreenwell/investr/issues/48)).

- Fixed the formula for `g` in the inversion-interval vignette, which
  had been broken by a regression in a prior “typo fix”
  ([\#51](https://github.com/bgreenwell/investr/issues/51)).

- Fixed
  [`predFit()`](https://bgreenwell.github.io/investr/reference/predFit.md)
  with `se.fit = TRUE` failing on `"nls"` models with more than one
  predictor variable
  ([\#53](https://github.com/bgreenwell/investr/issues/53)).

- Fixed
  [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  and
  [`predFit()`](https://bgreenwell.github.io/investr/reference/predFit.md)
  failing when the fitted model was created inside a function, or with
  no `data` argument at all; training/prediction data is now
  reconstructed via
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) or the
  environment captured by the model’s formula, instead of evaluating
  `object$call$data` in the caller’s frame
  ([\#41](https://github.com/bgreenwell/investr/issues/41),
  [\#42](https://github.com/bgreenwell/investr/issues/42),
  [\#45](https://github.com/bgreenwell/investr/issues/45)).

- [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  now validates that `newdata` columns have classes (and, for factors,
  levels) compatible with the fitted model’s training data, rather than
  failing with a confusing downstream error or silently returning a
  nonsensical result
  ([\#36](https://github.com/bgreenwell/investr/issues/36)).

- The “Use plotFit for guidance” hint on a failed confidence-limit
  search no longer appears for multi-predictor `lm`/`glm` models, since
  [`plotFit()`](https://bgreenwell.github.io/investr/reference/plotFit.md)
  doesn’t support more than one predictor variable
  ([\#35](https://github.com/bgreenwell/investr/issues/35)).

## investr 1.4.2

CRAN release: 2022-03-31

- Fixed the output of
  [`predFit()`](https://bgreenwell.github.io/investr/reference/predFit.md)
  in situations whenever standard errors and confidence/predictions
  bands are both requested.

- New data sets `bladder` (a repeated measures data set) and `whisky`.

- Fixed typos in documentation throughout the package.

- Updated URLs throughout the package.

- Better y-axis limits when x-axis is on log scale (e.g., using log =
  “x”). Fixes issue
  [\#27](https://github.com/bgreenwell/investr/issues/27).

- Better default y-axis label when using
  [`plotFit()`](https://bgreenwell.github.io/investr/reference/plotFit.md)
  on a model with transformed response. For example, calling
  `plotFit(fit)` where `fit <- lm(sqrt(dist) ~ speed, data = cars)` will
  have a default y-axis label of `"sqrt(dist)"`.

- [`plotFit()`](https://bgreenwell.github.io/investr/reference/plotFit.md)
  has been completely re-written using much less code.

- [`predFit()`](https://bgreenwell.github.io/investr/reference/predFit.md)
  (and hence
  [`plotFit()`](https://bgreenwell.github.io/investr/reference/plotFit.md))
  now works for `"nls"` objects fit using the Golub-Pereyra algorithm
  (i.e., `algorithm = "plinear"`); however, confidence/prediction bands
  are still not available.

- New introductory vignette.

## investr 1.4.0

CRAN release: 2016-04-09

- Multiple predictor variables are allowed for `"lm"` and `"glm"`
  objects.

- All non-base package functions are now imported.

- The generic function
  [`predFit()`](https://bgreenwell.github.io/investr/reference/predFit.md)
  is now exported. This function is used by **investr** to obtain
  predictions, and hence, inverse predictions. For example,
  [`predFit()`](https://bgreenwell.github.io/investr/reference/predFit.md)
  can be used to obtain prediction intervals for nonlinear least-squares
  fits (i.e., models of class `"nls"`).

- Improved tests and test coverage.

- [`plotFit()`](https://bgreenwell.github.io/investr/reference/plotFit.md)
  gained methods for `"rlm"` and `"lqs"` objects from package **MASS**.

## investr 1.3.0

CRAN release: 2015-03-25

- [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  now accepts objects of class `"glm"` (experimental).

- Functions
  [`calibrate()`](https://bgreenwell.github.io/investr/reference/calibrate.md)
  and
  [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  now return an object of class `"invest"`.

## investr 1.2.1

CRAN release: 2015-01-06

- Cleaned up documentation.

- Added `AnyNA()` function for those using older versions of R.

## investr 1.2.0

CRAN release: 2014-12-23

- Cleaned up examples.

- Added `bootstrap` option to
  [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md).

## investr 1.1.2

CRAN release: 2014-09-15

- Changed tests to satisfy CRAN check.

## investr 1.1.1

CRAN release: 2014-09-09

- Updated citation file.

- Minor code changes.

- [`plotFit()`](https://bgreenwell.github.io/investr/reference/plotFit.md)
  should now plot models with transformed responses correctly.

- Fixed error causing
  [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  to fail because of a missing data argument.

- Added more tests.

## investr 1.1.0

CRAN release: 2014-07-14

- [`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
  now accepts objects of class `"lme"` (experimental).

- A few minor bug fixes and code improvements.

- Added more tests.

## investr 1.0.1

CRAN release: 2014-04-19

- A few minor bug fixes.

- Slightly better documentation.

## investr 1.0

CRAN release: 2013-08-09

- Initial release.
