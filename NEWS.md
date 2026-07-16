# investr (development version)

* `nlme` moved from Imports to Suggests; `"lme"`-specific code paths now check for it with `requireNamespace()` (#47).

* Fixed `calibrate.lm()` silently mishandling models with `poly()`-transformed predictors; it now reads the predictor from the model frame and errors informatively instead (#48).

* Fixed the formula for `g` in the inversion-interval vignette, which had been broken by a regression in a prior "typo fix" (#51).

* Fixed `predFit()` with `se.fit = TRUE` failing on `"nls"` models with more than one predictor variable (#53).

* Fixed `invest()` and `predFit()` failing when the fitted model was created inside a function, or with no `data` argument at all; training/prediction data is now reconstructed via `model.frame()` or the environment captured by the model's formula, instead of evaluating `object$call$data` in the caller's frame (#41, #42, #45).

* `invest()` now validates that `newdata` columns have classes (and, for factors, levels) compatible with the fitted model's training data, rather than failing with a confusing downstream error or silently returning a nonsensical result (#36).

* The "Use plotFit for guidance" hint on a failed confidence-limit search no longer appears for multi-predictor `lm`/`glm` models, since `plotFit()` doesn't support more than one predictor variable (#35).

# investr 1.4.2

* Fixed the output of `predFit()` in situations whenever standard errors and confidence/predictions bands are both requested.

* New data sets `bladder` (a repeated measures data set) and `whisky`.

* Fixed typos in documentation throughout the package.

* Updated URLs throughout the package.

* Better y-axis limits when x-axis is on log scale (e.g., using log = "x"). Fixes issue #27.

* Better default y-axis label when using `plotFit()` on a model with transformed response. For example, calling `plotFit(fit)` where `fit <- lm(sqrt(dist) ~ speed, data = cars)` will have a default y-axis label of `"sqrt(dist)"`.

* `plotFit()` has been completely re-written using much less code.

* `predFit()` (and hence `plotFit()`) now works for `"nls"` objects fit using the Golub-Pereyra algorithm (i.e., `algorithm = "plinear"`); however, confidence/prediction bands are still not available.

* New introductory vignette.

# investr 1.4.0

* Multiple predictor variables are allowed for `"lm"` and `"glm"` objects.

* All non-base package functions are now imported.

* The generic function `predFit()` is now exported. This function is used by **investr**
  to obtain predictions, and hence, inverse predictions. For example, `predFit()` 
  can be used to obtain prediction intervals for nonlinear least-squares fits 
  (i.e., models of class `"nls"`).

* Improved tests and test coverage.

* `plotFit()` gained methods for `"rlm"` and `"lqs"` objects from package **MASS**.

# investr 1.3.0

* `invest()` now accepts objects of class `"glm"` (experimental).

* Functions `calibrate()` and `invest()` now return an object of class `"invest"`.

# investr 1.2.1

* Cleaned up documentation.

* Added `AnyNA()` function for those using older versions of R.

# investr 1.2.0

* Cleaned up examples.

* Added `bootstrap` option to `invest()`.

# investr 1.1.2

* Changed tests to satisfy CRAN check.

# investr 1.1.1

* Updated citation file.

* Minor code changes.

* `plotFit()` should now plot models with transformed responses correctly.

* Fixed error causing `invest()` to fail because of a missing data argument.

* Added more tests.

# investr 1.1.0

* `invest()` now accepts objects of class `"lme"` (experimental).

* A few minor bug fixes and code improvements.

* Added more tests.

# investr 1.0.1

* A few minor bug fixes.

* Slightly better documentation.

# investr 1.0

* Initial release.
