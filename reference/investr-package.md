# investr: a package for inverse estimation in R

Inverse estimation, also referred to as the calibration problem, is a
classical and well-known problem in regression. In simple terms, it
involves the use of an observed value of the response (or specified
value of the mean response) to make inference on the corresponding
unknown value of the explanatory variable.

## Details

A detailed [introduction to
investr](https://journal.r-project.org/articles/RJ-2014-009/index.html)
has been published in The R Journal: "investr: An R Package for Inverse
Estimation." You can track development at
<https://github.com/bgreenwell/investr>. To report bugs or issues,
contact the main author directly or submit them to
<https://github.com/bgreenwell/investr/issues>.

As of right now, `investr` supports (univariate) inverse estimation with
objects of class:

- `lm` — linear models (multiple predictor variables allowed)

- `glm` — generalized linear models (multiple predictor variables
  allowed)

- `nls` — nonlinear least-squares models

- `lme` — linear mixed-effects models (fit using the `nlme` package)

- `survfit`/`Surv` — Kaplan-Meier survival curves (fit using the
  `survival` package); see
  [`invest.survfit()`](https://bgreenwell.github.io/investr/reference/invest.survfit.md)

## See also

Useful links:

- <https://github.com/bgreenwell/investr>

- <https://bgreenwell.github.io/investr/>

- <https://bgreenwell.r-universe.dev/investr>

- Report bugs at <https://github.com/bgreenwell/investr/issues>

## Author

**Maintainer**: Brandon M. Greenwell <greenwell.brandon@gmail.com>
([ORCID](https://orcid.org/0000-0002-8120-0084))
