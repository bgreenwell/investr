# Inverse estimation for Kaplan-Meier survival curves

Inverts a Kaplan-Meier estimate of the survival function: given a
survival probability `y0`, these methods estimate the survival time
\\t\\ at which \\S(t) = y_0\\, along with a confidence interval obtained
by inverting the pointwise confidence band for the survival curve (the
same inversion idea used by
[`invest()`](https://bgreenwell.github.io/investr/reference/invest.md)
for regression models). For example, the default `y0 = 0.5` gives the
median survival time.

## Usage

``` r
# S3 method for class 'survfit'
invest(object, y0 = 0.5, level = object$conf.int, ...)

# S3 method for class 'Surv'
invest(object, y0 = 0.5, level = 0.95, ...)
```

## Arguments

- object:

  An object that inherits from class `"survfit"` (a fitted Kaplan-Meier
  curve from
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html))
  or class `"Surv"` (a censored response object from
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html), for
  which a Kaplan-Meier curve is fit internally).

- y0:

  A numeric scalar in (0, 1) giving the survival probability at which to
  invert the curve. Default is `0.5` (i.e., the median survival time).

- level:

  A numeric scalar between 0 and 1 giving the confidence level for the
  interval. For `"survfit"` objects, this must match the `conf.int`
  level the curve was fit with (the default); refit with
  `survival::survfit(..., conf.int = level)` for a different level. For
  `"Surv"` objects, the internal fit uses `level` directly.

- ...:

  Additional optional arguments. At present, no optional arguments are
  used.

## Value

An object of class `"invest"` containing the components `estimate`,
`lower`, `upper`, and `interval` (always `"inversion"`). Components are
`NA` whenever the survival curve (or a confidence limit) never drops
below `y0`; in particular, `upper` is frequently `NA` in small samples,
meaning the upper confidence limit is indeterminate (infinite).

## Note

These methods only support a single survival curve; for a `"survfit"`
object with multiple strata, invert each stratum's curve separately
(e.g., via `survfit[i]` subscripting).

## See also

The `quantile.survfit` method in the survival package, which these
methods use for the inversion and which expresses the same calculation
in terms of quantiles of the event-time distribution (i.e.,
`probs = 1 - y0`).

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {

  # Median survival time (with 95% confidence limits) for the acute
  # myelogenous leukemia data
  km <- survival::survfit(survival::Surv(time, status) ~ 1,
                          data = survival::aml)
  invest(km)  # y0 = 0.5 is the median

  # Same, but fitting the curve internally from a censored response
  invest(survival::Surv(survival::aml$time, survival::aml$status))

  # Time at which 75% survival is reached
  invest(km, y0 = 0.75)

}
#> estimate    lower    upper 
#>       12        8       30 
```
