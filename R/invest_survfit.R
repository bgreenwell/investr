#' Inverse estimation for Kaplan-Meier survival curves
#'
#' Inverts a Kaplan-Meier estimate of the survival function: given a survival
#' probability `y0`, these methods estimate the survival time \eqn{t} at which
#' \eqn{S(t) = y_0}, along with a confidence interval obtained by inverting
#' the pointwise confidence band for the survival curve (the same inversion
#' idea used by [invest()] for regression models). For example, the default
#' `y0 = 0.5` gives the median survival time.
#'
#' @param object An object that inherits from class
#' `"survfit"` (a fitted Kaplan-Meier curve from
#' `survival::survfit()`) or class `"Surv"` (a censored response
#' object from `survival::Surv()`, for which a Kaplan-Meier curve is fit
#' internally).
#'
#' @param y0 A numeric scalar in (0, 1) giving the survival probability at
#' which to invert the curve. Default is `0.5` (i.e., the median
#' survival time).
#'
#' @param level A numeric scalar between 0 and 1 giving the confidence level
#' for the interval. For `"survfit"` objects, this must match the
#' `conf.int` level the curve was fit with (the default); refit with
#' `survival::survfit(..., conf.int = level)` for a different level. For
#' `"Surv"` objects, the internal fit uses `level` directly.
#'
#' @param ... Additional optional arguments. At present, no optional
#' arguments are used.
#'
#' @return An object of class `"invest"` containing the components
#' `estimate`, `lower`, `upper`, and `interval` (always
#' `"inversion"`). Components are `NA` whenever the survival curve
#' (or a confidence limit) never drops below `y0`; in particular,
#' `upper` is frequently `NA` in small samples, meaning the upper
#' confidence limit is indeterminate (infinite).
#'
#' @note These methods only support a single survival curve; for a
#' `"survfit"` object with multiple strata, invert each stratum's curve
#' separately (e.g., via `survfit[i]` subscripting).
#'
#' @seealso The `quantile.survfit` method in the \pkg{survival} package,
#' which these methods use for the inversion and which expresses the same
#' calculation in terms of quantiles of the event-time distribution (i.e.,
#' `probs = 1 - y0`).
#'
#' @rdname invest.survfit
#'
#' @export
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'
#'   # Median survival time (with 95% confidence limits) for the acute
#'   # myelogenous leukemia data
#'   km <- survival::survfit(survival::Surv(time, status) ~ 1,
#'                           data = survival::aml)
#'   invest(km)  # y0 = 0.5 is the median
#'
#'   # Same, but fitting the curve internally from a censored response
#'   invest(survival::Surv(survival::aml$time, survival::aml$status))
#'
#'   # Time at which 75% survival is reached
#'   invest(km, y0 = 0.75)
#'
#' }
invest.survfit <- function(object, y0 = 0.5, level = object$conf.int, ...) {

  # "survfit" objects require the survival package to be available
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package \"survival\" is required for objects of class ",
         "\"survfit\". Please install it.", call. = FALSE)
  }
  if (length(y0) != 1L || !is.numeric(y0) || y0 <= 0 || y0 >= 1) {
    stop("'y0' must be a single survival probability in (0, 1).",
         call. = FALSE)
  }
  if (!is.null(object$strata)) {
    stop("Only a single survival curve is supported; invert each stratum ",
         "separately (e.g., invest(object[1], ...)).", call. = FALSE)
  }

  # The confidence limits are baked into the survfit object at fit time, so
  # a different level requires refitting the curve
  if (!isTRUE(all.equal(level, object$conf.int))) {
    stop("'level' (", level, ") does not match the confidence level this ",
         "curve was fit with (", object$conf.int, "). Refit with ",
         "survival::survfit(..., conf.int = ", level, ").", call. = FALSE)
  }

  # Invert the curve and its pointwise confidence band; in quantile terms,
  # S(t) = y0 corresponds to the (1 - y0) quantile of the event-time
  # distribution
  q <- stats::quantile(object, probs = 1 - y0, conf.int = TRUE)
  res <- list("estimate" = unname(q$quantile),
              "lower"    = unname(q$lower),
              "upper"    = unname(q$upper),
              "interval" = "inversion")
  if (is.na(res$estimate)) {
    warning("The survival curve never drops below y0 = ", y0, "; the ",
            "estimated survival time is indeterminate.", call. = FALSE)
  }
  class(res) <- "invest"
  res

}


#' @rdname invest.survfit
#'
#' @export
invest.Surv <- function(object, y0 = 0.5, level = 0.95, ...) {

  # "Surv" objects require the survival package to be available
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package \"survival\" is required for objects of class \"Surv\". ",
         "Please install it.", call. = FALSE)
  }

  # Fit a single Kaplan-Meier curve at the requested confidence level, then
  # invert it
  km <- survival::survfit(object ~ 1, conf.int = level)
  invest.survfit(km, y0 = y0, level = level, ...)

}
