#' @keywords internal
computeWaldInterval <- function(object, ...) {
  UseMethod("computeWaldInterval")
}


#' @keywords internal
computeWaldInterval.lm <- function(object, multi, x0.name, var.pooled, m, p,
                                   eta, crit, x0.est, mean.response, newdata,
                                   lower, upper, extendInt, tol, maxiter) {
  
  # Make a copy of the fitted model
  object.copy <- object # FIXME: Is a copy really needed?
  
  # Function of parameters whose gradient is required
  parFun <- function(params) {
    
    # Replace appropriate parameters
    if (mean.response) {  # regulation
      object.copy$coefficients <- params
      z <- eta
    } else {  # calibration (last elements of params corresponds to eta)
      object.copy$coefficients <- params[-length(params)]
      z <- params[length(params)]
    }
    
    # Find inverese estimate for current value of params
    stats::uniroot(function(x) { 
      nd <- if (multi) {
        cbind(newdata, makeData(x, x0.name))  # append newdata
      } else {
        makeData(x, x0.name)
      }
      stats::predict(object.copy, newdata = nd) - z
    }, interval = c(lower, upper), extendInt = extendInt, 
    tol = tol, maxiter = maxiter)$root
    
  }
  
  # Set up variance-covariane matrix
  if (mean.response) {  # regulation
    params <- stats::coef(object)
    covmat <- stats::vcov(object)
  } else {  # calibration (includes Y0)
    params <- c(stats::coef(object), eta)
    covmat <- diag(p + 1)
    covmat[p + 1, p + 1] <- var.pooled / m
    covmat[1:p, 1:p] <- stats::vcov(object)
  }
  
  # Calculate gradient and standard error
  gv <- attr(stats::numericDeriv(quote(parFun(params)), "params"), "gradient")
  se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))  # standard formula
  
  # Store results in a list
  res <- list("estimate" = x0.est,           # original point estimate
              "lower" = x0.est - crit * se,  # standard Wald-based interval
              "upper" = x0.est + crit * se,  # standard Wald-based interval
              "se" = se,                     # large sample standard error
              "interval" = "Wald")           # type of interval
  
  
}


#' @keywords internal
computeWaldInterval.glm <- function(object, multi, x0.name, eta, crit, x0.est, 
                                    newdata, lower, upper, extendInt, tol, 
                                    maxiter) {
  
  # Make a copy of the fitted model
  object.copy <- object # FIXME: Is a copy really needed?
  
  # Function of parameters whose gradient is required
  parFun <- function(params) {
    object.copy$coefficients <- params
    z <- eta
    stats::uniroot(function(x) { 
      nd <- if (multi) {
        cbind(newdata, makeData(x, x0.name))  # append newdata
      } else {
        makeData(x, x0.name)
      }
      stats::predict(object.copy, newdata = nd, type = "link") - z
    }, interval = c(lower, upper), extendInt = extendInt, tol = tol, 
    maxiter = maxiter)$root
  }
  
  # Variance-covariane matrix
  params <- stats::coef(object)
  covmat <- stats::vcov(object)
  
  # Calculate gradient, and return standard error
  gv <- attr(stats::numericDeriv(quote(parFun(params)), "params"), "gradient")
  se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
  
  # Store results in a list
  res <- list("estimate" = x0.est,           # original point estimate
              "lower" = x0.est - crit * se,  # standard Wald-based interval
              "upper" = x0.est + crit * se,  # standard Wald-based interval
              "se" = se,                     # large sample standard error
              "interval" = "Wald")           # type of interval

}


#' @keywords internal
computeWaldInterval.nls <- function(object, x0.name, var.pooled, m, p, eta, 
                                    crit, x0.est, mean.response, newdata,
                                    lower, upper, extendInt, tol, maxiter) {
  
  # Make a copy of the fitted model
  object.copy <- object # FIXME: Is a copy really needed?
  
  # Function of parameters whose gradient is required
  parFun <- function(params) {
    if (mean.response) {
      object.copy$m$setPars(params)
      z <- eta
    } else {
      object.copy$m$setPars(params[-length(params)])
      z <- params[length(params)]
    }
    stats::uniroot(function(x) { 
      stats::predict(object.copy, newdata = makeData(x, x0.name)) - z
    }, interval = c(lower, upper), extendInt = extendInt, tol = tol, 
    maxiter = maxiter)$root
  }
  
  # Variance-covariance matrix
  if (mean.response) {
    params <- stats::coef(object)
    covmat <- stats::vcov(object)
  } else {
    params <- c(stats::coef(object), eta)
    covmat <- diag(p + 1)
    covmat[p + 1, p + 1] <- var.pooled/m
    covmat[1:p, 1:p] <- stats::vcov(object)
  }
  
  # Calculate gradient and standard error
  gv <- attr(stats::numericDeriv(quote(parFun(params)), "params"), "gradient")
  se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))  # standard formula
  
  # Store results in a list
  res <- list("estimate" = x0.est,           # original point estimate
              "lower" = x0.est - crit * se,  # standard Wald-based interval
              "upper" = x0.est + crit * se,  # standard Wald-based interval
              "se" = se,                     # large sample standard error
              "interval" = "Wald")           # type of interval

}


#' @keywords internal
computeWaldInterval.lme <- function(object, x0.name, p, eta, q1, q2, x0.est, 
                                    mean.response, var.y0, lower, upper, 
                                    extendInt, tol, maxiter) {
  
  # Make a copy of the fitted model
  object.copy <- object # FIXME: Is a copy really needed?
  
  # Function of parameters whose gradient is required
  parFun <- function(params) {
    fun <- function(x) {
      X <- stats::model.matrix(eval(object$call$fixed)[-2], 
                               data = makeData(x, x0.name))
      if (mean.response) {
        X %*% params - eta
      } else {
        X %*% params[-length(params)] - params[length(params)]
      }
    }
    stats::uniroot(fun, lower = lower, upper = upper, tol = tol, 
                   maxiter = maxiter)$root
  }
  
  # Variance-covariance matrix
  if (mean.response) {
    params <- nlme::fixef(object)
    covmat <- stats::vcov(object)
  } else {
    params <- c(nlme::fixef(object), eta)
    covmat <- diag(p + 1)
    covmat[p + 1, p + 1] <- var.y0
    covmat[1:p, 1:p] <- stats::vcov(object)
  }
  
  # Calculate gradient, and return standard error
  gv <- attr(stats::numericDeriv(quote(parFun(params)), "params"), "gradient")
  se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
  
  # Store results in a list
  res <- list("estimate" = x0.est, 
              "lower" = x0.est - q2 * se, 
              "upper" = x0.est + q2 * se, 
              "se" = se,
              "interval" = "Wald")
  
} 
