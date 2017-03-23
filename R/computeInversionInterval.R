#' @keywords internal
computeInversionInterval <- function(object, ...) {
  UseMethod("computeInversionInterval")
}


#' @keywords internal
computeInversionInterval.lm <- function(object, multi, x0.name, var.pooled, m, 
                                        rat, eta, crit, x0.est, mean.response, 
                                        newdata, lower, upper, extendInt, tol, 
                                        maxiter) {
  
  # Pivotal quantity for linear model (pg. 95)
  rootfun <- function(x) {
    nd <- if (multi) {
      cbind(newdata, makeData(x, x0.name))  # append newdata
    } else {
      makeData(x, x0.name)
    }
    pred <- stats::predict(object, newdata = nd, se.fit = TRUE)
    denom <- if (mean.response) {
      pred$se.fit ^ 2 
    } else {
      var.pooled / m + rat * pred$se.fit ^ 2
    }
    (eta - pred$fit) ^ 2 / denom - crit ^ 2
  }
  
  # Compute lower and upper confidence limits
  lwr <- try(stats::uniroot(rootfun, interval = c(lower, x0.est), 
                            extendInt = extendInt, tol = tol, 
                            maxiter = maxiter)$root, silent = TRUE)
  upr <- try(stats::uniroot(rootfun, interval = c(x0.est, upper), 
                            extendInt = extendInt, tol = tol, 
                            maxiter = maxiter)$root, silent = TRUE)
  
  # Provide (informative) error message if confidence limits not found
  if (inherits(lwr, "try-error")) {
    stop(paste("Lower confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper. ", 
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  if (inherits(upr, "try-error")) {
    stop(paste("Upper confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper. ", 
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  # Return list of results
  list("estimate" = x0.est, "lower" = lwr, "upper" = upr, 
       "interval" = "inversion")
  
}
  

#' @keywords internal
computeInversionInterval.glm <- function(object, multi, x0.name, eta, crit, 
                                         x0.est, mean.response, newdata, lower, 
                                         upper, extendInt, tol, maxiter) {
  
  # Pivotal quantity for generalized linear model
  rootfun <- function(x) {
    nd <- if (multi) {
      cbind(newdata, makeData(x, x0.name))  # append newdata
    } else {
      makeData(x, x0.name)
    }
    pred <- stats::predict(object, newdata = nd, se.fit = TRUE, type = "link")
    (eta - pred$fit) ^ 2 / pred$se.fit ^ 2 - crit ^ 2
  }
  
  # Compute lower and upper confidence limits (i.e., the roots of the 
  # inversion function)
  lwr <- try(stats::uniroot(rootfun, interval = c(lower, x0.est), 
                            extendInt = extendInt, tol = tol, 
                            maxiter = maxiter)$root, silent = TRUE)
  upr <- try(stats::uniroot(rootfun, interval = c(x0.est, upper), 
                            extendInt = extendInt, tol = tol, 
                            maxiter = maxiter)$root, silent = TRUE)
  
  # Provide (informative) error message if confidence limits not found
  if (inherits(lwr, "try-error")) {
    stop(paste("Lower confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper. ", 
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  if (inherits(upr, "try-error")) {
    stop(paste("Upper confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper. ", 
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  # Store results in a list
  res <- list("estimate" = x0.est, "lower" = lwr, "upper" = upr, 
              "interval" = "inversion")
  
}


#' @keywords internal
computeInversionInterval.nls <- function(object, x0.name, var.pooled, m, eta, 
                                         crit, x0.est, mean.response, newdata,
                                         lower, upper, extendInt, tol, 
                                         maxiter) {
  
  # Pivotal quantity for linear model (pg. 95)
  rootfun <- function(x) {
    pred <- predFit(object, newdata = makeData(x, x0.name), se.fit = TRUE) 
    denom <- if (mean.response) {
      pred$se.fit ^ 2 
    } else {
      var.pooled / m + pred$se.fit ^ 2
    }
    (eta - pred$fit) ^ 2 / denom - crit ^ 2
  }
  
  # Compute lower and upper confidence limits
  lwr <- try(stats::uniroot(rootfun, interval = c(lower, x0.est), 
                            extendInt = extendInt, tol = tol, 
                            maxiter = maxiter)$root, silent = TRUE)
  upr <- try(stats::uniroot(rootfun, interval = c(x0.est, upper), 
                            extendInt = extendInt, tol = tol, 
                            maxiter = maxiter)$root, silent = TRUE)
  
  # Provide (informative) error message if confidence limits not found
  if (inherits(lwr, "try-error")) {
    stop(paste("Lower confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper. ", 
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  if (inherits(upr, "try-error")) {
    stop(paste("Upper confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper. ", 
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  # Return list of results
  list("estimate" = x0.est, "lower" = lwr, "upper" = upr, 
       "interval" = "inversion")
  
}


#' @keywords internal
computeInversionInterval.lme <- function(object, x0.name, m, eta, q1, q2, 
                                         x0.est, mean.response, var.y0, lower, 
                                         upper, extendInt, tol, maxiter) {
  
  # Inversion function
  rootfun <- function(x, bound = c("lower", "upper")) {
    pred <- predFit(object, newdata = makeData(x, x0.name), se.fit = TRUE)
    denom <- if (mean.response) {
      pred[, "se.fit"] 
    } else {
      sqrt(var.y0 + pred[, "se.fit"]^2)
    }
    bound <- match.arg(bound)
    if (bound == "upper") {
      (eta - pred[, "fit"]) / denom - q1
    } else {
      (eta - pred[, "fit"]) / denom - q2
    }
  }
  
  # Compute lower and upper confidence limits (i.e., the roots of the 
  # inversion function)
  lwr <- try(stats::uniroot(rootfun, interval = c(lower, x0.est), 
                            bound = "lower", tol = tol, maxiter = maxiter)$root, 
             silent = TRUE)
  upr <- try(stats::uniroot(rootfun, interval = c(x0.est, upper), 
                            bound = "upper", tol = tol, maxiter = maxiter)$root, 
             silent = TRUE)
  
  # Provide (informative) error message if confidence limits not found
  if (inherits(lwr, "try-error")) {
    stop(paste("Lower confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper.", 
               sep = ""), 
         call. = FALSE)
  }
  if (inherits(upr, "try-error")) {
    stop(paste("Upper confidence limit not found in the search interval (", 
               lower, ", ", upper, 
               "). ", "Try tweaking the values of lower and upper.",
               sep = ""), 
         call. = FALSE)
  }
  
  # Store results in a list
  res <- list("estimate" = x0.est, "lower" = lwr, "upper" = upr, 
              "interval" = "inversion")
  
}
