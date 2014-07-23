## Warning: Using getData2 for now to see if scoping issues go away!

##' Calibration for Linear and Nonlinear Regression Models.
##' 
##' The function \code{invest} computes the inverse estimate and a condfidence 
##' interval for the unknown predictor value that corresponds to an observed 
##' value of the response (or vector thereof) or specified value of the mean 
##' response. See the references listed below for more details. 
##' 
##' @rdname invest
##' @export
##' 
##' @param object An object that inherits from class \code{lm} or \code{nls}.
##' @param y0 The value of the observed response(s) or specified value of the 
##'           mean response.
##' @param interval The type of interval required.
##' @param level A numeric scalar between 0 and 1 giving the confidence level for 
##'              the interval to be calculated. 
##' @param mean.response Logical indicating whether confidence intervals should
##'                      correspond to an individual response (\code{FALSE}) or a 
##'                      mean response
##'        (\code{TRUE}).
##' @param lower The lower endpoint of the interval to be searched.
##' @param upper The upper endpoint of the interval to be searched.
##' @param q1 Optional lower cutoff to be used in forming confidence intervals. 
##'           Only used when object inherits from class \code{lme}. Defaults to
##'           \code{qnorm((1+level)/2)}.
##' @param q2 Optional upper cutoff to be used in forming confidence intervals. 
##'           Only used when object inherits from class \code{lme}. Defaults to
##'           \code{qnorm((1-level)/2)}.
##' @param tol The desired accuracy passed on to \code{uniroot}. Recommend a 
##'            minimum of 1e-10.
##' @param maxiter The maximum number of iterations passed on to \code{uniroot}. 
##' @param adjust A logical value indicating if an adjustment should be made to
##'               the critical value used in calculating the confidence interval.
##'               This is useful for when the calibration curve is to be used 
##'               multiple, say k, times.
##' @param k The number times the calibration curve is to be used for computing a 
##'          confidence interval. Only needed when \code{adjust = "Bonferroni"}.
##' @param ... Additional optional arguments. At present, no optional arguments 
##'            are used.
##' @return An object of class \code{calibrate} containing the following 
##'         components:
##' \describe{
##'   \item{\code{estimate}}{The estimate of x0.}
##'   \item{\code{lwr}}{The lower confidence limit for x0.}
##'   \item{\code{upr}}{The upper confidence limit for x0.}
##'   \item{\code{se}}{An estimate of the standard error (Wald interval only).}
##'   \item{\code{interval}}{The method used for calculating \code{lower} and 
##'                   \code{upper} (only used by \code{print} method).}
##' }
##' @references
##' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
##' Applications. Belmont, Calif: Duxbury Press, 1994. 
##'
##' Huet, S., Bouvier, A., Poursat, M-A., and Jolivet, E. Statistical Tools for 
##' Nonlinear Regression: A Practical Guide with S-PLUS and R Examples. New York: 
##' Springer, 2004. 
##' 
##' Seber, G. A. F., and Wild, C. J.. Nonlinear regression. New York: Wiley, 
##' 1989.
##' @examples
##' data(Puromycin, package = "datasets")
##' Puromycin2 <- Puromycin[Puromycin$state == "treated", ]
##' Puro2.nls <- nls(rate ~ (theta1 * conc) / (theta2 + conc), 
##'                  data = Puromycin2, start = list(theta1 = 200, theta2 = 1))
##' plotFit(Puro2.nls, interval = "both")
##' invest(Puro2.nls, y0 = 100, interval = "inversion")
##' invest(Puro2.nls, y0 = 100, interval = "inversion", mean.response = TRUE)
invest <- function(object, ...) {
  UseMethod("invest")
} 

##' @rdname invest
##' @export
##' @method invest lm
invest.lm <- function(object, y0, interval = c("inversion", "Wald", "none"), 
                      level = 0.95, mean.response = FALSE, lower, upper, 
                      tol = .Machine$double.eps^0.25, maxiter = 1000,  
                      adjust = c("none", "Bonferroni"), k,  ...) {
  

  vars <- getVarInfo(object)  # extract data and variable information
  xname <- vars$x.names  # predictor variable name
  yname <- vars$y.names  # response variable name
  if (missing(lower)) lower <- min(vars$x)  # search grid lower limit default
  if (missing(upper)) upper <- max(vars$x)  # search grid upper limit default
  eta <- mean(y0)  # mean unknown
  m <- length(y0)  # number of unknowns 
  n <- length(resid(object))  # in case of missing values
  p <- length(coef(object))   # number of regression coefficients
  df1 <- n - p  # stage 1 degrees of freedom
  df2 <- m - 1  # stage 2 degrees of freedom
  var1 <- Sigma(object)^2  # stage 1 variance estimate
  var2 <- if (m == 1) 0 else var(y0)  # stage 2 variance estimate
  var.pooled <- (df1*var1 + df2*var2)/(df1 + df2)  # pooled estimate of variance
  rat <- var.pooled/var1  # right variance?
  
  ## Try to catch errors
  if (length(xname) != 1) stop("Only one independent variable allowed.") 
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  
  ## Calculate point estimate by inverting fitted model
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(object, x)) - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  ## Provide (informative) error message if point estimate is not found
  if (inherits(x0.est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  ## Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.est)
  
  ## Critical value for confidence interval computations
  adjust <- match.arg(adjust)
  crit <- if (adjust == "Bonferroni" && m == 1) qt((level + 2*k - 1) / (2*k), n+m-p-1)
          else qt((level+1) / 2, n+m-p-1)
  
  ## inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    ## Inversion function
    inversionFun <- function(x) {
      pred <- predict(object, makeData(object, x), se.fit = TRUE)
      denom <- if (mean.response) pred$se.fit^2 else var.pooled/m + 
        rat*pred$se.fit^2
      (eta - pred$fit)^2/denom - crit^2
    }
    
    ## Compute lower and upper confidence limits (i.e., the roots of the 
    ## inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0.est), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0.est, upper), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    
    ## Provide (informative) error message if confidence limits not found
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
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
  } 
  
  ## Wald interval -------------------------------------------------------------
  if (interval == "Wald") { 
    
    ## Function of parameters whose gradient is required
    object.copy <- object # FIXME: Is a copy really needed?
    dmFun <- function(params) {
      if (mean.response) {
        object.copy$coefficients <- params
        z <- eta
      } else {
        object.copy$coefficients <- params[-length(params)]
        z <- params[length(params)]
      }
      uniroot(function(x) { 
        predict(object.copy, makeData(object.copy, x)) - z
      }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root
    }
    
    ## Variance-covariane matrix
    if (mean.response) {
      params <- coef(object)
      covmat <- vcov(object)
    } else {
      params <- c(coef(object), eta)
      covmat <- diag(p + 1)
      covmat[p + 1, p + 1] <- var.pooled/m
      covmat[1:p, 1:p] <- vcov(object)
    }
    
    ## Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = x0.est - crit * se, 
                "upper" = x0.est + crit * se, 
                "se" = se,
                "interval" = interval)
  }
  
  ## Assign class label and return results
  class(res) <- "calibrate"
  res
  
}

##' @rdname invest
##' @export
##' @method invest nls
invest.nls <- function(object, y0, interval = c("inversion", "Wald", "none"),  
                       level = 0.95, mean.response = FALSE, lower, upper, 
                       tol = .Machine$double.eps^0.25, maxiter = 1000, 
                       adjust = c("none", "Bonferroni"), k, ...) {
  
  ## No support for the Golub-Pereyra algorithm for partially linear 
  ## least-squares models
  if (object$call$algorithm == "plinear") {
    stop(paste("The Golub-Pereyra algorithm for partially linear least-squares 
               models is currently not supported."))
  }
  
  vars <- getVarInfo(object)  # extract data and variable information
  xname <- vars$x.names  # predictor variable name
  yname <- vars$y.names  # response variable name
  if (missing(lower)) lower <- min(vars$x)  # search grid lower limit default
  if (missing(upper)) upper <- max(vars$x)  # search grid upper limit default
  eta <- mean(y0)  # mean response
  m <- length(y0)  # number of unknowns 
  n <- length(resid(object))  # sample size
  p <- length(coef(object))  # number of parameters
  var.pooled <- Sigma(object)^2  # residual variance
  
  ## Try to catch errors
  if (vars$x.dim != 1) stop("Only one independent variable allowed.")
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  
  ## Critical value for confidence interval computations
  adjust <- match.arg(adjust)
  crit <- if (adjust == "Bonferroni" && m == 1) qt((level + 2*k - 1) / (2*k), n+m-p-1)
          else qt((level + 1) / 2, n+m-p-1)

  ## Calculate point estimate by inverting fitted model
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(object, x)) - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  ## Provide (informative) error message if point estimate is not found
  if (inherits(x0.est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  ## Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.est)
  
  ## Inversion interval --------------------------------------------------------
  if (interval == "inversion") {
    
    ## Inversion function
    inversionFun <- function(x) {
      pred <- predict2(object, makeData(object, x)) # FIXME:, se.fit = TRUE)
      denom <- if (mean.response) pred$se.fit^2 else (var.pooled/m + pred$se.fit^2)
      (eta - pred$fit)^2/denom - crit^2
    }
    
    ## Compute lower and upper confidence limits (i.e., the roots of the 
    ## inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0.est), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0.est, upper), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    
    ## Provide (informative) error message if confidence limits not found
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
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
  } 
  
  # Wald interval --------------------------------------------------------------
  if (interval == "Wald") { 
    
    ## Function of parameters whose gradient is required
    object.copy <- object # FIXME: Is a copy really needed?
    dmFun <- function(params) {
      if (mean.response) {
        object.copy$m$setPars(params)
        z <- eta
      } else {
        object.copy$m$setPars(params[-length(params)])
        z <- params[length(params)]
      }
      uniroot(function(x) { 
        predict(object.copy, makeData(object.copy, x)) - z
      }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root
    }
    
    ## Variance-covariance matrix
    if (mean.response) {
      params <- coef(object)
      covmat <- vcov(object)
    } else {
      params <- c(coef(object), eta)
      covmat <- diag(p + 1)
      covmat[p + 1, p + 1] <- var.pooled/m
      covmat[1:p, 1:p] <- vcov(object)
    }
    
    ## Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = x0.est - crit * se, 
                "upper" = x0.est + crit * se, 
                "se" = se,
                "interval" = interval)
  }
  
  ## Assign class label and return results
  class(res) <- "calibrate"
  res
  
}

##' @rdname invest
##' @export
##' @method invest lme
invest.lme <- function(object, y0, interval = c("inversion", "Wald", "none"),  
                       level = 0.95, mean.response = FALSE, lower, upper, q1, 
                       q2, tol = .Machine$double.eps^0.25, maxiter = 1000, ...) 
{
  
  vars <- getVarInfo(object)  # extract data and variable information
  xname <- vars$x.names  # predictor variable name
  yname <- vars$y.names  # response variable name
  if (missing(lower)) lower <- min(vars$x)  # search grid lower limit default
  if (missing(upper)) upper <- max(vars$x)  # search grid upper limit default
  eta <- mean(y0)
  m <- length(y0)
  if (m != 1) stop('Only a single unknown allowed for objects of class "lme".')
  N <- length(resid(object)) 
  p <- length(fixef(object))
#   res.var <- Sigma(object)^2  # residual variance
  
  ## Try to catch errors
  if (vars$x.dim != 1) stop("Only one independent variable allowed.")
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  
  ## Critical value. Oman (1998. pg. 445) suggests a t(1-alpha/2, N-1) dist.
  if (missing(q1)) q1 <- qnorm((1-level) / 2)
  if (missing(q2)) q2 <- qnorm((1+level) / 2)
  
  ## Calculate point estimate by inverting fitted model
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(object, x), level = 0) - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  ## Provide (informative) error message if point estimate is not found
  if (inherits(x0.est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper.", 
               sep = ""), 
         call. = FALSE)
  }
  
  ## Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.est)
  
  ## Estimate variance of new response
  if (!mean.response) var.y0 <- varY(object, makeData(object, x0.est))
  
  ## Inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    ## Inversion function
    inversionFun <- function(x, bound = c("lower", "upper")) {
      bound <- match.arg(bound)
      pred <- predict2(object, makeData(object, x), se.fit = TRUE)
      denom <- if (mean.response) pred$se.fit else sqrt(var.y0 + pred$se.fit^2)
      if (bound == "upper") {
        (eta - pred$fit)/denom - q1
      } else {
        (eta - pred$fit)/denom - q2
      }
    }
    
    ## Compute lower and upper confidence limits (i.e., the roots of the 
    ## inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0.est), 
                       bound = "lower", tol = tol, maxiter = maxiter)$root, 
               silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0.est, upper), 
                       bound = "upper", tol = tol, maxiter = maxiter)$root, 
               silent = TRUE)
    
    ## Provide (informative) error message if confidence limits not found
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
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
    
  }
  
  ## Wald interval -------------------------------------------------------------
  if (interval == "Wald") { 
    
    ## Function of parameters whose gradient is required
    dmFun <- function(params) {
      fun <- function(x) {
        X <- model.matrix(eval(object$call$fixed)[-2], 
                          data = makeData(object, x))
        if (mean.response) {
          X %*% params - eta
        } else {
          X %*% params[-length(params)] - params[length(params)]
        }
      }
      uniroot(fun, lower = lower, upper = upper, tol = tol, 
              maxiter = maxiter)$root
    }
    
    ## Variance-covariance matrix
    if (mean.response) {
      params <- fixef(object)
      covmat <- vcov(object)
    } else {
      params <- c(fixef(object), eta)
      covmat <- diag(p + 1)
      covmat[p + 1, p + 1] <- var.y0
      covmat[1:p, 1:p] <- vcov(object)
    }
    
    ## Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = x0.est - q2 * se, 
                "upper" = x0.est + q2 * se,
                "se" = se,
                "interval" = interval)
    
  } 
  
  ## Assign class label and return results
  class(res) <- "calibrate"
  res
  
}
