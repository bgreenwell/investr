#' Calibration for the simple linear regression model.
#' 
#' The function \code{calibrate} computes the maximum likelihood estimate and a
#' condfidence interval for the unknown predictor value that corresponds to an 
#' observed value of the response (or vector thereof) or specified value of the 
#' mean response. See the reference listed below for more details.
#'  
#' @param object An object that inherits from class \code{lm}, a matrix, a list, 
#' or a data frame.
#' @param formula A formula of the form \code{y ~ x}.
#' @param data an optional data frame, list or environment (or object coercible 
#' by \code{as.data.frame} to a data frame) containing the variables in the 
#' model. If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which \code{lm}
#' is called. 
#' @param subset An optional vector specifying a subset of observations to be 
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data 
#' contain \code{NA}s. 
#' @param y0 The value of the observed response(s) or specified value of the
#'           mean response.
#' @param interval The method to use for forming a confidence interval.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for 
#'              the interval to be calculated. 
#' @param mean.response Logicial indicating whether confidence intervals should 
#' correspond to an observed response(s) (\code{FALSE}) or a specified value of 
#' the mean response (\code{TRUE}). Default is \code{FALSE}.
#' @param adjust A logical value indicating if an adjustment should be made to
#'               the critical value used in calculating the confidence interval.
#'               This useful for when the calibration curve is to be used 
#'               multiple, say k, times.
#' @param k The number times the calibration curve is to be used for computing a 
#'          confidence interval. Only needed when \code{adjust = TRUE}.
#' @param ... Additional optional arguments. At present, no optional arguments 
#'            are used.
#' @return An object of class \code{calibrate} containing the following 
#'         components:
#' \itemize{
#'   \item{estimate}{The maximum likelihood estimate of x0.}
#'   \item{lwr}{The lower confidence bound on x0.}
#'   \item{upr}{The upper confidence bound on x0.}
#'   \item{se}{An estimate of the standard error (Wald interval only).}
#'   \item{interval}{The method used for calculating \code{lower} and 
#'                   \code{upper}.}
#' }
#' @references 
#' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
#' Applications. Belmont, Calif: Duxbury Press, 1994.
#' @rdname calibrate
#' @aliases print.calibrate
#' @export
#' @section Warning:
#'   You must not call this function unless ...
#' @note The function \code{invest} is more general, but based on numerical
#' techniques to find the solution. When the underlying model is that of the 
#' simple linear regression model with normal errors, closed-form expressions
#' exist which are utilized by the function \code{calibrate}.
#' @examples
#' \donttest{
#' ## Inverting a prediction interval for an individual response
#' arsenic.lm <- lm(measured ~ actual, data = arsenic)
#' plotFit(arsenic.lm, interval = "prediction", shade = TRUE, 
#'         col.pred = "lightblue")
#' calibrate(arsenic.lm, y0 = 3, interval = "inversion")
#' 
#' ## Inverting a confidence interval for the mean response
#' crystal.lm <- lm(weight ~ time, data = crystal)
#' plotFit(crystal.lm, interval = "confidence", shade = TRUE,
#'         col.conf = "lightblue")
#' calibrate(crystal.lm, y0 = 8, interval = "inversion", mean.response = TRUE)
#'
#' ## Wald interval and approximate standard error based on the delta method
#' calibrate(crystal.lm, y0 = 8, interval = "Wald", mean.response = TRUE)
#' 
#' ## Alterntively, we can use the car package to compute the standard error (this 
#' ## is trickier though when mean.respone = FALSE, hence, it is better to use the
#' ## calibrate function).
#' library(car)
#' deltaMethod(crystal.lm, g = "(8 - b0) / b1", parameterNames = c("b0", "b1"))
#' }
# calibrate <- function(z, ...) {
#   ## TODO:
#   ##   (1) Add option for maximum modulus intervals; adjust = "maximum".
#   if (is.null(class(z))) class(z) <- data.class(z)
#   UseMethod("calibrate")
# } 

calibrate <- function(object, ...) {
  UseMethod("calibrate")
}

#' @rdname calibrate
#' @export
#' @method calibrate default
calibrate.default <- function(object, y0, interval = c("inversion", "Wald"), 
                              level = 0.95, mean.response = FALSE, 
                              adjust = c("none", "Bonferroni", "Scheffe"), k, 
                              ...) {

  ## Extract needed components from fitted model
  if (inherits(object, "matrix")) {
    x <- object[, 1]
    y <- object[, 2]
  } else if (inherits(object, "data.frame")) {
    object <- data.matrix(object)
    x <- object[, 1]
    y <- object[, 2]
  } else if (inherits(object, "list")) {
    x <- object[[1]]
    y <- object[[2]]
    if (length(x) != length(y)) {
      stop(paste("Components of '", deparse(substitute(object)), 
                 "' not of same length.", sep = ""))
    }
  } else {
    stop(paste("The object '", deparse(substitute(object)), 
               "' is not a valid matrix, list, or data frame.", sep = ""))
  }
  eta <- mean(y0)             # mean of new observations
  m <- length(y0)             # number of new observations
  z <- lm(y ~ x, data = data.frame(x, y)) # fit simple linear regression model
  e <- resid(z)               # residuals
  n <- length(e)              # sample size
  v1 <- n - 2                 # stage I degrees of freedom
  v2 <- m - 1                 # stage II degrees of freedom
  v <- v1 + v2                # total degrees of freedom
  u1 <- summary(z)$sigma^2    # stage I variance estimate
  u2 <- if (m == 1) 0 else var(y0) # stage II variance estimate
  u <- (v1*u1 + v2*u2)/v      # pooled estimate of variance
  sigma <- sqrt(u)            # sqrt of pooled variance estimate
  b <- as.numeric(coef(z))    # regression coefficients
  x0.mle <- (eta - b[1L])/b[2L] # MLE of x0
  ssx <- sum((x - mean(x))^2)   # sum-of-squares for x, Sxx
  alpha <- 1 - level            # significance level
  interval <- match.arg(interval)
  
  ## Try to catch errors
  if (mean.response && m > 1) {
    stop("Only one mean response value allowed.")
  }
  
  # Adjustment for simultaneous intervals
  adjust <- match.arg(adjust)
  w <- if (adjust == "Bonferroni" && m == 1) {
         qt(1 - alpha/(2*k), n+m-3)
       } else if (adjust == "Scheffe" && m == 1) {
         sqrt(k * qf(1 - alpha, k, n+m-3))
       } else {
         qt(1 - alpha/2, n+m-3)
       }

  ## Calculate inversion interval
  if (interval == "inversion") { 

    c1 <- b[2L]^2 - (sigma^2 * w^2)/ssx
    c2 <- if (mean.response) {
      c1/n + (eta - mean(y))^2/ssx
    } else {
      c1*(1/m + 1/n) + (eta - mean(y))^2/ssx
    }
    c3 <- b[2L] * (eta - mean(y))
    c4 <- w * sigma
    
    ## FIXME: catch errors and throw an appropriate warning
    if (c1 < 0 && c2 <= 0) {
      
      warning("The calibration line is not well determined.", call. = FALSE)
      lwr <- -Inf
      upr <- Inf
      
    } else {
      
      lwr <- mean(x) + (c3 - c4*sqrt(c2))/c1
      upr <- mean(x) + (c3 + c4*sqrt(c2))/c1
      if (c1 < 0 && c2 > 0) {
        
        stop(paste("The calibration line is not well determined.\nReturning two semi-infinite intervals:\n(", -Inf, ",", 
                   round(upr, 4), ") and (", round(lwr, 4), ",", Inf, ")"), 
             call. = FALSE)
        
      }
      
    }
    res <- list("estimate" = x0.mle,
                "lower" = lwr,
                "upper" = upr,
                "interval" = interval)
    
  ## Calculate Wald interval  
  } else if (interval == "Wald") { 
    
    ## Compute standard error for Wald interval
    se <- if (mean.response) {
      abs((sigma/b[2]))*sqrt((1/n + (x0.mle - mean(x))^2/ssx))
    } else {
      abs((sigma/b[2]))*sqrt((1/m + 1/n + (x0.mle - mean(x))^2/ssx))
    }

    ## Store results in a list
    res <- list("estimate" = x0.mle,
                "lower" = x0.mle - w*se,
                "upper" = x0.mle + w*se,
                "se" = se, 
                "interval" = interval)
    
  } 
  
  ## Assign class label and return results
  class(res) <- "calibrate"
  return(res)
  
} 

#' @rdname calibrate
#' @export
#' @method calibrate formula
calibrate.formula <- function(formula, data = NULL, ..., subset, 
                              na.action = na.fail) {
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, sys.parent()))) {
    m$data <- as.data.frame(data)
  }
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  y <- model.extract(m, "response")
  mm <- model.matrix(Terms, m)
  if (ncol(mm) > 1) stop("This function only works for the simple linear regression model (i.e., y ~ x).")
  x <- as.numeric(mm)
  calibrate(cbind(x, y), ...)
} 

#' @rdname calibrate
#' @export
#' @method calibrate lm
calibrate.lm <- function(object, ...) {
  d <- eval(object$call$data) # may be null
  calibrate(formula(object), data = d, ...)
} 

#' @keywords internal
print.calibrate <- function(x, digits = 4, ...) {
  if (x$interval == "inversion") print(round(unlist(x[1:3]), digits))
  if (x$interval == "Wald") print(round(unlist(x[1:4]), digits))
  invisible(x)
} 