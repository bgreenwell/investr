#' Calibration for Linear and Nonlinear Regression Models.
#' 
#' The function \code{invest} computes the inverse estimate and a condfidence 
#' interval for the unknown predictor value that corresponds to an observed 
#' value of the response (or vector thereof) or specified value of the mean 
#' response. See the references listed below for more details. 
#' 
#' @rdname invest
#' @export
#' 
#' @param object An object that inherits from class \code{lm} or \code{nls}.
#' @param y0 The value of the observed response(s) or specified value of the 
#'           mean response.
#' @param interval The type of interval required.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for 
#'              the interval to be calculated. 
#' @param mean.response Logical indicating whether confidence intervals should
#'                      correspond to an individual response (\code{FALSE}) or a 
#'                      mean response
#'        (\code{TRUE}).
#' @param lower The lower endpoint of the interval to be searched.
#' @param upper The upper endpoint of the interval to be searched.
#' @param tol The desired accuracy passed on to \code{uniroot}.
#' @param maxiter The maximum number of iterations passed on to \code{uniroot}. 
#' (\code{TRUE}).
#' @param adjust A logical value indicating if an adjustment should be made to
#'               the critical value used in calculating the confidence interval.
#'               This is useful for when the calibration curve is to be used 
#'               multiple, say k, times.
#' @param k The number times the calibration curve is to be used for computing a 
#'          confidence interval. Only needed when \code{adjust = "Bonferroni"}.
#' @param ... Additional optional arguments. At present, no optional arguments 
#'            are used.
#' @return An object of class \code{calibrate} containing the following 
#'         components:
#' \itemize{
#'   \item{estimate}{The estimate of x0.}
#'   \item{lwr}{The lower confidence bound on x0.}
#'   \item{upr}{The upper confidence bound on x0.}
#'   \item{se}{An estimate of the standard error (Wald interval only).}
#'   \item{interval}{The method used for calculating \code{lower} and 
#'                   \code{upper}.}
#' }
#' @references
#' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
#' Applications. Belmont, Calif: Duxbury Press, 1994. 
#'
#' Huet, S., Bouvier, A., Poursat, M-A., and Jolivet, E. Statistical Tools for 
#' Nonlinear Regression: A Practical Guide with S-PLUS and R Examples. New York: 
#' Springer, 2004. 
#' 
#' Seber, G. A. F., and Wild, C. J.. Nonlinear regression. New York: Wiley, 
#' 1989.
#' @examples
#' data(Puromycin, package = "datasets")
#' Puromycin2 <- Puromycin[Puromycin$state == "treated", ]
#' Puro2.nls <- nls(rate ~ (theta1 * conc) / (theta2 + conc), 
#'                  data = Puromycin2, start = list(theta1 = 200, theta2 = 1))
#' plotFit(Puro2.nls, interval = "both")
#' invest(Puro2.nls, y0 = 100, interval = "inversion")
#' invest(Puro2.nls, y0 = 100, interval = "inversion", mean.response = TRUE)
invest <- function(object, ...) {
  UseMethod("invest")
} 

#' @rdname invest
#' @export
#' @method invest lm
invest.lm <- function(object, y0, interval = c("inversion", "Wald"), 
                      level = 0.95, mean.response = FALSE, lower, upper, 
                      tol = .Machine$double.eps^0.25, maxiter = 1000,  
                      adjust = c("none", "Bonferroni"), k,  ...) {
  
  ## Extract data, variables, etc.
  d <- eval(object$call$data, sys.frame())
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (missing(lower)) {
    lower <- min(d[, xvar])
  }
  if (missing(upper)) {
    upper <- max(d[, xvar])
  }
  interval <- match.arg(interval)
  alpha <- 1 - level
  eta <- mean(y0)
  m <- length(y0)
  n <- length(resid(object)) # in case of missing values
  p <- length(coef(object))
  
  ## Calculate correct variance
  ## FIXME: is this the correct variance to use for all univariate linear 
  ## models? For example, is this the correct variance for a quadratic fit?
  v1 <- n - p                      # stage I degrees of freedom
  v2 <- m - 1                      # stage II degrees of freedom
  u1 <- summary(object)$sigma^2    # stage I variance estimate
  u2 <- if (m == 1) 0 else var(y0) # stage II variance estimate
  u <- (v1*u1 + v2*u2)/(v1 + v2)   # pooled estimate of variance
  rat <- u/u1                      # temporary fix 
  
  ## Try to catch errors
  if (length(xvar) != 1) {
    stop("only one independent variable allowed")
  }
  if (mean.response && m > 1) {
    stop("only one mean response value allowed")
  }
  
  # Adjustment for simultaneous intervals
  adjust <- match.arg(adjust)
  w <- if (adjust == "Bonferroni" && m == 1) {
         qt(1 - alpha/(2 * k), n+m-p-1)
       } else {
         qt(1 - alpha/2, n+m-p-1)
       }
  
  ## Compute point estimate by "inverting" the fitted model at y = eta
  invFun.est <- function(x) {
    z <- list(x)
    names(z) <- xvar
    predict(object, newdata = z) - eta
  }
  x0.est <- uniroot(invFun.est, interval = c(lower, upper), tol = tol, 
                    maxiter = maxiter)$root

  ## Compute interval estimate
  if (interval == "inversion") { ## inversion interval
    
    ## "Invert" confidence/prediction band at y = eta
    if (mean.response) { ## "Invert" confidence band (regulation)
      invFun <- function(x) {
        z <- list(x); names(z) <- xvar
        pred <- predict(object, newdata = z, se.fit = TRUE)
        (eta - pred$fit)^2/(pred$se.fit^2) - w^2
      }
    } else { ## "Invert" prediction band (calibration)
      invFun <- function(x) {
        z <- list(x); names(z) <- xvar
        pred <- predict(object, newdata = z, se.fit = TRUE)
        (eta - pred$fit)^2/(u/m + rat*pred$se.fit^2) - w^2
      }
    }
    
    ## Compute lower and upper endpoints of confidence interval
    lwr <- uniroot(invFun, interval = c(lower, x0.est), tol = tol, 
                   maxiter = maxiter)$root
    upr <- uniroot(invFun, interval = c(x0.est, upper), tol = tol, 
                   maxiter = maxiter)$root
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
  } else { ## Wald interval
    
    ## Delta method based on a modification of deltaMethod function from the car
    ## package
    se <- if (mean.response) { ## regulation
      
            ## Function of parameters whose gradient is required
            dmFun <- function(params) {
              object$coefficients <- params
              invFun <- function(x) {
                z <- list(x)
                names(z) <- xvar
                predict(object, z) - eta
              }
              uniroot(invFun, interval = c(lower, upper), tol = tol, 
                      maxiter = maxiter)$root
            }
          
            ## Assign parameter names, calculate gradient, and return standard  
            ## error
            params <- coef(object)
            gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
            #gv <- t(grad(dmFun, params))
            as.numeric(sqrt(gv %*% vcov(object) %*% t(gv)))

            } else { ## calibration
            
            ## Function of parameters whose gradient is required
            dmFun <- function(params) {
              object$coefficients <- params ## FIXME: params[1:length(params)-1]
              invFun <- function(x) {
                z <- list(x)
                names(z) <- xvar 
                predict(object, z) - params[length(params)]
              }
              uniroot(invFun, interval = c(lower, upper), tol = tol, 
                      maxiter = maxiter)$root
            }
            
            ## Assign parameter names, calculate gradient, and return standard 
            ## error
            params <- c(coef(object), eta)
            covmat <- diag(p + 1)
            covmat[p + 1, p + 1] <- u/m
            covmat[1:p, 1:p] <- vcov(object)
            gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
            #gv <- t(grad(dmFun, params))
            as.numeric(sqrt(gv %*% covmat %*% t(gv)))
            
          }
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = x0.est - se * w, 
                "upper" = x0.est + se * w, 
                "se" = se,
                "interval" = interval)
  }
  
  ## Assign class label and return results
  class(res) <- "calibrate"
  res
  
}

#' @rdname invest
#' @export
#' @method invest nls
invest.nls <- function(object, y0, interval = c("inversion", "Wald"),  
                       level = 0.95, mean.response = FALSE, lower, upper, 
                       tol = .Machine$double.eps^0.25, maxiter = 1000, 
                       adjust = c("none", "Bonferroni"), k, ...) 
{
  
  ## Extract data, variables, etc.
  d <- eval(object$call$data, sys.frame())
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (missing(lower)) {
    lower <- min(d[, xvar])
  }
  if (missing(upper)) {
    upper <- max(d[, xvar])
  }
  interval <- match.arg(interval)
  alpha <- 1 - level
  eta <- mean(y0)
  m <- length(y0)
  n <- length(resid(object)) 
  p <- length(coef(object))
  
  ## Calculate variance
  u <- summary(object)$sigma^2 
  
  ## Try to catch errors
  if (length(xvar) != 1) {
    stop("only one independent variable allowed")
  }
  if(mean.response && m > 1) {
    stop("only one value of the mean response is allowed")
  }
  
  # Adjustment for simultaneous intervals
  adjust <- match.arg(adjust)
  w <- if (adjust == "Bonferroni" && m == 1) {
         qt(1 - alpha/(2 * k), n+m-p-1)
       } else {
         qt(1 - alpha/2, n+m-p-1)
       }
  
  ## Compute point estimate by "inverting" the fitted model at y = eta
  invFun.est <- function(x) {
    z <- list(x)
    names(z) <- xvar
    predict(object, newdata = z) - eta
  }
  x0.est <- uniroot(invFun.est, interval = c(lower, upper), tol = tol, 
                    maxiter = maxiter)$root
  
  ## Compute interval estimate
  if (interval == "inversion") { # inversion interval
    
    ## "Invert" confidence/prediction band at y = eta
    if (mean.response) { # "Invert" confidence band (regulation)
      invFun <- function(x) {
        z <- list(x); names(z) <- xvar
        pred <- predict2(object, newdata = z)
        (eta - pred$fit)^2/(pred$se.fit^2) - w^2
      }
    } else { # calibration
      invFun <- function(x) { # "Invert" prediction band (calibration)
        z <- list(x); names(z) <- xvar
        pred <- predict2(object, newdata = z)
        (eta - pred$fit)^2/(u/m + pred$se.fit^2) - w^2
      }
    }
    
    ## Compute lower and upper endpoints of confidence interval
    lwr <- uniroot(invFun, interval = c(lower, x0.est), tol = tol, 
                   maxiter = maxiter)$root
    upr <- uniroot(invFun, interval = c(x0.est, upper), tol = tol, 
                   maxiter = maxiter)$root
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
  } else { # Wald interval
    
    ## Calculate standard error based on the delta method 
    object.copy <- object
    se <- if (mean.response) { # regulation
      
            ## Function of parameters whose gradient is required
            dmFun <- function(params) {
              object.copy$m$setPars(params)
              invFun <- function(x) {
                z <- list(x)
                names(z) <- xvar
                predict2(object.copy, z)$fit - eta
              }
              uniroot(invFun, interval = c(lower, upper), tol = tol, 
                      maxiter = maxiter)$root
            }
            
            ## Assign parameter names, calculate gradient, and return standard 
            ## error
            params <- coef(object)
            gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
            #gv <- t(grad(dmFun, params))
            as.numeric(sqrt(gv %*% vcov(object) %*% t(gv)))
            
          } else { # calibration

            ## Function of parameters whose gradient is required
            dmFun <- function(params) {
              object.copy$m$setPars(params)
              invFun <- function(x) {
                z <- list(x)
                names(z) <- xvar
                predict2(object.copy, z)$fit - params[length(params)]
              }
              uniroot(invFun, interval = c(lower, upper), tol = tol, 
                      maxiter = maxiter)$root
            }
            
            ## Assign parameter names, calculate gradient, and return standard 
            ## error
            params <- c(coef(object), eta)
            covmat <- diag(p + 1)
            covmat[p + 1, p + 1] <- u/m
            covmat[1:p, 1:p] <- vcov(object)
            gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
            #gv <- t(grad(dmFun, params))
            as.numeric(sqrt(gv %*% covmat %*% t(gv)))
            
          }
    
    ## Store results in a list
    res <- list("estimate" = x0.est, 
                "lower" = x0.est - w * se, 
                "upper" = x0.est + w * se, 
                "se" = se,
                "interval" = interval)
  }

  ## Assign class label and return results
  class(res) <- "calibrate"
  res
  
}

#' @rdname invest
#' @export
#' @method invest lme
invest.lme <- function(object, y0, interval = c("inversion", "Wald", "pboot"),  
                       level = 0.95, mean.response = FALSE, lower, upper, 
                       tol = .Machine$double.eps^0.25, maxiter = 1000, ...) 
{
  
  ## Extract data, variables, etc.
  d <- eval(object$call$data, sys.frame())
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (missing(lower)) {
    lower <- min(d[, xvar])
  }
  if (missing(upper)) {
    upper <- max(d[, xvar])
  }
  interval <- match.arg(interval)
  alpha <- 1 - level
  eta <- mean(y0)
  m <- length(y0)
  n <- length(resid(object)) 
  p <- length(fixef(object))
  
  ## Calculate variance
  u <- summary(object)$sigma^2 # residual variance
  
  ## Try to catch errors
  if (length(xvar) != 1) {
    stop("only one independent variable allowed")
  }
  if(mean.response && m > 1) {
    stop("only one value of the mean response is allowed")
  }
  
  # Critical value
  w <- qnorm(1 - alpha/2)
  
  ## Compute point estimate by "inverting" the fitted model at y = eta
  invFun.est <- function(x) {
    z <- list(x)
    names(z) <- xvar
    predict(object, newdata = z, level = 0) - eta
  }
  x0.est <- uniroot(invFun.est, interval = c(lower, upper), tol = tol, 
                    maxiter = maxiter)$root
  
  x0.est
  
#   ## Compute interval estimate
#   if (interval == "inversion") { # inversion interval
#     
#     ## "Invert" confidence/prediction band at y = eta
#     if (mean.response) { # "Invert" confidence band (regulation)
#       invFun <- function(x) {
#         z <- list(x); names(z) <- xvar
#         pred <- predict2(object, newdata = z)
#         (eta - pred$fit)^2/(pred$se.fit^2) - w^2
#       }
#     } else { # calibration
#       invFun <- function(x) { # "Invert" prediction band (calibration)
#         z <- list(x); names(z) <- xvar
#         pred <- predict2(object, newdata = z)
#         (eta - pred$fit)^2/(u/m + pred$se.fit^2) - w^2
#       }
#     }
#     
#     ## Compute lower and upper endpoints of confidence interval
#     lwr <- uniroot(invFun, interval = c(lower, x0.est), tol = tol, 
#                    maxiter = maxiter)$root
#     upr <- uniroot(invFun, interval = c(x0.est, upper), tol = tol, 
#                    maxiter = maxiter)$root
#     
#     ## Store results in a list
#     res <- list("estimate" = x0.est, 
#                 "lower" = lwr, 
#                 "upper" = upr, 
#                 "interval" = interval)
#   } else { # Wald interval
#     
#     ## Calculate standard error based on the delta method 
#     object.copy <- object
#     se <- if (mean.response) { # regulation
#       
#       ## Function of parameters whose gradient is required
#       dmFun <- function(params) {
#         object.copy$m$setPars(params)
#         invFun <- function(x) {
#           z <- list(x)
#           names(z) <- xvar
#           predict2(object.copy, z)$fit - eta
#         }
#         uniroot(invFun, interval = c(lower, upper), tol = tol, 
#                 maxiter = maxiter)$root
#       }
#       
#       ## Assign parameter names, calculate gradient, and return standard 
#       ## error
#       params <- coef(object)
#       gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
#       #gv <- t(grad(dmFun, params))
#       as.numeric(sqrt(gv %*% vcov(object) %*% t(gv)))
#       
#     } else { # calibration
#       
#       ## Function of parameters whose gradient is required
#       dmFun <- function(params) {
#         object.copy$m$setPars(params)
#         invFun <- function(x) {
#           z <- list(x)
#           names(z) <- xvar
#           predict2(object.copy, z)$fit - params[length(params)]
#         }
#         uniroot(invFun, interval = c(lower, upper), tol = tol, 
#                 maxiter = maxiter)$root
#       }
#       
#       ## Assign parameter names, calculate gradient, and return standard 
#       ## error
#       params <- c(coef(object), eta)
#       covmat <- diag(p + 1)
#       covmat[p + 1, p + 1] <- u/m
#       covmat[1:p, 1:p] <- vcov(object)
#       gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
#       #gv <- t(grad(dmFun, params))
#       as.numeric(sqrt(gv %*% covmat %*% t(gv)))
#       
#     }
#     
#     ## Store results in a list
#     res <- list("estimate" = x0.est, 
#                 "lower" = x0.est - w * se, 
#                 "upper" = x0.est + w * se, 
#                 "se" = se,
#                 "interval" = interval)
#   }
#   
#   ## Assign class label and return results
#   class(res) <- "calibrate"
#   res
  
}
