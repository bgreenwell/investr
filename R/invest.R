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
##' @param object An object that inherits from class \code{"lm"}, \code{"nls"}, 
##'               or \code{"lme"}.
##' @param y0 The value of the observed response(s) or specified value of the 
##'           mean response.
##' @param interval The type of interval required.
##' @param level A numeric scalar between 0 and 1 giving the confidence level for 
##'              the interval to be calculated. 
##' @param mean.response Logical indicating whether confidence intervals should
##'                      correspond to an individual response (\code{FALSE}) or a 
##'                      mean response
##'        (\code{TRUE}).
##' @param data An optional data frame. This is required if \code{object$data} 
##'             is \code{NULL}.
##' @param boot Logical indicating whether to carry out a bootstrap simulation.
##' @param nsim Positive integer specifying the number of bootstrap simulations; 
##'             the bootstrap B (or R).
##' @param type Character string specifying the type of bootstrap, 
##'             \code{"parametric"} or \code{"nonparametric"}.
##' @param seed Optional argument to \code{set.seed}.
##' @param progress Logical indicating whether to display a text-based progress
##'                 bar during the bootstrap simulation.
##' @param lower The lower endpoint of the interval to be searched.
##' @param upper The upper endpoint of the interval to be searched.
##' @param q1 Optional lower cutoff to be used in forming confidence intervals. 
##'           Only used when \code{object} inherits from class \code{"lme"}. Defaults to
##'           \code{qnorm((1+level)/2)}.
##' @param q2 Optional upper cutoff to be used in forming confidence intervals. 
##'           Only used when \code{object} inherits from class \code{"lme"}. Defaults to
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
##' @return If \code{boot = FALSE}, then an object of class \code{"calibrate"} 
##'         containing the following components:
##'   \itemize{
##'     \item \code{estimate} The estimate of x0.
##'     \item \code{lwr} The lower confidence limit for x0.
##'     \item \code{upr} The upper confidence limit for x0.
##'     \item \code{se} An estimate of the standard error (Wald interval only).
##'     \item \code{interval} The method used for calculating \code{lower} and 
##'           \code{upper} (only used by \code{print} method).
##'   }
##'         Otherwise, an object of class \code{"bootCal"} containing the 
##'         following components:
##'   \itemize{
##'     \item \code{original} The estimate of x0.
##'     \item \code{bootreps} The lower confidence limit for x0.
##'     \item \code{nsim} The number of bootstrap replicates.
##'     \item \code{level} The desired confidence level.
##'   }
##' @references
##' Greenwell, B. M., and Schubert Kabban, C. M. (2014). investr: An R Package 
##' for Inverse Estimation. \emph{The R Journal}, \bold{6}(1), 90--100. 
##' URL http://journal.r-project.org/archive/2014-1/greenwell-kabban.pdf.
##'
##' Graybill, F. A., and Iyer, H. K. (1994).
##' \emph{Regression analysis: Concepts and Applications}. Duxbury Press. 
##'
##' Huet, S., Bouvier, A., Poursat, M-A., and Jolivet, E.  (2004)
##' \emph{Statistical Tools for Nonlinear Regression: A Practical Guide with S-PLUS and R Examples}. Springer. 
##' 
##' Seber, G. A. F., and Wild, C. J. (1989)
##' \emph{Nonlinear regression}. Wiley.
##' 
##' Oman, Samuel D. (1998).
##' Calibration with Random Slopes.
##' \emph{Biometrics} \bold{85}(2): 439--449.
##' doi:10.1093/biomet/85.2.439.
##' @examples
##' ##
##' ## Nasturtium example (nonlinear regression with replication)
##' ##
##' 
##' ## Log-logistic model
##' mod <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
##'            start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
##'            data = nasturtium)
##' plotFit(mod, lwd.fit = 2)
##'            
##' ## Compute approximate 95% calibration intervals
##' invest(mod, y0 = c(309, 296, 419), interval = "inversion")
##' invest(mod, y0 = c(309, 296, 419), interval = "Wald")  
##' 
##' ## Bootstrap calibration intervals. In general, nsim should be as large as 
##' ## reasonably possible (say, nsim = 9999).
##' boo <- invest(mod, y0 = c(309, 296, 419), boot = TRUE, nsim = 999, 
##'               seed = 101)
##' boo  # print bootstrap summary
##' plot(boo)  # plot results
invest <- function(object, ...) {
  UseMethod("invest")
} 

##' @rdname invest
##' @export
##' @method invest lm
invest.lm <- function(object, y0, interval = c("inversion", "Wald", "none"), 
                      level = 0.95, mean.response = FALSE, data, boot = FALSE,
                      type = c("parametric", "nonparametric"), nsim = 1, 
                      seed = NULL, progress = FALSE, lower, upper, 
                      tol = .Machine$double.eps^0.25, maxiter = 1000, 
                      adjust = c("none", "Bonferroni"), k,  ...) {
  
  ## Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (missing(lower)) lower <- min(.data[, xname])  # lower limit default
  if (missing(upper)) upper <- max(.data[, xname])  # upper limit default
  
  ## Set up for inverse estimation
  m <- length(y0)  # number of unknowns 
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  eta <- mean(y0)  # mean unknown
  n <- length(resid(object))  # in case of missing values
  p <- length(coef(object))   # number of regression coefficients
  df1 <- n - p  # stage 1 degrees of freedom
  df2 <- m - 1  # stage 2 degrees of freedom
  var1 <- Sigma(object)^2  # stage 1 variance estimate
  var2 <- if (m == 1) 0 else var(y0)  # stage 2 variance estimate
  var.pooled <- (df1*var1 + df2*var2) / (df1 + df2)  # pooled estimate
  rat <- var.pooled / var1  # right variance?
  
  ## Calculate point estimate by inverting fitted model
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(x, xname)) - eta
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
  
  ## Parametric bootstrap ------------------------------------------------------
  if (boot) {
    
    ## Sanity check
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    
    ## Set up progress bar (if requested)
    if (progress) { 
      pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    }
    
    ## Initialize random number generator
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    if (is.null(seed)) 
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    
    ## Simulate new response vectors
    ftd <- fitted(object)  # fitted values
    res <- residuals(object) # redisuals
    type <- match.arg(type)
    if (type == "parametric") {  
      ss <- simulate(object, nsim = nsim)
    } else {
      ss <- replicate(nsim, ftd + sample(res, replace = TRUE), simplify = FALSE)
    }
    
    ## Function to calculate inverse estimate
    x0Fun <- function(i) {
      
      ## Update model using simulated response data
      boot.data <- eval(object$call$data)  # copy data
      boot.data[, yname] <- ss[[i]]  # simulated response vector
      boot.object <- tryCatch(update(object, data = boot.data),
                              error = function(e) NULL)
      
      ## If updating the model fails, then return value is NA
      if (is.null(boot.object)) {
        ret <- NA
      } else {
        
        ## Simulate new response (different from simulated response vector)
        if (mean.response) {  # regulation
          y0.star <- y0  # hold constant in bootstrap replications
        } else {  # calibration
          if (type == "parametric") {
            y0.star <- y0 + rnorm(length(y0), sd = Sigma(object))
          } else {
            y0.star <- y0 + sample(res, size = length(y0), replace = TRUE)
          }
        }
        
        ## Calculate point estimate
        ret <- tryCatch(uniroot(function(x) {
          predict(boot.object, newdata = makeData(x, xname)) - mean(y0.star)
        }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
        error = function(e) NA)
      }
      
      ## Update progress bar
      if (progress) { 
        setTxtProgressBar(pb, i) 
      }
      
      ## Return estimate
      ret
      
    }
    
    ## Calculate bootstrap replicates
    x0.star <- sapply(seq_len(nsim), x0Fun)
    
    ## Check for errors and return the runs that did not fail
    if (AnyNA(x0.star)) {
      num_fail <- sum(is.na(x0.star))
      warning("some bootstrap runs failed (", num_fail, "/", nsim, 
              ")")
      x0.star <- na.omit(x0.star)  # remove runs that failed
      attributes(x0.star) <- NULL  # remove attributes
    } else {
      num_fail <- NULL
    }
    
    ## Create and return a bootCal object (essentially a list that can later be 
    ## converted to an object of class boot)
    res <- list(original  = x0.est,  # original estimate
                bootreps  = x0.star,  # bootstrap replicates
                nsim      = nsim,  # number of simulations
                level     = level)  # desired confidence level
    class(res) = "bootCal"
    attr(res, "bootFail") <- num_fail
    return(res)
     
  }
  
  ## Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.est)
  
  ## Critical value for confidence interval computations
  adjust <- match.arg(adjust)
  crit <- if (adjust == "Bonferroni" && m == 1) {
            qt((level + 2*k - 1) / (2*k), n+m-p-1)  # Bonferroni adjustment
          } else {
            qt((level+1) / 2, n+m-p-1)  # no adjustment
          }
  
  ## inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    ## Inversion function
    inversionFun <- function(x) {
      pred <- predict(object, newdata = makeData(x, xname), se.fit = TRUE)
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
#         predict(object.copy, makeData(object.copy, x)) - z
        predict(object.copy, newdata = makeData(x, xname)) - z
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
                       level = 0.95, mean.response = FALSE, data, boot = FALSE,
                       type = c("parametric", "nonparametric"), nsim = 1, 
                       seed = NULL, progress = FALSE, lower, upper, 
                       tol = .Machine$double.eps^0.25, maxiter = 1000, 
                       adjust = c("none", "Bonferroni"), k, ...) {
  
  ## No support for the Golub-Pereyra algorithm for partially linear 
  ## least-squares models
  if (object$call$algorithm == "plinear") {
    stop(paste("The Golub-Pereyra algorithm for partially linear least-squares 
               models is currently not supported."))
  }
  
  ## Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (missing(lower)) lower <- min(.data[, xname])  # lower limit default
  if (missing(upper)) upper <- max(.data[, xname])  # upper limit default
  
  ## Set up for inverse estimation
  m <- length(y0)  # number of unknowns 
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  eta <- mean(y0)  # mean response
  n <- length(resid(object))  # sample size
  p <- length(coef(object))  # number of parameters
  var.pooled <- Sigma(object)^2  # residual variance
  
  ## Critical value for confidence interval computations
  adjust <- match.arg(adjust)
  crit <- if (adjust == "Bonferroni" && m == 1) {
            qt((level + 2*k - 1) / (2*k), n+m-p-1)  # Bonferroni adjustment
          } else {
            qt((level + 1) / 2, n+m-p-1)  # no adjustment
          }

  ## Calculate point estimate by inverting fitted model
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(x, xname)) - eta
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
  
  ## Parametric bootstrap ------------------------------------------------------
  if (boot) {
    
    ## Sanity check
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    
    ## Set up progress bar (if requested)
    if (progress) { 
      pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    }
    
    ## Initialize random number generator
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    if (is.null(seed)) 
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    
    ## Simulate new response vectors
    ftd <- fitted(object)  # fitted values
    res <- residuals(object) # redisuals
    type <- match.arg(type)
    if (type == "parametric") {  
      ss <- simulate(object, nsim = nsim)
    } else {
      ss <- replicate(nsim, ftd + sample(res, replace = TRUE), simplify = FALSE)
    }
    
    ## Function to calculate inverse estimate
    x0Fun <- function(i) {
      
      ## Update model using simulated response data
      boot.data <- eval(object$call$data)  # copy data
      boot.data[, yname] <- ss[[i]]  # simulated response vector
      boot.object <- tryCatch(update(object, data = boot.data),
                              error = function(e) NULL)
      
      ## If updating the model fails, then return value is NA
      if (is.null(boot.object)) {
        ret <- NA
      } else {
        
        ## Simulate new response (different from simulated response vector)
        if (mean.response) {  # regulation
          y0.star <- y0  # hold constant in bootstrap replications
        } else {  # calibration
          if (type == "parametric") {
            y0.star <- y0 + rnorm(length(y0), sd = Sigma(object))
          } else {
            y0.star <- y0 + sample(res, size = length(y0), replace = TRUE)
          }
        }
        
        ## Calculate point estimate
        ret <- tryCatch(uniroot(function(x) {
            predict(boot.object, newdata = makeData(x, xname)) - mean(y0.star)
          }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
          error = function(e) NA)
      }
      
      ## Update progress bar
      if (progress) { 
        setTxtProgressBar(pb, i) 
      }
      
      ## Return estimate
      ret
      
    }
    
    ## Calculate bootstrap replicates
    x0.star <- sapply(seq_len(nsim), x0Fun)
    
    ## Check for errors and return the runs that did not fail
    if (AnyNA(x0.star)) {
      num_fail <- sum(is.na(x0.star))
      warning("some bootstrap runs failed (", num_fail, "/", nsim, 
              ")")
      x0.star <- na.omit(x0.star)  # remove runs that failed
      attributes(x0.star) <- NULL  # remove attributes
    } else {
      num_fail <- NULL
    }

    ## Create and return a bootCal object (essentially a list that can later be 
    ## converted to an object of class boot)
    res <- list(original  = x0.est,  # original estimate
                bootreps  = x0.star,  # bootstrap replicates
                nsim      = nsim,  # number of simulations
                level     = level)  # desired confidence level
    class(res) = "bootCal"
    attr(res, "bootFail") <- num_fail
    return(res)

  }
  
  ## Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.est)
  
  ## Inversion interval --------------------------------------------------------
  if (interval == "inversion") {
    
    ## Inversion function
    inversionFun <- function(x) {
      pred <- predict2(object, newdata = makeData(x, xname)) 
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
        predict(object.copy, newdata = makeData(x, xname)) - z
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
                       level = 0.95, mean.response = FALSE, data, lower, upper, 
                       q1, q2, tol = .Machine$double.eps^0.25, maxiter = 1000, 
                       ...) 
{
  
  ## Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else object$data
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (missing(lower)) lower <- min(.data[, xname])  # lower limit default
  if (missing(upper)) upper <- max(.data[, xname])  # upper limit default
  
  ## Set up for inverse estimation
#   if(m > 1) stop("Only one response value allowed.")
  m <- length(y0)
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  eta <- mean(y0)
  if (m != 1) stop('Only a single unknown allowed for objects of class "lme".')
  N <- length(resid(object)) 
  p <- length(fixef(object))
#   res.var <- Sigma(object)^2  # residual variance
  
  ## Critical value. Oman (1998. pg. 445) suggests a t(1-alpha/2, N-1) dist.
  if (missing(q1)) q1 <- qnorm((1-level) / 2)
  if (missing(q2)) q2 <- qnorm((1+level) / 2)
  
  ## Calculate point estimate by inverting fitted model
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(x, xname), level = 0) - eta
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
  if (!mean.response) var.y0 <- varY(object, newdata = makeData(x0.est, xname))
  
  ## Inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    ## Inversion function
    inversionFun <- function(x, bound = c("lower", "upper")) {
      pred <- predict2(object, newdata = makeData(x, xname), se.fit = TRUE)
      denom <- if (mean.response) pred$se.fit else sqrt(var.y0 + pred$se.fit^2)
      bound <- match.arg(bound)
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
                          data = makeData(x, xname))
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

##' @rdname calibrate
##' @export
##' @method calibrate glm
invest.glm <- function(object, y0, interval = c("inversion", "Wald", "none"), 
                       level = 0.95, lower, upper, data, linkinv = FALSE,
                       tol = .Machine$double.eps^0.25, maxiter = 1000, ...) {
  
  ## NOTE: Currently, this function only works for the case 
  ##       mean.response = TRUE. 
  
  ## Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (missing(lower)) lower <- min(.data[, xname])  # lower limit default
  if (missing(upper)) upper <- max(.data[, xname])  # upper limit default
  
  ## Perform all calculation on the link scale and then convert back to response
  ## scale if linkinv = TRUE.
  trans <- if (linkinv) family(object)$linkinv else I
  
  ## Calculate point estimate by inverting fitted model
  eta <- family(object)$linkfun(y0)
  x0.est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(x, xname), type = "link") - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  ## Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.est)
  
  ## Critical value
  crit <- qnorm((level + 1) / 2)  # quantile from standard normal
  
  ## Inversion interval
  ##
  ## Based on exercise 5.31 on pg. 207 of Categorical Data Analysis (2nd ed.) by 
  ## Alan Agresti.
  if (interval == "inversion") { 
    
    covmat <- vcov(object)
    
    ## Inversion function
    inversionFun <- function(x) {
      pred <- predict(object, newdata = makeData(x, xname), se.fit = TRUE,
                      type = "link")
      ((eta - pred$fit) ^ 2) / (pred$se.fit ^ 2) - crit^2
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
    res <- list("estimate" = trans(x0.est), 
                "lower" = trans(lwr), 
                "upper" = trans(upr), 
                "interval" = interval)
    
  }
  
  ## Wald interval -------------------------------------------------------------
  if (interval == "Wald") { 
    
    ## Function of parameters whose gradient is required
    object.copy <- object # FIXME: Is a copy really needed?
    dmFun <- function(params) {
      object.copy$coefficients <- params
      z <- eta
      uniroot(function(x) { 
        predict(object.copy, newdata = makeData(x, xname), type = "link") - z
      }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root
    }
    
    ## Variance-covariane matrix
    params <- coef(object)
    covmat <- vcov(object)
    
    ## Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    ## Store results in a list
    res <- list("estimate" = trans(x0.est), 
                "lower" = trans(x0.est - crit * se), 
                "upper" = trans(x0.est + crit * se), 
                "se" = se,  ## FIXME: How should this get transformed?
                "interval" = interval)
  }
  
  ## Assign class label and return results
  class(res) <- "calibrate"
  res
  
}

##' @keywords internal
print.bootCal <- function(x, digits = 4, ...) { 
  ci <- round(quantile(x$bootreps, probs = c((1-x$level)/2, (1+x$level)/2)),
              digits = digits)
  names(ci) <- c("lower", "upper")
  op <- c("estimate" = x$original, 
          "se"       = sd(x$bootreps), 
          "bias"     = mean(x$bootreps) - x$original)
  print(round(op, digits = digits))
  cat("\n", "Percentile bootstrap interval: (", ci[1], ",", ci[2], ")", "\n")
  invisible(x)
} 

##' Plots of the Output of a Bootstrap Calibration Simulation
##' 
##' This takes a bootstrap calibration object and produces plots for the 
##' bootstrap replicates of the inverse estimate.
##' 
##' @rdname plot.bootCal
##' @export
##' @method plot bootCal
##' 
##' @param x An object that inherits from class \code{"bootCal"}.
##' @param ... Additional optional arguments. At present, no optional arguments 
##'            are used.
plot.bootCal <- function(x, ...) {
  
  t <- x$bootreps  # bootstrap replicates
  t0 <- x$original  # original estimate
  
  ## Calculate number of histogram breaks (better than default)
  if (!is.null(t0)) {
    nclass <- min(max(ceiling(length(t)/25), 10), 100)
    rg <- range(t)
    if (t0 < rg[1L]) 
      rg[1L] <- t0
    else if (t0 > rg[2L]) 
      rg[2L] <- t0
    rg <- rg + 0.05 * c(-1, 1) * diff(rg)
    lc <- diff(rg)/(nclass - 2)
    n1 <- ceiling((t0 - rg[1L])/lc)
    n2 <- ceiling((rg[2L] - t0)/lc)
    bks <- t0 + (-n1:n2) * lc
  }
  
  ## Plots
  par(mfrow = c(1, 2))
  hist(t, breaks = bks, probability = TRUE, main = "",
       xlab = "Bootstrap replicate")
  qqnorm(t, xlab = "Standard normal quantile", ylab = "Bootstrap quantile",
         main = "")
  qqline(t)
  par(mfrow = c(1, 1))
  invisible(x)
  
}