#' Predictions from a Fitted Model
#'
#' Generic prediction method for various types of fitted models. \code{predFit} 
#' can be used to obtain standard errors of fitted values and 
#' adjusted/unadjusted confidence/prediction intervals for objects of class 
#' \code{"lm"}, \code{"nls"}, and \code{"glm"}. 
#' 
#' @param object An object that inherits from class \code{"lm"}, \code{"glm"},
#'   \code{"nls"}, or \code{"lme"}.
#' @param newdata An optional data frame in which to look for variables with 
#'   which to predict. If omitted, the fitted values are used.    
#' @param type Character string specifying the type of prediction. Current 
#'   options are \code{type = "link"} (the default) and 
#'   \code{type = "response"}.
#' @param se.fit A logical vaue indicating if standard errors are required. 
#'   Default is \code{FALSE}.
#' @param interval Type of interval to be calculated. Can be one of "none" 
#'   (default), "confidence", or "prediction". Default is \code{"none"}.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for 
#'   the intervals (if any) to be calculated. Default is \code{0.95}.
#' @param adjust A logical value indicating if an adjustment should be made to
#'   the critical value used in calculating the confidence interval. This is 
#'   useful for when the calibration curve is to be used multiple, say k, times.
#'   Default is \code{FALSE}.
#' @param k The number times the calibration curve is to be used for computing 
#'   a confidence/prediction interval. Only needed when 
#'   \code{adjust = "Bonferroni"}.
#' @param ... Additional optional arguments. At present, no optional arguments 
#'   are used.
#'  
#' @returns If \code{se.fit = FALSE}, then \code{predFit()} returns a vector of 
#' predictions or a matrix of predictions and bounds with column names 
#' \code{fit}, \code{lwr}, and \code{upr} if \code{interval} is not 
#' \code{"none"}. (This function is more so meant for internal use.)
#' 
#' If \code{se.fit = TRUE}, then a list with the following components is 
#' returned:
#' 
#' \itemize{
#'   
#'   \item \code{fit} a vector or matrix as described above;
#'   
#'   \item \code{se.fit} a vector containing the standard errors of the 
#'   predicted means;
#'   
#'   \item \code{residual.scale} the residual standard deviations;
#'   
#'   \item \code{df} the residual degrees of freedom.
#'   
#' }
#'   
#' @details 
#' Confidence and prediction intervals for linear models (i.e., \code{"lm"} 
#' objects) are obtained according to the usual formulas. Nonlinear and 
#' generalized linear models (i.e., \code{"nls"} and \code{"glm"} objects), on 
#' the other hand, rely on Taylor-series approximations for the standard errors 
#' used in forming the intervals. Approximate standard errors for the fitted 
#' values in linear mixed-effects models (i.e., \code{"lme"} objects) can also 
#' be computed; however, these rely on the approximate variance-covariance 
#' matrix of the fixed-effects estimates and often under estimate the true
#' standard error. More accurate standard errors can be obtained using the 
#' parametric bootstrap; see \code{\link[lme4]{bootMer}} for details.
#' 
#' For linear and nonlinear models, it is possible to request \emph{adjusted}
#' confidence or prediction intervals using the Bonferroni method 
#' (\code{adjust = "Bonferroni"}) or Scheffe's method 
#' (\code{adjust = "Scheffe"}). For the Bonferroni adjustment, you must specify 
#' a value for \code{k}, the number of intervals for which the coverage is to 
#' hold simultaneously. For the Scheffe adjustment, specifying a value for 
#' \code{k} is only required when \code{interval = "prediction"}; if 
#' \code{interval = "confidence"}, \code{k} is set equal to \eqn{p}, the number 
#' of regression parameters. For example, calling \code{plotFit} on \code{"lm"}
#' objects with \code{interval = "confidence"} and \code{adjust = "Scheffe"} 
#' will plot the Working-Hotelling band.
#'   
#' @export
#' @examples 
#' # A linear regression example (see ?datasets::cars)
#' cars.lm <- lm(dist ~ speed + I(speed^2), data = cars)
#' predFit(cars.lm, interval = "confidence")
#' 
#' # A nonlinear least squares example (see ?datasets::Puromycin)
#' data(Puromycin, package = "datasets")
#' Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
#' Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
#'                 start = c(Vm = 200, K = 0.05))
#' conc <- seq(from = 0.02, to = 1.10, length = 101)
#' pred <- predFit(Puro.nls, newdata = data.frame(conc), interval = "prediction")
#' plot(Puromycin2, ylim = c(min(pred[, "lwr"]), max(pred[, "upr"])))
#' lines(conc, pred[, "fit"], lwd = 2)
#' lines(conc, pred[, "lwr"], lty = 2)
#' lines(conc, pred[, "upr"], lty = 2)
predFit <- function(object, ...) {
  UseMethod("predFit")
} 


#' @rdname predFit
#' @export
predFit.default <- function(object, ...) {
  stats::predict(object, ...)
}


#' @rdname predFit
#' @export
predFit.lm <- function(object, newdata, se.fit = FALSE,
                       interval = c("none", "confidence", "prediction"), 
                       level = 0.95, 
                       adjust = c("none", "Bonferroni", "Scheffe"), k, 
                       ...) {
  
  # Make sure se.fit is set to TRUE if intervals are requested
  interval <- match.arg(interval)
  compute.se.fit <- if (se.fit || (interval != "none")) TRUE else FALSE
  
  # Predicted values and, if requested, standard errors too
  if (missing(newdata)) {
    pred <- stats::predict(object, se.fit = compute.se.fit) 
  } else {
    # Suppress warning message printed when calling predict.lm with new data
    suppressWarnings(
      pred <- stats::predict(object, newdata = as.data.frame(newdata), 
                             se.fit = compute.se.fit)
    )
  } 

  # Compute results
  if (interval == "none") {
    
    # Vector of fitted/predicted values
    res <- pred
    
  } else { 
    
    # Critical value for interval computations
    adjust <- match.arg(adjust)
    crit <- if (adjust == "Bonferroni") {  # Bonferroni adjustment
      
      stats::qt((level + 2*k - 1) / (2*k), pred$df)
      
    } else if (adjust == "Scheffe") {  # Scheffe adjustment
      
      # Working-Hotelling band or adjusted prediction band for k predictions
      if (interval == "confidence") {
        p <- length(stats::coef(object))
        sqrt(p * stats::qf(level, p, pred$df))  # Working-Hotelling band
      } else {
        sqrt(k * stats::qf(level, k, pred$df))  # need k for prediction
      }  
      
    } else {   # no adjustment
      
      stats::qt((level + 1) / 2, pred$df)    
      
    }
    
    # Interval calculations
    if (interval == "confidence") {  # confidence interval for mean response
      lwr <- pred$fit - crit * pred$se.fit
      upr <- pred$fit + crit * pred$se.fit
    } else {  # prediction interval for individual response
      lwr <- pred$fit - crit * sqrt(stats::sigma(object)^2 + pred$se.fit^2)
      upr <- pred$fit + crit * sqrt(stats::sigma(object)^2 + pred$se.fit^2)
      # warning("predictions on current data refer to _future_ responses")
    }
    
    # Store results in a matrix
    res <- cbind("fit" = pred$fit, "lwr" = lwr, "upr" = upr)
    
  }
  
  # If standard errors of fitted values are requested, convert results to a list
  # and store additional information
  if (se.fit) {
    res <- list("fit" = if (interval != "none") res else pred$fit,
                "se.fit" = pred$se.fit,
                "df" = pred$df,
                "residual.scale" = pred$residual.scale)
  }
  
  # Return results
  return(res)
  
}


#' @rdname predFit
#' @export
predFit.glm <- function(object, newdata, type = c("link", "response"),
                        se.fit = FALSE, interval = c("none", "confidence"), 
                        level = 0.95, ...) {
  
  # Match arguments
  type <- match.arg(type)
  interval <- match.arg(interval)

  # Make sure se.fit is set to TRUE if intervals are requested
  if (se.fit || (interval != "none")) {
    compute.se.fit <- TRUE 
  } else {
    compute.se.fit <-  FALSE
  }
  
  # Predicted values and, if requested, standard errors too
  if (missing(newdata)) {
    pred <- stats::predict(object, se.fit = compute.se.fit) 
  } else {
    # Suppress warning message printed when calling predict.lm with new data
    suppressWarnings(
      pred <- stats::predict(object, newdata = as.data.frame(newdata), 
                             type = type, se.fit = compute.se.fit)
    )
  } 
  
  # Compute results
  if (interval == "none") {
    res <- pred
  } else { 
    res <- cbind("fit" = pred$fit, 
                 "lwr" = pred$fit - pred$se.fit * stats::qnorm((level + 1) / 2), 
                 "upr" = pred$fit + pred$se.fit * stats::qnorm((level + 1) / 2))
    if (type == "response") {
      res <- apply(res, MARGIN = 2, FUN = function(x) {
        stats::family(object)$linkinv(x)
      })
    }
  }
  
  # If standard errors of fitted values are requested, convert results to a list
  # and store addional information
  if (se.fit) {
    res <- list("fit" = if (interval != "none") res else pred$fit,
                "se.fit" = pred$se.fit,
                "df" = pred$df,
                "residual.scale" = pred$residual.scale)
  }
  
  # Return results
  return(res)
  
}


#' @rdname predFit
#' @export
predFit.nls <- function(object, newdata, se.fit = FALSE,
                        interval = c("none", "confidence", "prediction"), 
                        level = 0.95, 
                        adjust = c("none", "Bonferroni", "Scheffe"), k, 
                        ...) {
  
  # Match arguments
  interval <- match.arg(interval)
  adjust <- match.arg(adjust)
  
  # Make sure se.fit is set to TRUE if intervals are requested
  compute.se.fit <- if (se.fit || (interval != "none")) {
    TRUE
  } else {
    FALSE
  }
  
  # No support for the Golub-Pereyra algorithm for partially linear 
  # least-squares models
  if (interval != "none") {
    if (!is.null(object$call$algorithm) && object$call$algorithm == "plinear") {
      stop(paste0("The Golub-Pereyra algorithm for partially linear least-", 
                  "squares models is currently not supported."), call. = FALSE)
    }
  }

  # Prediction data
  newdata <- if (missing(newdata)) {
    eval(stats::getCall(object)$data, envir = parent.frame()) 
  } else {
    as.data.frame(newdata) 
  }
  if (is.null(newdata)) {
    stop("No data available for predictions.", call. = FALSE)
  }
  
  # Name of independent variable
  xname <- intersect(all.vars(stats::formula(object)[[3]]), colnames(newdata)) 
  
  # Predicted values
  pred <- object$m$predict(newdata)
  
  # Compute standard error
  if (compute.se.fit) {
    
    # Assign values to parameter names in current environment
    param.names <- names(stats::coef(object))  
    for (i in 1:length(param.names)) { 
      assign(param.names[i], stats::coef(object)[i])  
    }
    
    # Assign values to independent variable name
    assign(xname, newdata[, xname])  
    
    # Calculate gradient (numerically)
    form <- object$m$formula()
    rhs <- eval(form[[3]])
    if (is.null(attr(rhs, "gradient"))) {
      f0 <- attr(stats::numericDeriv(form[[3]], param.names), "gradient")
    } else {  # self start models should have gradient attribute
      f0 <- attr(rhs, "gradient")
    }
    
    # Calculate standard error
    R1 <- object$m$Rmat()
    # v0 <- diag(f0 %*% solve(t(R1) %*% R1) %*% t(f0))
    v0 <- diag(f0 %*% tcrossprod(solve(crossprod(R1)), f0))  # slightly faster
    se_fit <- sqrt(stats::sigma(object)^2 * v0)
    
  }
  
  # Compute results
  if (interval == "none") {
    
    # Vector of fitted/predicted values
    res <- pred    
    
  } else { 
    
    # Adjustment for simultaneous inference
    crit <- if (adjust == "Bonferroni") {  # Bonferroni adjustment 
      
      stats::qt((level + 2*k - 1) / (2*k), stats::df.residual(object))
      
    } else if (adjust == "Scheffe") {  # Scheffe adjustment
      
      if (interval == "confidence") {
        p <- length(stats::coef(object))  # number of regression parameters
        # sqrt(p * stats::qf((level + 1) / 2, p, stats::df.residual(object))) 
        sqrt(p * stats::qf(level, p, stats::df.residual(object))) 
      } else {
        # sqrt(k * stats::qf((level + 1) / 2, k, stats::df.residual(object))) 
        sqrt(k * stats::qf(level, k, stats::df.residual(object))) 
      }     
      
    } else {  # no adjustment   
      
      stats::qt((level + 1) / 2, stats::df.residual(object))   
      
    }
    
    # Interval calculations
    if (interval == "confidence") {  # confidence limits for mean response
      lwr <- pred - crit * se_fit  # lower limits
      upr <- pred + crit * se_fit  # upper limits
    } else {  # prediction limits for individual response
      lwr <- pred - crit * sqrt(stats::sigma(object)^2 + se_fit^2)  # lower limits
      upr <- pred + crit * sqrt(stats::sigma(object)^2 + se_fit^2)  # upper limits
    }
    
    # Store results in a matrix
    res <- cbind("fit" = pred, "lwr" = lwr, "upr" = upr)
    
  }
  
  # If standard errors of fitted values are requested, convert results to a list
  # and store additional information
  if (se.fit) {
    res <- list("fit" = if (interval != "none") res else pred, #res,
                "se.fit" = se_fit,
                "df" = stats::df.residual(object),
                "residual.scale" = stats::sigma(object))
  }
  
  # Return results
  return(res)
  
}


#' @rdname predFit
#' @export
predFit.lme <- function(object, newdata, se.fit = FALSE, ...) {
  
  # Prediction data
  newdata <- if (missing(newdata)) {
    object$data 
  } else {
    as.data.frame(newdata) 
  }  
  
  # Names of independent variables
  xname <- intersect(all.vars(stats::formula(object)[[3L]]), colnames(newdata)) 
  
  # Population predicted values
  pred <- stats::predict(object, newdata = newdata, level = 0)
  
  # Approximate standard errors of fitted values
  if (se.fit) {
    Xmat <- makeX(object, newdata)  # fixed-effects design matrix
    se <- sqrt(diag(Xmat %*% stats::vcov(object) %*% t(Xmat)))
    cbind("fit" = pred, "se.fit" = se)
  } else {
    pred
  }
  
}
