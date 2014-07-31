##' Predict method for (Single-Regressor) Linear, Nonlinear, and (Linear) Mixed
##' Model Fits
##'
##' Generic prediction method for various types of fitted models. (For internal 
##' use only.)
##' 
##' @keywords internal
predict2 <- function(object, ...) {
  UseMethod("predict2")
} 

##' @keywords internal
predict2.lm <- function(object, newdata, 
                        interval = c("none", "confidence", "prediction"), 
                        level = 0.95, 
                        adjust = c("none", "Bonferroni", "Scheffe"), k, 
                        ...) {
  
  newdata <- if (missing(newdata)) getData(object) else as.data.frame(newdata) 
  n <- length(resid(object))  # sample size
  p <- length(coef(object))  # number of regression coefficients
  pred <- suppressWarnings(predict(object, newdata = newdata, se.fit = TRUE))
  
  ## Compute results
  interval <- match.arg(interval)
  if (interval == "none") {
    res <- pred  
  } else { 
    ## Critical value for interval computations
    adjust <- match.arg(adjust)
    crit <- if (adjust == "Bonferroni") {
      qt((level + 2*k - 1)/(2*k), pred$df)
    } else if (adjust == "Scheffe") {
      if (interval == "confidence") {
        sqrt(p*qf(level, p, pred$df))  # Working-Hotelling band
      } else {
        sqrt(k*qf(level, k, pred$df))  # need k for prediction
      }     
    } else {      
      qt((level + 1)/2, pred$df)      
    }
    ## Calculate intervals
    if (interval == "confidence") {  # confidence interval for mean response
      lwr <- pred$fit - crit * pred$se.fit
      upr <- pred$fit + crit * pred$se.fit
    } else {  # prediction interval for individual response
      lwr <- pred$fit - crit * sqrt(Sigma(object)^2 + pred$se.fit^2)
      upr <- pred$fit + crit * sqrt(Sigma(object)^2 + pred$se.fit^2)
    }
    ## Store results in a list
    res <- list(fit = as.numeric(pred$fit), 
                lwr = as.numeric(lwr), 
                upr = as.numeric(upr))
  }
  
  ## Return list of results
  return(res)
  
}

##' @keywords internal
predict2.nls <- function(object, newdata, 
                         interval = c("none", "confidence", "prediction"), 
                         level = 0.95, 
                         adjust = c("none", "Bonferroni", "Scheffe"), k, 
                         ...) {
  
  ## TODO:
  ##  * Improve code for estimating standard error of fitted values.
  
  ## No support for the Golub-Pereyra algorithm for partially linear 
  ## least-squares models
  if (object$call$algorithm == "plinear") {
    stop(paste("The Golub-Pereyra algorithm for partially linear least-squares 
               models is currently not supported."))
  }
  
  newdata <- if (missing(newdata)) getData(object) else as.data.frame(newdata) 
  xname <- getVarInfo(object)$x.names  # extract covariate label
  n <- length(resid(object))  # sample size
  p <- length(coef(object))  # number of regression parameters
  
  ## Compute standard error
  param.names <- names(coef(object))  # FIXME: Are these always named? 
  for (i in 1:length(param.names)) { 
    assign(param.names[i], coef(object)[i])  # FIXME: Can we assign all at once?
  }
  assign(xname, newdata[, xname])  # FIXME: Why is this here?
  form <- object$m$formula()
  rhs <- eval(form[[3]])
  if (is.null(attr(rhs, "gradient"))) {
    f0 <- attr(numericDeriv(form[[3]], param.names), "gradient")
  } else {  # self start models should have gradient attribute
    f0 <- attr(rhs, "gradient")
  }
  R1 <- object$m$Rmat()
  v0 <- diag(f0 %*% solve(t(R1) %*% R1) %*% t(f0))
  se.fit <- sqrt(Sigma(object)^2 * v0)
  pred <- list(fit = object$m$predict(newdata), se.fit = se.fit) 
  
  ## Compute results
  interval <- match.arg(interval)
  if (interval == "none") {
    res <- pred                
  } else { 
    # Adjustment for simultaneous inference
    adjust <- match.arg(adjust)
    alpha <- 1 - level
    crit <- if (adjust == "Bonferroni") {
      qt(1 - alpha/(2*k), df.residual(object))
    } else if (adjust == "Scheffe") {
      if (interval == "confidence") {
        sqrt(p * qf(1 - alpha, p, df.residual(object))) 
      } else {
        sqrt(k * qf(1 - alpha, k, df.residual(object))) 
      }     
    } else {      
      qt(1 - alpha/2, df.residual(object))      
    }
    
    
    if (interval == "confidence") {  # confidence limits for mean response
      lwr <- pred$fit - crit * pred$se.fit  # lower limits
      upr <- pred$fit + crit * pred$se.fit  # upper limits
    } else {  # prediction limits for individual response
      lwr <- pred$fit - crit * sqrt(Sigma(object)^2 + pred$se.fit^2)  # lower limits
      upr <- pred$fit + crit * sqrt(Sigma(object)^2 + pred$se.fit^2)  # upper limits
    }
    
    ## Store results in a list
    res <- list(fit = as.numeric(pred$fit), 
                lwr = as.numeric(lwr), 
                upr = as.numeric(upr),
                se.fit = as.numeric(se.fit))
  }
  
  ## Return list of results
  return(res)
  
  }

##' @keywords internal
predict2.lme <- function(object, newdata, se.fit = FALSE, ...) {
  if (missing(newdata)) newdata <- getData(object)
  xname <- getVarInfo(object)$x.names
  fit <- predict(object, newdata = newdata, level = 0)  # population predictions
  ## Approximate standard error of fitted values
  if (se.fit) {
    Xmat <- makeX(object, makeData(x, xname))
    se.fit <- sqrt(diag(Xmat %*% vcov(object) %*% t(Xmat)))
    list(fit = fit, se.fit = se.fit)
  } else fit
}