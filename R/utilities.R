#' Predict method for (Single-Regressor) Linear and Nonlinear Model Fits
#'
#' use.
#' 
#' @rdname predict2
#' @keywords internal
predict2 <- function(object, ...) {
  UseMethod("predict2")
} 

#' @rdname predict2
#' @keywords internal
predict2.lm <- function(object, newdata, 
                        interval = c("none", "confidence", "prediction"), 
                        level = 0.95, 
                        adjust = c("none", "Bonferroni", "Scheffe"), k, 
                        ...) {
  
  ## Extract data, variables, etc.
  if (missing(newdata)) {
    d <- eval(object$call$data, sys.frame())
    if (!is.null(object$call$subset)) {
      dsub <- with(d, eval(object$call$subset))
      d <- d[dsub, ]
    }
  } else {
    d <- as.data.frame(newdata) 
  }
  yname <- all.vars(formula(object)[[2]])
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  xx <- list(d[[xname]])
  names(xx) <- xname
  alpha <- 1 - level
  n <- length(resid(object))
  p <- length(coef(object))
  alpha <- 1 - level
  
  ## FIXME: Why does this throw a warning?
  pred <- suppressWarnings(predict(object, newdata = xx, se.fit = TRUE))

  ## Compute results
  interval <- match.arg(interval)
  if (interval == "none") {
    
    ## Store results in a list
    res <- pred
  
  } else { 
    
    ## Adjustment for simultaneous inference.
    adjust <- match.arg(adjust)
    w <- if (adjust == "Bonferroni") {
           qt(1 - alpha/(2*k), pred$df)
         } else if (adjust == "Scheffe") {
           if (interval == "confidence") {
             sqrt(p*qf(1 - alpha, p, pred$df)) # Working-Hotelling band
           } else {
             sqrt(k*qf(1 - alpha, k, pred$df)) # need k for prediction
           }     
         } else {      
           qt(1 - alpha/2, pred$df)      
         }
    
    ## Confidence interval for mean response
    if (interval == "confidence") {

        lwr <- pred$fit - w * pred$se.fit
        upr <- pred$fit + w * pred$se.fit
    
    ## Prediction interval for individual response
    } else {
      
      lwr <- pred$fit - w * sqrt(summary(object)$sigma^2 + pred$se.fit^2)
      upr <- pred$fit + w * sqrt(summary(object)$sigma^2 + pred$se.fit^2)

    }
    
    ## Store results in a list
    res <- list(fit = as.numeric(pred$fit), 
                lwr = as.numeric(lwr), 
                upr = as.numeric(upr))
  }
    
  ## Return list of results
  return(res)
  
}

#' @rdname predict2
#' @keywords internal
predict2.nls <- function(object, newdata, 
                         interval = c("none", "confidence", "prediction"), 
                         level = 0.95, 
                         adjust = c("none", "Bonferroni", "Scheffe"), k, 
                         ...) {
  
  ## Extract data, variables, etc.
  if (missing(newdata)) {
    d <- eval(object$call$data, sys.frame())
    if (!is.null(object$call$subset)) {
      dsub <- with(d, eval(object$call$subset))
      d <- d[dsub, ]
    }
  } else {
    d <- as.data.frame(newdata) 
  }
  yname <- all.vars(formula(object)[[2]])
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  xx <- d[[xname]]
  alpha <- 1 - level
  n <- length(resid(object))
  p <- length(coef(object))
  alpha <- 1 - level
  
  ## Stop if object$call$algorithm == "plinear"
  if (object$call$algorithm == "plinear") {
    stop(paste("predict2.nls not yet implemented for 'nls' objects fit with the 
                plinear algorithm"))
  }
  
  ## Compute gradient
  param.names <- names(coef(object)) 
  for (i in 1:length(param.names)) { 
    assign(param.names[i], coef(object)[i]) 
  }
  assign(xname, d[, xname])
  form <- object$m$formula()
  rhs <- eval(form[[3]])
  if (is.null(attr(rhs, "gradient"))) {
    f0 <- attr(numericDeriv(form[[3]], param.names), "gradient")
  } else {
    f0 <- attr(rhs, "gradient")
  }
  
  ## Compute standard error of fitted values
  R1 <- object$m$Rmat()
  v0 <- diag(f0 %*% solve(t(R1) %*% R1) %*% t(f0))
  se.fit <- sqrt(summary(object)$sigma^2 * v0)
  pred <- list(fit = object$m$predict(d), 
               se.fit = se.fit)  

  ## Compute results
  interval <- match.arg(interval)
  if (interval == "none") {
    
    ## Store results in a list
    res <- pred
                
  } else { 
    
    # Adjustment for simultaneous inference
    adjust <- match.arg(adjust)
    w <- if (adjust == "Bonferroni") {
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
    
    ## Confidence interval for mean response
    if (interval == "confidence") {
      
      lwr <- pred$fit - w * pred$se.fit
      upr <- pred$fit + w * pred$se.fit
      
      ## Prediction interval for individual response
    } else {
      
      lwr <- pred$fit - w * sqrt(summary(object)$sigma^2 + pred$se.fit^2)
      upr <- pred$fit + w * sqrt(summary(object)$sigma^2 + pred$se.fit^2)
      
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

#' Crystal weight data
#' 
#' The data give the growing time and final weight of crystals.
#' 
#' \itemize{
#'   \item time time taken to grow (hours)  
#'   \item weight final weight of the crystal (grams) 
#' }
#' @docType data
#' @keywords datasets
#' @format A data frame with 14 rows and 2 variables
#' @name crystal  
#' @references
#' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
#' Applications. Belmont, Calif: Duxbury Press, 1994.                      
NULL

#' Concentrations of arsenic in water samples
#' 
#' The data give the actual and measures concentration of arsenic present in 
#' water samples.
#' 
#' \itemize{
#'   \item actual True amount of arsenic present
#'   \item measured Measured amount of arsenic present 
#' }
#' @docType data
#' @keywords datasets
#' @format A data frame with 32 rows and 2 variables
#' @name arsenic  
#' @references 
#' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
#' Applications. Belmont, Calif: Duxbury Press, 1994.
NULL 
