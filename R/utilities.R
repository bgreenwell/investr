## Function to make new data frame for a given x
makeData <- function(object, x) {
  
  ## FIXME: What if object$call$data is NULL?
  d <- if (inherits(object, "lme")) {
    getData(object) # object$data
  } else {
    eval(object$call$data, sys.frame())
  }
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (length(xname) != 1) stop("Only a single predictor variable is allowed")
  newdata <- data.frame(x)
  names(newdata) <- xname
  newdata
  
}

##' Standard deviation function
##'
##' Function to extract the estimated error standard deviation
##' 
##' @rdname sigma
##' @keyword internal
sigma <- function(object, ...) {
  UseMethod("sigma")
} 
sigma.lm <- function(object, ...) summary(object)$sigma
sigma.nls <- function(object, ...) summary(object)$sigma
sigma.lme <- function(object, ...) object$sigma
sigma.merMod <- function(object, ...) getME(object, "sigma")

##' Construct design matrix for random effects
##'
##' \code{makeZ} creates a new design matrix for the random effects based on the
##'   data in \code{newdata}.
##' 
##' @rdname makeZ
##' @keyword internal
makeZ <- function(object, newdata) {
  Q <- object$dims$Q
  mCall <- object$call
  fixed <- eval(eval(mCall$fixed)[-2])
  reSt <- object$modelStruct$reStruct
  mfArgs <- list(formula = asOneFormula(formula(reSt), fixed),
                 data = newdata, na.action = na.fail,
                 drop.unused.levels = TRUE)
  dataMix <- do.call("model.frame", mfArgs)
  model.matrix(reSt, dataMix)
}

##' Construct design matrix for fixed effects
##'
##' \code{makeZ} creates a new design matrix for the fixed effects based on the
##'   data in \code{newdata}.
##' 
##' @rdname makeZ
##' @keyword internal
makeX <- function(object, newdata) {
  model.matrix(eval(object$call$fixed)[-2], data = newdata)
}

## Function to evaluate var(Y) = ZGZ' + sigma^2
varY <- function(object, newdata) {
  Z <- makeZ(object, newdata)
  G <- getVarCov(object)
  as.numeric(Z %*% G %*% t(Z) + sigma(object)^2)
}

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
  
  ## TODO:
  ##   (1) How should missing values be handled?
  
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
  
  ## TODO:
  ##   (1) Add option se.fit
  ##   (2) How should missing values be handled?
  
  ## Extract data, variables, etc.
  if (missing(newdata)) {
    .data <- eval(object$call$data, sys.frame())
    if (!is.null(object$call$subset)) {
      .data <- .data[with(.data, eval(object$call$subset)), ] # subset data
    }
  } else {
    .data <- as.data.frame(newdata) 
  }
  yname <- all.vars(formula(object)[[2]])
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data))
  xx <- .data[[xname]]       # predictor values
  alpha <- 1 - level         # alpha level
  n <- length(resid(object)) # sample size
  p <- length(coef(object))  # number of regression parameters
  
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
  assign(xname, .data[, xname])
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
  pred <- list(fit = object$m$predict(.data), se.fit = se.fit)  

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


##' @rdname predict2
##' @keywords internal
predict2.lme <- function(object, newdata, se.fit = FALSE, ...) {
  
  ## TODO:
  ##   (1) Check output from se.fit = TRUE against SAS
  
  if (missing(newdata)) newdata <- object$data
  fit <- predict(object, newdata = newdata, level = 0)
  if (se.fit) {
#     warning("se.fit ignores the variability of the estimated variance components!")
    X <- makeX(object, makeData(object, newdata))
    se.fit <- sqrt(diag(X %*% vcov(object) %*% t(X)))
    list(fit = fit, se.fit = se.fit)
  } else fit
  
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
