##' Extract Data from an Object
##'
##' If present in the calling sequence used to produce \code{object}, the data 
##' frame used to fit the model is obtained.
##' @name getData
##'
##' @param object An object inheriting from class \code{lm}, representing a 
##'               fitted linear model.
##' @return if a \code{data} argument is present in the calling sequence that 
##'         produced \code{object}, the corresponding data frame (with 
##'         \code{na.action} and \code{subset} applied to it, if also present in 
##'         the call that produced \code{object}) is returned; else, \code{NULL}
##'          is returned.
##'
##' @importFrom nlme getData
##' @export getData
##' @method getData lm
##' @export
getData.lm <- function(object) {
  mCall <- object$call
  data <- eval(mCall$data)
  if (is.null(data)) return(data)
  ## Handle missing data
  naAct <- object[["na.action"]]
  if (!is.null(naAct)) {
    ## Guessing here: known values (omit, exclude) work.
    data <- if (inherits(naAct, "omit")) data[-naAct, ]
    else if (inherits(naAct, "exclude")) data
    else eval(mCall$na.action)(data)
  }
  ## Handle subsetted data
  subset <- mCall$subset
  if (!is.null(subset)) {
    subset <- eval(asOneSidedFormula(subset)[[2]], data)
    data <- data[subset, ]
  }
  data
}

##' Extract Variables Information from a Fitted Model
##' 
##' Extract the variables used to fit a particular model.
##' 
##' @rdname getVars
##' @export
##' 
##' @param object A fitted model object.
##' 
##' @return A list containing the following components:
##' \describe{
##'   \item{\code{x.names}}{The names of the predictor variables.}
##'   \item{\code{x.names}}{The names of the response variables.}
##'   \item{\code{x}}{The values of the predictor variables.}
##'   \item{\code{y}}{The values of the response variables.}
##'   \item{\code{x.dim}}{The number of predictor variables used.}
##'   \item{\code{y.dim}}{The number of response variables used.}
##' }
# 
# @examples
# data(cars, package = "datasets")
# cars.lm <- lm(dist ~ speed + I(speed^2), data = cars)
# getVars(cars.lm)
getVars <- function(object) {
  data <- getData(object)  # extract data from object
  x.names <- intersect(all.vars(formula(object)[[3]]), colnames(data))  
  y.names <- all.vars(formula(object)[[2]])
  x <- data[, x.names]  # extract predictor columns
  y <- data[, y.names]  # extract response columns
  list(x.names = x.names, y.names = y.names, x = x, y = y,
       x.dim = length(x.names), y.dim = length(y.names))
}

##' Extract residual standard error
##' 
##' Extract residual standard error from a fitted model. (For internal use 
##' only.)
##' 
##' @keywords internal
Sigma <- function(object, ...) {
  UseMethod("Sigma")
} 
Sigma.lm <- function(object, ...) summary(object)$sigma
Sigma.nls <- function(object, ...) summary(object)$sigma
Sigma.lme <- function(object, ...) object$sigma

##' Make new data frame
##' 
##' Create a new data frame from a specified x value that has the same structure 
##' as the data frame used to create \code{object}. 
##' 
##' @keywords internal
makeData <- function(object, x) {
  
  ## FIXME: What if object$call$data is NULL?
  .data <- if (inherits(object, "lme")) {
    getData(object) # object$data
  } else {
    eval(object$call$data, sys.frame())
  }
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data))
  if (length(xname) != 1) stop("Only a single predictor variable is allowed.")
  newdata <- data.frame(x)
  names(newdata) <- xname
  newdata
  
}

##' Construct design matrix for random effects
##'
##' Create a random effects design matrix from \code{newdata} based on a fitted 
##' model. (For internal use only.)
##' 
##' @rdname makeZ
##' @keywords internal
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
##' Create a fixed effects design matrix from \code{newdata} based on a fitted 
##' model. (For internal use only.)
##' 
##' @keywords internal
makeX <- function(object, newdata) {
  model.matrix(eval(object$call$fixed)[-2], data = newdata)
}

##' Evaluate response variance
##'
##' Evaluate response variance at a given value of the predictor variable. (For 
##' internal use only.)
##' 
##' @keywords internal
varY <- function(object, newdata) {
  Z <- makeZ(object, newdata)  # random effects design matrix
  G <- getVarCov(object)  # random effects variance covariance matrix
  as.numeric(Z %*% G %*% t(Z) + Sigma(object)^2)  # ZGZ' + (sigma^2)I
}

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
  
  ## TODO:
  ##   (1) How should missing values be handled?
  
  ## Extract data, variables, etc.
  if (missing(newdata)) {
    d <- eval(object$call$data, envir = parent.frame())
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
  alpha <- 1 - level  # FIXME: Could just use (level+1)/2 instead of 1-alpha/2
                      #        and (1-level)/2 instead of alpha/2.
  n <- length(resid(object))
  p <- length(coef(object))
  
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
             sqrt(p*qf(1 - alpha, p, pred$df))  # Working-Hotelling band
           } else {
             sqrt(k*qf(1 - alpha, k, pred$df))  # need k for prediction
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
      
      lwr <- pred$fit - w * sqrt(Sigma(object)^2 + pred$se.fit^2)
      upr <- pred$fit + w * sqrt(Sigma(object)^2 + pred$se.fit^2)

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
  ##   (1) Add option se.fit
  ##   (2) How should missing values be handled?
  
  ## Extract data, variables, etc.
  if (missing(newdata)) {
    .data <- eval(if("data" %in% names(object)) object$data else object$call$data,
                  envir = parent.frame())
    if (!is.null(object$call$subset)) {
      .data <- .data[with(.data, eval(object$call$subset)), ]  # subset data
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
    assign(param.names[i], coef(object)[i])  # FIXME: Should assign be used here?
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
  se.fit <- sqrt(Sigma(object)^2 * v0)
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
      
      lwr <- pred$fit - w * sqrt(Sigma(object)^2 + pred$se.fit^2)
      upr <- pred$fit + w * sqrt(Sigma(object)^2 + pred$se.fit^2)
      
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
  
  ## TODO:
  ##   (1) Check output from se.fit = TRUE against SAS
  
  if (missing(newdata)) newdata <- object$data
  fit <- predict(object, newdata = newdata, level = 0)
  if (se.fit) {
    X <- makeX(object, makeData(object, newdata))
    se.fit <- sqrt(diag(X %*% vcov(object) %*% t(X)))
    list(fit = fit, se.fit = se.fit)
  } else fit
  
}

##' Crystal weight data
##' 
##' The data give the growing time and final weight of crystals.
##' 
##' \itemize{
##'   \item time time taken to grow (hours)  
##'   \item weight final weight of the crystal (grams) 
##' }
##' @docType data
##' @keywords datasets
##' @format A data frame with 14 rows and 2 variables
##' @name crystal  
##' @references
##' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
##' Applications. Belmont, Calif: Duxbury Press, 1994.                      
NULL

##' Concentrations of arsenic in water samples
##' 
##' The data give the actual and measures concentration of arsenic present in 
##' water samples.
##' 
##' \itemize{
##'   \item actual True amount of arsenic present
##'   \item measured Measured amount of arsenic present 
##' }
##' @docType data
##' @keywords datasets
##' @format A data frame with 32 rows and 2 variables
##' @name arsenic  
##' @references 
##' Graybill, F. A., and Iyer, H. K. Regression analysis: Concepts and 
##' Applications. Belmont, Calif: Duxbury Press, 1994.
NULL 
