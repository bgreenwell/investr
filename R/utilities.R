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
getData.lm <- function(object, envir) {
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
getVars <- function(object, return.data = FALSE) {
  data <- getData(object)  # extract data from object
  x.names <- intersect(all.vars(formula(object)[[3]]), colnames(data)) 
  y.names <- all.vars(formula(object)[[2]])
  x <- data[, x.names]  # extract predictor columns
  y <- data[, y.names]  # extract response columns
  res <- list(x.names = x.names, y.names = y.names, x = x, y = y,
              x.dim = length(x.names), y.dim = length(y.names))
  if (return.data) res$data <- data
  res
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
##' as the data frame used to create \code{object}. (For internal use only.)
##' 
##' @keywords internal
makeData <- function(object, x) {
  vars <- getVars(object)
#   if (is.null(dim(x))) x <- t(x)
#   if (ncol(x) != vars$x.dim) stop("Must supply values for each covariate.")
  if (vars$x.dim != 1) stop("Only objects with a single covariate are allowed.")
  setNames(data.frame(x), vars$x.names)
}

##' Construct design matrix for random effects
##'
##' Create a random effects design matrix from \code{newdata} based on a fitted 
##' model. (For internal use only.)
##' 
##' @rdname makeZ
##' @keywords internal
makeZ <- function(object, newdata) {
  Q <- object$dims$Q  # number of grouping levels
  mCall <- object$call  # list containing image of the nlme call
  fixed <- eval(eval(mCall$fixed)[-2])  # fixed effects formula
  reSt <- object$modelStruct$reStruct  # random effects structure
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
  Zmat <- makeZ(object, newdata)  # random effects design matrix
  Gmat <- getVarCov(object)  # random effects variance-covariance matrix
  var.y <- Zmat %*% Gmat %*% t(Zmat) + Sigma(object)^2  # ZGZ' + (sigma^2)I
  if (is.matrix(var.y)) unname(diag(var.y)) else var.y
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
  
  newdata <- if (missing(newdata)) getData(object) else as.data.frame(newdata) 
  xname <- getVars(object)$x.names  # extract covariate label
  n <- length(resid(object))  # sample size
  p <- length(coef(object))  # number of regression parameters
  
  ## No support for the Golub-Pereyra algorithm for partially linear 
  ## least-squares models
  if (object$call$algorithm == "plinear") {
    stop(paste("The Golub-Pereyra algorithm for partially linear least-squares 
               models is currently not supported."))
  }
  
  ## Calculate standard error of fitted values
  f0 <- attr(predict(object, newdata = newdata), "gradient")
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
  fit <- predict(object, newdata = newdata, level = 0)  # population predictions
  ## Approximate standard error of fitted values
  if (se.fit) {
    Xmat <- makeX(object, makeData(object, newdata))
    se.fit <- sqrt(diag(Xmat %*% vcov(object) %*% t(Xmat)))
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
