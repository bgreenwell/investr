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
  # FIXME: Need to fix scoping issue caused by eval!
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
##' @rdname getVarInfo
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
getVarInfo <- function(object, data) {
  # FIXME: Need to fix scoping issue caused by eval within getData!
  if (missing(data)) {
    data <- try(getData(object), silent = TRUE)  # extract data from object
    if (inherits(data, "try-error")) {
      stop("Scoping issuse with generic function 'getData'.")
    }
  }
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
##' as the data frame used to create \code{object}. (For internal use only.)
##' 
##' @keywords internal
makeData <- function(x, label) {
  setNames(data.frame(x), label)
}

# makeData <- function(object, x) {
#   vars <- getVarInfo(object)
#   #   if (is.null(dim(x))) x <- t(x)
#   #   if (ncol(x) != vars$x.dim) stop("Must supply values for each covariate.")
#   if (vars$x.dim != 1) stop("Only objects with a single covariate are allowed.")
#   setNames(data.frame(x), vars$x.names)
# }

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