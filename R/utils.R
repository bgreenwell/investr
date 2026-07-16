#' @keywords internal
makeData <- function(x, label) {
  stats::setNames(data.frame(x), label)
}


#' @keywords internal
checkNewdataClasses <- function(newdata, data) {
  for (nm in intersect(names(newdata), names(data))) {
    old <- data[[nm]]
    new <- newdata[[nm]]
    if (is.factor(old)) {
      if (!is.factor(new) && !is.character(new)) {
        stop("Column \"", nm, "\" in newdata must be a factor or character ",
             "vector (the fitted model treats it as a factor).",
             call. = FALSE)
      }
      bad <- setdiff(as.character(new), levels(old))
      if (length(bad) > 0) {
        stop("Column \"", nm, "\" in newdata has level(s) not seen when ",
             "fitting the model: ", paste(sQuote(bad), collapse = ", "), ".",
             call. = FALSE)
      }
    } else if (is.numeric(old) && !is.numeric(new)) {
      stop("Column \"", nm, "\" in newdata must be numeric (the fitted ",
           "model treats it as numeric).", call. = FALSE)
    }
  }
  invisible(NULL)
}


#' @keywords internal
makeX <- function(object, newdata) {
  # Fixed effects model matrix
  stats::model.matrix(eval(object$call$fixed)[-2], data = newdata)
}


#' @keywords internal
makeZ <- function(object, newdata) {
  # Random effects model matrix
  Q <- object$dims$Q  # number of grouping levels
  mCall <- object$call  # list containing image of the nlme call
  fixed <- eval(eval(mCall$fixed)[-2])  # fixed effects formula
  reSt <- object$modelStruct$reStruct  # random effects structure
  mfArgs <- list(formula = nlme::asOneFormula(stats::formula(reSt), fixed),
                 data = newdata, na.action = stats::na.fail,
                 drop.unused.levels = TRUE)
  dataMix <- do.call("model.frame", mfArgs)
  stats::model.matrix(reSt, dataMix)
}


#' @keywords internal
sigma.lme <- function(object, ...) {
  object$sigma  # estimated standard deviation of the within-group error
}


#' @keywords internal
varY <- function(object, newdata) {
  
  # FIXME: What if object$call$correlation is not NULL?
  
  # Unconditional response variance: Var[Y] = Var[X*beta + Z*alpha + error]
  Zmat <- makeZ(object, newdata)  # random effects design matrix
  Gmat <- nlme::getVarCov(object)  # random effects variance-covariance matrix
  var.y <- Zmat %*% Gmat %*% t(Zmat) + stats::sigma(object) ^ 2  # ZGZ' + (sigma^2)I
  if (is.matrix(var.y)) {
    unname(diag(var.y)) 
  } else {
    var.y
  }
}
