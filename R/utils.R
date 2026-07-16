#' @keywords internal
makeData <- function(x, label) {
  stats::setNames(data.frame(x), label)
}


#' @keywords internal
make_bootcal <- function(x0.est, x0.boot, nsim, level, call) {
  perc <- unname(stats::quantile(x0.boot,
                                 probs = c((1 - level) / 2, (1 + level) / 2)))
  # The first five components are printed positionally by print.invest().
  # The trailing t0/t/R/sim/call components make the result a valid "boot"
  # object (issue #32), so it works with, e.g., boot::boot.ci(). Both
  # boot.type options simulate responses from the fitted model (rather than
  # resampling cases), which is "parametric" in the boot package's taxonomy;
  # labeling it as such makes boot.ci() correctly refuse the bca/stud types,
  # which cannot be computed from model-based replicates.
  boo <- list("estimate" = x0.est,  # original estimate
              "lower"    = perc[1L],  # lower percentile
              "upper"    = perc[2L],  # upper percentile
              "se"       = stats::sd(x0.boot),  # standard error
              "bias"     = mean(x0.boot) - x0.est,  # estimated bias
              "bootreps" = x0.boot,  # bootstrap replicates
              "nsim"     = nsim,  # number of simulations
              "level"    = level,  # desired confidence level
              "interval" = "percentile",  # type of interval requested
              "t0"       = x0.est,
              "t"        = matrix(x0.boot, ncol = 1L),
              "R"        = length(x0.boot),
              "sim"      = "parametric",
              "call"     = call)

  # Assign number of failed bootstrap replications as an attribute
  attr(boo, "bootFail") <- attr(x0.boot, "bootFail")
  boo
}


#' @keywords internal
check_newdata_classes <- function(newdata, data) {
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
