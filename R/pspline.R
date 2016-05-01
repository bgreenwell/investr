#' Fitting Penalized Regression Splines
#' 
#' Fits a penalized regression spline to the supplied data using a linear mixed- 
#' effects model representation.
#' 
#' @param x A numeric vector giving the values of the explanatory variable.
#' @param y A numeric vector giving the values of the response variable.
#' @param degree Integer specifying the degree of the piecewise polynomial, 
#'   default is \code{1} for piecewise linear splines.
#' @param knots A numeric vector specifying the internal breakpoints that define 
#'   the spline.
#' @param num.knots Integer specifying the number of knots.
#' @param method Method used for comuting the variance components. Default is
#'   \code{"REML"}.
#' @param ... Additional optional arguments to be passed on to \code{nlme::lme}.
#' @return An object of class \code{pspline} (essentially a list) containing the 
#'   following components:
#'   \itemize{
#'     \item \code{x} The vector of predictor values.
#'     \item \code{y} The vector of response values.
#'     \item \code{xname}
#'     \item \code{yname}
#'     \item \code{fit} The vecotr of fitted values.
#'     \item \code{knots} Numeric vector containing the knot locations.
#'     \item \code{degree} The degree of the penalized regression spline.
#'     \item \code{beta.hat} Vector of estimated fixed effects.
#'     \item \code{u.hat} Vector of estimated random effects.
#'     \item \code{spar} The smoothing parameter.
#'     \item \code{var.rat} The ratio of variance components.
#'     \item \code{matrices} A list of key matrices.
#'     \item \code{var.components} Vector of variance components.
#'     \item \code{df.res} Residual degrees of freedom.
#'     \item \code{lmm} The original \code{"lme"} object.
#'     \item \code{call} The (matched) function call.
#'   }
#'   
#' @references 
#' Ruppert, D., Wand, M. & Carroll, R. (2003). Semiparametric regression. 
#' Cambridge New York: Cambridge University Press.
#' 
#' @rdname pspline
#' @aliases print.pspline
#' 
#' @export
#' 
#' @importFrom nlme lme pdIdent reStruct
#' @importFrom stats as.formula fitted quantile
#' 
#' @examples
#' # Proof of whisky stored in a charred oak barrel against time in years (see 
#' # Schoeneman, et.al. 1971)
#' whiskey <- data.frame(
#'   age = c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8),
#'   proof = c(104.6, 104.1, 104.4, 105.0, 106.0, 
#'             106.8, 107.7, 108.7, 110.6, 112.1)
#' )
#' plot(whiskey)
#' whiskey_ps <- with(whiskey, pspline(age, proof, degree = 1))
#' plotFit(whiskey_ps)
pspline <- function(x, ...) {
  UseMethod("pspline")
}


#' @rdname pspline
#' @export
pspline.default <- function(x, y, degree = 1, knots, num.knots, 
                            method = c("REML", "ML"), ...) {
  
  # Capture variable names
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))

  # Catch errors 
  if (!(degree %in% 1:3)) {
    stop("degree must be either 1, 2, or 3", call. = FALSE)
  }
    
  # Determine knot locations
  if (missing(knots)) {
    if (missing(num.knots)) {
      num.knots <- max(5, min(floor(length(unique(x)) / 4), 35))
    }
    knotseq <- seq(from = 0, to = 1, length = num.knots + 2)
    knots <- as.numeric(quantile(unique(x), knotseq)[-c(1, (num.knots + 2))])
  } else {
    num.knots <- length(knots)
  }
  
  # Spline basis matrix
  Z <- outer(x, knots, "-")
  Z <- (Z * (Z > 0)) ^ degree
  colnames(Z) <- paste("k", 1:ncol(Z), sep = "")
  
  # Linear mixed-effects model
  random.form <- as.formula(paste("~-1+", paste(colnames(Z), collapse = "+")))
  lmm <- lme(as.formula(paste("y ~ poly(x, degree =", degree, ", raw = TRUE)")), 
             data = data.frame(y = y, x = x, Z, 
                               subject = factor(rep(1, length(x)))), 
             method = match.arg(method),
             random = list(subject = pdIdent(random.form)), ...)                         
  beta.hat <- as.numeric(lmm$coef$fixed)  # fixed effects
  u.hat <- as.numeric(unlist(lmm$coef$random))  # random effects
  
  # Extract variance components and calculate smoothing parameter 
  var.e.hat <- lmm$sigma ^ 2  # residual variance
  var.u.hat <- as.numeric(var.e.hat * exp(2 * unlist(lmm$modelStruct)))
  var.rat <- var.e.hat / var.u.hat  # form ratio of variance components
  spar <- var.rat ^ (1 / (2 * degree))  # smoothing parameter
  
  # Calculate fitted values, key matrices, and residual degrees of freedom
  ftd <- as.numeric(fitted(lmm))  # fitted values
  X <- cbind(1, poly(x, degree = degree, raw = TRUE))  # design matrix
  C.mat <- cbind(X, Z)  # C matrix
  D.mat <- diag(c(rep(0, degree + 1), rep(1, num.knots))) 
  A.mat <- qr.solve(crossprod(C.mat) + var.rat * D.mat, tol = 1e-10)  
  S.mat <- C.mat %*% tcrossprod(A.mat, C.mat)  # smoother matrix
  df.res <- nrow(X) - 2 * sum(diag(S.mat)) + sum(diag(tcrossprod(S.mat)))
  
  # Return list of results
  cl <- match.call()
  cl[[1L]] <- as.name("pspline")
  res <- list("x" = x, 
              "y" = y, 
              "xname" = xname,
              "yname" = yname,
              "fit" = ftd, 
              "knots" = knots, 
              "degree" = degree, 
              "beta.hat" = beta.hat, 
              "u.hat" = u.hat, 
              "spar" = spar,
              "matrices" = list("A" = A.mat, "C" = C.mat, 
                                "D" = D.mat, "S" = S.mat),
              "var.rat" = var.rat,
              "var.components" = c("error" = var.e.hat, "u" = var.u.hat), 
              "df.res" = df.res, 
              "lmm" = lmm, 
              "call" = cl)
  class(res) <- "pspline"
  res
  
}


#' @rdname pspline
#' @export
pspline.data.frame <- function(x, ...) {
  res <- pspline.default(x = x[[1L]], y = x[[2L]], ...)
  cl <- match.call()
  cl[[1L]] <- as.name("pspline")
  res$call <- cl
  res
}


#' @rdname pspline
#' @export
pspline.matrix <- function(x, ...) {
  res <- pspline.default(x = x[, 1L], y = x[, 2L], ...)
  cl <- match.call()
  cl[[1L]] <- as.name("pspline")
  res$call <- cl
  res
}


#' @keywords internal
#' @export
print.pspline <- function (x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Smoothing parameter:\n")
  print.default(format(x$spar, digits = digits), print.gap = 2L, 
                quote = FALSE)
  cat("\n")
  if (length(x$beta.hat)) {
    cat("Fixed effects:\n")
    print.default(format(x$beta.hat, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No fixed effects\n")
  cat("\n")
  if (length(x$u.hat)) {
    cat("Random effects:\n")
    print.default(format(x$u.hat, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No random effects\n")
  cat("\n")
  invisible(x)
}


#' Predicting from Penalized Regression Spline Fits
#'
#' Generic prediction method for various types of fitted models. (For internal 
#' use only.)
#' 
#' @keywords internal
predict.pspline <- function(object, newdata, se.fit = FALSE, 
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95) {
  
  # Make sure se.fit is set to TRUE if intervals are requested
  interval <- match.arg(interval)
  compute.se.fit <- if (se.fit || (interval != "none")) TRUE else FALSE
  
  # Extract needed components from model
  degree <- object$degree
  beta.hat <- object$beta.hat
  u.hat <- object$u.hat
  
  # Create design matrices and compute fitted value(s)
  newx <- if (missing(newdata)) object$x else newdata
  Z <- outer(newx, object$knots, "-")
  Z <- (Z * (Z > 0)) ^ degree
  X <- cbind(1, poly(newx, degree = degree, raw = TRUE))
  f.hat <- as.numeric(X %*% beta.hat + Z %*% u.hat)
  
  res <- f.hat
  
  # Standard error of fitted value(s)
  if (compute.se.fit) {
    
    # Extract needed components from pspline object
    C.x <- cbind(X, Z)
    C.mat <- object$matrices$C
    A.mat <- object$matrices$A
      
    # Calculate standard error of fit
    CAC.mat <- C.x %*% A.mat %*% t(C.x)
    sigma.f <- sqrt(object$var.components[1L] * diag(CAC.mat))
        
  }
  
  # Confidence/prediction interval(s)
  if (interval == "none") {
    
    # Vector of fitted/predicted values
    res <- f.hat
    
  } else {
    
    # Standard error
    se <- if (interval == "confidence") {
      sigma.f
    } else {
      sqrt(object$var.components[1L] + sigma.f^2)
    }
    
    # Store results in a matrix
    res <- cbind("fit" = f.hat, 
                 "lwr" = f.hat - qt((1 + level) / 2, 
                                    df = df.residual(object)) * se,
                 "upr" = f.hat + qt((1 + level) / 2, 
                                    df = df.residual(object)) * se)
    
  }
  
  # If standard errors of fitted values are requested, convert results to a list
  # and store addional information
  if (se.fit) {
    res <- list("fit" = res,
                "se.fit" = sigma.f,
                "df" = df.residual(object),
                "residual.scale" = NULL)  # Sigma(object)
  }
  
  # Return results
  return(res)
  
}


# Returns the residual degrees-of-freedom extracted from a P-spline object
df.residual.pspline <- function(object) {
  object$df.res
}

