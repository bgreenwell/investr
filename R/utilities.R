#' Predict method for (Single-Regressor) Linear and Nonlinear Model Fits
#'
#' Convenience function to be called by \code{plotFit}. It is not for routine 
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
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  xx <- list(d[[ xvar ]])
  names(xx) <- xvar
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
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  xx <- d[[xvar]]
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
  assign(xvar, d[, xvar])
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

#' Plotting Confidence/Prediction Bands
#' 
#' Plots fitted model for an object of class \code{lm} or \code{nls}vwith the 
#' option of adding a confidence and/or prediction band. 
#'
#' @param object An object that inherits from class \code{lm} or \code{nls}.
#' @param interval A character string indicating if a prediction band, 
#'   confidence band, both, or none should be plotted.
#' @param level The desired confidence level.
#' @param adjust A character string indicating the type of adjustment (if any) 
#' to make to the confidence/prediction bands.
#' @param k An integer to be used in computing the critical value for the 
#' confidence/prediction bands. Only needed when \code{adjust = "Bonferroni"} or
#' when \code{adjust = "Scheffe"} and \code{interval = "prediction"}.
#' @param shade A logical value indicating if the band should be shaded.
#' @param extend.range A logical value indicating if the fitted regression line
#' and bands (if any) should extend to the edges of the plot. Default is 
#' \code{FALSE}.
#' @param col.conf Shade color for confidence band.
#' @param col.pred Shade color for prediction band.
#' @param col.fit The color to use for the fitted line.
#' @param border.conf The color to use for the confidence band border.
#' @param border.pred The color to use for the prediction band border. 
#' @param lty.conf Line type to use for confidence band border.
#' @param lty.pred Line type to use for prediction band border.
#' @param lty.fit Line type to use for the fitted regression line.
#' @param lwd.conf Line width to use for confidence band border.
#' @param lwd.pred Line width to use for prediction band border.
#' @param lwd.fit Line width to use for the fitted regression line.
#' @param n The number of predictor values at which to evaluate the fitted model
#' (larger implies a smoother plot).
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param xlim The x limits (x1, x2) of the plot.
#' @param ylim The y limits (y1, y2) of the plot. 
#' @param hide A logical value indicating if the fitted model should be 
#' plotted on top of the points (FALSE) or behind them (TRUE). Default is 
#' TRUE.
#' @param ... Additional optional arguments passed on to \code{plot}.
#' @rdname plotFit
#' @export
#' @note
#' By default, the plotted intervals are pointwise intervals. For simultaneous 
#' intervals use \code{adjust = "Bonferroni"} or \code{adjust = "Scheffe"}. For
#' the Bonferroni adjustment, you must specify value for \code{k}, the number
#' of intervals for which the coverage is to hold simultaneously. For the 
#' Scheffe adjustment, specifying a value for \code{k} is only required when
#' \code{interval = "prediction"}; if \code{interval = "confidence"}, \code{k} 
#' is set equal to \eqn{p}, the number of regression parameters. For example,
#' if \code{object} is a simple linear regression model, then calling 
#' \code{plotFit} with \code{interval = "confidence"} and 
#' \code{adjust = "Scheffe"} will plot the Working-Hotelling band.
#' 
#' Confidence/prediction bands for nonlinear regression (i.e., objects of class
#' \code{nls}) are based on a linear approximation as described in Bates & Watts 
#' (2007). This funtion was inpired by the \code{plotfit} function in the 
#' \code{nlstools} package.
#' @references
#' Bates, D. M., and Watts, D. G. Nonlinear Regression Analysis and its 
#' Applications. New York: Wiley, 2007.
#' 
#' F. Baty and M. L. Delignette-Muller (2012), nlstools: Tools for Nonlinear 
#' Regression Diagnostics.
#' @examples
#' \donttest{
#' ## A linear regression example
#' data(cars, package = "datasets")
#' library(splines)
#' cars.lm1 <- lm(dist ~ speed, data = cars)
#' cars.lm2 <- lm(dist ~ speed + I(speed^2), data = cars)
#' cars.lm3 <- lm(dist ~ speed + I(speed^2) + I(speed^3), data = cars)
#' cars.lm4 <- lm(dist ~ ns(speed, df = 3), data = cars)
#' par(mfrow = c(2, 2))
#' plotFit(cars.lm1, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
#'         main = "linear")
#' plotFit(cars.lm2, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
#'         main = "quadratic")
#' plotFit(cars.lm3, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
#'         main = "cubic")
#' plotFit(cars.lm4, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
#'         main = "cubic spline")
#'   
#' ## A nonlinear regression example
#' par(mfrow = c(1, 1))
#' library(RColorBrewer) # requires that RColorBrewer be installed
#' blues <- brewer.pal(9, "Blues")
#' data(Puromycin, package = "datasets")
#' Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
#' Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
#'                 start = c(Vm = 200, K = 0.05))
#' plotFit(Puro.nls, interval = "both", pch = 19, shade = TRUE, 
#'         col.conf = blues[4], col.pred = blues[2])
#' }     
plotFit <- function(object, ...) {
  UseMethod("plotFit")
} 

#' @rdname plotFit
#' @export
#' @method plotFit lm
plotFit.lm <- function(object,  
  interval = c("none", "both", "confidence", "prediction"), level = 0.95,
  adjust = c("none", "Bonferroni", "Scheffe"), k, ...,
  shade = FALSE, extend.range = FALSE, hide = TRUE,
  col.conf = if (shade) grey(0.7) else "black", 
  col.pred = if (shade) grey(0.9) else "black",  
  border.conf = col.conf, border.pred = col.pred, col.fit = "black", 
  lty.conf = if (shade) 1 else 2, lty.pred = if (shade) 1 else 3, lty.fit = 1, 
  lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, xlab, ylab, xlim, ylim = NULL)
{

  ## Extract data, variable names, etc.
  d <- eval(object$call$data, sys.frame())
  if (!is.null(object$call$subset)) {
    dsub <- with(d, eval(object$call$subset))
    d <- d[dsub, ]
  }
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (length(xvar) != 1) {
    stop("only a single predictor is allowed")
  }
  if (missing(xlim)) {
    xlim <- c(min(d[, xvar]), max(d[, xvar]))
  }
  xx <- if (extend.range) {
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xx) <- xvar
  if (missing(xlab)) xlab <- xvar
  if (missing(ylab)) ylab <- yvar
  interval = match.arg(interval)
  adjust <- match.arg(adjust)
  
  ## Maximum and minimum of fitted values
  if (interval == "none") {
    fitvals <- predict2(object, newdata = xx)$fit
    fit.ymin <- min(fitvals)
    fit.ymax <- max(fitvals)
  }
  
  ## Confidence interval for mean response
  if (interval == "confidence" || interval == "both") {
    conf <- predict2(object, newdata = xx, interval = "confidence", 
                     level = level, adjust = adjust, k = k)
    conf.lwr <- conf$lwr
    conf.upr <- conf$upr
    conf.ymin <- min(conf.lwr)
    conf.ymax <- max(conf.upr)
  }
  
  ## Prediction interval for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predict2(object, newdata = xx, interval = "prediction", 
                     level = level, adjust = adjust, k = k)
    pred.lwr <- pred$lwr
    pred.upr <- pred$upr
    pred.ymin <- min(pred.lwr)
    pred.ymax <- max(pred.upr)
  }
  
  ## Automatic limits for y-axis
  ylim <- if(interval == "prediction" || interval == "both") {
            c(min(c(pred.ymin, d[, yvar])), max(c(pred.ymax, d[, yvar])))
          } else if (interval == "confidence") {
            c(min(c(conf.ymin, d[, yvar])), max(c(conf.ymax, d[, yvar])))
          } else if (interval == "none") {
            c(min(c(fit.ymin, d[, yvar])), max(c(fit.ymax, d[, yvar])))
          } else {
            ylim
          }
  
  ## Plot data, fit, etc.
  if (hide) { ## Draw band behind points
    
    plot(d[, c(xvar, yvar)], xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         
      panel.first = {
        
        if (shade) {
          
          ## Draw (hidden) shaded prediction band
          if (interval == "prediction" || interval == "both") {
            polygon(c(xx[[1]], rev(xx[[1]])), c(pred.lwr, rev(pred.upr)), 
                    col = col.pred, border = border.pred, lty = lty.pred, 
                    lwd = lwd.pred)
          }
          
          ## Draw (hidden) shaded confidence band
          if (interval == "confidence" || interval == "both") {
            polygon(c(xx[[1]], rev(xx[[1]])), c(conf.lwr, rev(conf.upr)), 
                    col = col.conf, border = border.conf, lty = lty.conf, 
                    lwd = lwd.conf)
          }
          
        } else {
          
          ## Draw (hidden) unshaded prediction band
          if (interval == "prediction" || interval == "both") {
            lines(xx[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                  lwd = lwd.pred)
            lines(xx[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                  lwd = lwd.pred)
          }
          
          ## Draw (hidden) unshaded confidence band
          if (interval == "confidence" || interval == "both") {
            lines(xx[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                  lwd = lwd.conf)
            lines(xx[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                  lwd = lwd.conf)
          }
          
        }
        
        ## Draw (hidden) fitted response curve
        lines(xx[[1]], suppressWarnings(predict(object, newdata = xx)), 
              lty = lty.fit, lwd = lwd.fit, col = col.fit)
        
        }, ...)
    
  } else { ## Draw band on top of points
    
    plot(d[, c(xvar, yvar)], xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         
      panel.last = {
        
        if (shade) {
          
          ## Draw shaded prediction band
          if (interval == "prediction" || interval == "both") {
            polygon(c(xx[[1]], rev(xx[[1]])), c(pred.lwr, rev(pred.upr)), 
                    col = col.pred, border = border.pred, 
                    lty = lty.pred, lwd = lwd.pred)
          }
          
          ## Draw shaded confidence band
          if (interval == "confidence" || interval == "both") {
            polygon(c(xx[[1]], rev(xx[[1]])), c(conf.lwr, rev(conf.upr)), 
                    col = col.conf, border = border.conf, 
                    lty = lty.conf, lwd = lwd.conf)
          }
          
        } else {
          
          ## Draw unshaded prediction band
          if (interval == "prediction" || interval == "both") {
            lines(xx[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                  lwd = lwd.pred)
            lines(xx[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                 lwd = lwd.pred)
          }
          
          ## Draw unshaded confidence band
          if (interval == "confidence" || interval == "both") {
            lines(xx[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                  lwd = lwd.conf)
            lines(xx[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                  lwd = lwd.conf)
          }
          
        }
        
        ## Draw fitted response curve
        lines(xx[[1]], suppressWarnings(predict(object, newdata = xx)), 
              lty = lty.fit, lwd = lwd.fit, col = col.fit)
        
      }, ...)
    
  }
  
}

#' @rdname plotFit
#' @export
#' @method plotFit nls
plotFit.nls <- function(object, 
  interval = c("none", "both", "confidence", "prediction"), level = 0.95,
  adjust = c("none", "Bonferroni", "Scheffe"), k, ..., 
  shade = FALSE, extend.range = FALSE, hide = TRUE,
  col.conf = if (shade) grey(0.7) else "black", 
  col.pred = if (shade) grey(0.9) else "black",  
  border.conf = col.conf, border.pred = col.pred, col.fit = "black", 
  lty.conf = if (shade) 1 else 2, lty.pred = if (shade) 1 else 3, lty.fit = 1, 
  lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, xlab, ylab, xlim, ylim = NULL)
{
  
  ## Extract data, variable names, etc.
  d <- eval(object$call$data, sys.frame())
  if (!is.null(object$call$subset)) {
    dsub <- with(d, eval(object$call$subset))
    d <- d[dsub, ]
  }
  yvar <- all.vars(formula(object)[[2]])
  xvar <- intersect(all.vars(formula(object)[[3]]), colnames(d))
  if (length(xvar) != 1) {
    stop("only a single predictor is allowed")
  }
  if (missing(xlim)) {
    xlim <- c(min(d[, xvar]), max(d[, xvar]))
  }
  xx <- if (extend.range) {
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xx) <- xvar
  if (missing(xlab)) xlab <- xvar
  if (missing(ylab)) ylab <- yvar
  interval = match.arg(interval)
  adjust <- match.arg(adjust)

  ## Confidence interval for mean response
  if (interval == "confidence" || interval == "both") {
    conf <- predict2(object, newdata = xx, interval = "confidence", 
                     level = level, adjust = adjust, k = k)
    conf.lwr <- conf$lwr
    conf.upr <- conf$upr
    conf.ymin <- min(conf.lwr)
    conf.ymax <- max(conf.upr)
  }

  ## Prediction interval for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predict2(object, newdata = xx, interval = "prediction", 
                     level = level, adjust = adjust, k = k)
    pred.lwr <- pred$lwr
    pred.upr <- pred$upr
    pred.ymin <- min(pred.lwr)
    pred.ymax <- max(pred.upr)
  }
  
  ## Automatic limits for y-axis
  ylim <- if(interval == "prediction" || interval == "both") {
    c(min(c(pred.ymin, d[, yvar])), max(c(pred.ymax, d[, yvar])))
  } else if (interval == "confidence") {
    c(min(c(conf.ymin, d[, yvar])), max(c(conf.ymax, d[, yvar])))
  } else ylim

  ## Plot data, fit, etc.
  if (hide) { ## Draw band behind points
    
    plot(d[, c(xvar, yvar)], xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         
         panel.first = {
           
           if (shade) {
             
             ## Draw (hidden) shaded prediction band
             if (interval == "prediction" || interval == "both") {
               polygon(c(xx[[1]], rev(xx[[1]])), c(pred.lwr, rev(pred.upr)), 
                       col = col.pred, border = border.pred, lty = lty.pred, 
                       lwd = lwd.pred)
             }
             
             ## Draw (hidden) shaded confidence band
             if (interval == "confidence" || interval == "both") {
               polygon(c(xx[[1]], rev(xx[[1]])), c(conf.lwr, rev(conf.upr)), 
                       col = col.conf, border = border.conf, lty = lty.conf, 
                       lwd = lwd.conf)
             }
             
           } else {
             
             ## Draw (hidden) unshaded prediction band
             if (interval == "prediction" || interval == "both") {
               lines(xx[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                     lwd = lwd.pred)
               lines(xx[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                     lwd = lwd.pred)
             }
             
             ## Draw (hidden) unshaded confidence band
             if (interval == "confidence" || interval == "both") {
               lines(xx[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                     lwd = lwd.conf)
               lines(xx[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                     lwd = lwd.conf)
             }
             
           }
           
           ## Draw (hidden) fitted response curve
           lines(xx[[1]], suppressWarnings(predict(object, newdata = xx)), 
                 lty = lty.fit, lwd = lwd.fit, col = col.fit)
           
         }, ...)
    
  } else { ## Draw band on top of points
    
    plot(d[, c(xvar, yvar)], xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
         
         panel.last = {
           
           if (shade) {
             
             ## Draw shaded prediction band
             if (interval == "prediction" || interval == "both") {
               polygon(c(xx[[1]], rev(xx[[1]])), c(pred.lwr, rev(pred.upr)), 
                       col = col.pred, border = border.pred, 
                       lty = lty.pred, lwd = lwd.pred)
             }
             
             ## Draw shaded confidence band
             if (interval == "confidence" || interval == "both") {
               polygon(c(xx[[1]], rev(xx[[1]])), c(conf.lwr, rev(conf.upr)), 
                       col = col.conf, border = border.conf, 
                       lty = lty.conf, lwd = lwd.conf)
             }
             
           } else {
             
             ## Draw unshaded prediction band
             if (interval == "prediction" || interval == "both") {
               lines(xx[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                     lwd = lwd.pred)
               lines(xx[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                     lwd = lwd.pred)
             }
             
             ## Draw unshaded confidence band
             if (interval == "confidence" || interval == "both") {
               lines(xx[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                     lwd = lwd.conf)
               lines(xx[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                     lwd = lwd.conf)
             }
             
           }
           
           ## Draw fitted response curve
           lines(xx[[1]], suppressWarnings(predict(object, newdata = xx)), 
                 lty = lty.fit, lwd = lwd.fit, col = col.fit)
           
         }, ...)
    
  }

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
