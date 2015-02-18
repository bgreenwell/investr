##' Plotting Confidence/Prediction Bands
##' 
##' Plots fitted model for an object of class \code{"lm"} or \code{"nls"} with 
##' the option of adding a confidence and/or prediction band. 
##'
##' @param object An object that inherits from class \code{"lm"} or \code{"nls"}.
##' @param interval A character string indicating if a prediction band, 
##'   confidence band, both, or none should be plotted.
##' @param level The desired confidence level.
##' @param data An optional data frame containing the variables in the model. 
##' @param adjust A character string indicating the type of adjustment (if any) 
##' to make to the confidence/prediction bands.
##' @param k An integer to be used in computing the critical value for the 
##' confidence/prediction bands. Only needed when \code{adjust = "Bonferroni"} or
##' when \code{adjust = "Scheffe"} and \code{interval = "prediction"}.
##' @param shade A logical value indicating if the band should be shaded.
##' @param extend.range A logical value indicating if the fitted regression line
##' and bands (if any) should extend to the edges of the plot. Default is 
##' \code{FALSE}.
##' @param col.conf Shade color for confidence band.
##' @param col.pred Shade color for prediction band.
##' @param col.fit The color to use for the fitted line.
##' @param border.conf The color to use for the confidence band border.
##' @param border.pred The color to use for the prediction band border. 
##' @param lty.conf Line type to use for confidence band border.
##' @param lty.pred Line type to use for prediction band border.
##' @param lty.fit Line type to use for the fitted regression line.
##' @param lwd.conf Line width to use for confidence band border.
##' @param lwd.pred Line width to use for prediction band border.
##' @param lwd.fit Line width to use for the fitted regression line.
##' @param n The number of predictor values at which to evaluate the fitted model
##' (larger implies a smoother plot).
##' @param xlab A title for the x axis.
##' @param ylab A title for the y axis.
##' @param xlim The x limits (x1, x2) of the plot.
##' @param ylim The y limits (y1, y2) of the plot. 
##' @param hide A logical value indicating if the fitted model should be 
##' plotted on top of the points (FALSE) or behind them (TRUE). Default is 
##' TRUE.
##' @param ... Additional optional arguments passed on to \code{plot}.
##' @rdname plotFit
##' @export
##' @note
##' By default, the plotted intervals are pointwise intervals. For simultaneous 
##' intervals use \code{adjust = "Bonferroni"} or \code{adjust = "Scheffe"}. For
##' the Bonferroni adjustment, you must specify a value for \code{k}, the number
##' of intervals for which the coverage is to hold simultaneously. For the 
##' Scheffe adjustment, specifying a value for \code{k} is only required when
##' \code{interval = "prediction"}; if \code{interval = "confidence"}, \code{k} 
##' is set equal to \eqn{p}, the number of regression parameters. For example,
##' if \code{object} is a simple linear regression model, then calling 
##' \code{plotFit} with \code{interval = "confidence"} and 
##' \code{adjust = "Scheffe"} will plot the Working-Hotelling band.
##' 
##' Confidence/prediction bands for nonlinear regression (i.e., objects of class
##' \code{nls}) are based on a linear approximation as described in Bates & Watts 
##' (2007). This funtion was inpired by the \code{\link[nlstools]{plotfit}} function
##' from the \code{nlstools} package.
##' @references
##' Bates, D. M., and Watts, D. G. (2007)
##' \emph{Nonlinear Regression Analysis and its Applications}. Wiley.
##' 
##' F. Baty and M. L. Delignette-Muller (2012), 
##' A Toolbox for Nonlinear Regression in R: The Package nlstools.
##' \emph{Journal of Statistical Software} \bold{(under revision)}.
##' @examples
##' ##
##' ## A nonlinear regression example
##' ##
##' data(Puromycin, package = "datasets")
##' Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
##' Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
##'                 start = c(Vm = 200, K = 0.05))
##' plotFit(Puro.nls, interval = "both", pch = 19, shade = TRUE, 
##'         col.conf = "darkgrey", col.pred = "lightgrey")  
plotFit <- function(object, ...) {
  UseMethod("plotFit")
} 

##' @rdname plotFit
##' @export
##' @method plotFit lm
plotFit.lm <- function(object, 
                       interval = c("none", "both", "confidence", "prediction"), 
                       level = 0.95, data,
                       adjust = c("none", "Bonferroni", "Scheffe"), k, ...,
                       shade = FALSE, extend.range = FALSE, hide = TRUE,
                       col.conf = if (shade) grey(0.7) else "black", 
                       col.pred = if (shade) grey(0.9) else "black",  
                       border.conf = col.conf, border.pred = col.pred, 
                       col.fit = "black", lty.conf = if (shade) 1 else 2, 
                       lty.pred = if (shade) 1 else 3, lty.fit = 1, 
                       lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, 
                       xlab, ylab, xlim, ylim) 
{
  
  ## Preliminary (extract data, variable names, etc.)
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  ## Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default plot domain
  xgrid <- if (extend.range) {  # the x values at which to evaluate
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis

  ## Maximum and minimum of fitted values
  interval = match.arg(interval)
  if (interval == "none") {
    fitvals <- predict2(object, newdata = xgrid)$fit
    fit.ymin <- min(fitvals)
    fit.ymax <- max(fitvals)
  }
  
  ## Confidence interval for mean response
  adjust <- match.arg(adjust)
  if (interval == "confidence" || interval == "both") {
    conf <- predict2(object, newdata = xgrid, interval = "confidence", 
                     level = level, adjust = adjust, k = k)
    conf.lwr <- conf$lwr
    conf.upr <- conf$upr
    conf.ymin <- min(conf.lwr)
    conf.ymax <- max(conf.upr)
  }
  
  ## Prediction interval for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predict2(object, newdata = xgrid, interval = "prediction", 
                     level = level, adjust = adjust, k = k)
    pred.lwr <- pred$lwr
    pred.upr <- pred$upr
    pred.ymin <- min(pred.lwr)
    pred.ymax <- max(pred.upr)
  }
  
  ## Automatic limits for y-axis
  if (missing(ylim)) {
    if(interval == "prediction" || interval == "both") {
      ylim <- c(min(c(pred.ymin, yvals)), max(c(pred.ymax, yvals)))
    }
    if (interval == "confidence") {
      ylim <- c(min(c(conf.ymin, yvals)), max(c(conf.ymax, yvals)))
    } 
    if (interval == "none") {
      ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))
    }
  }
  
  ## Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {  # draw points last
         if (shade) {
           ## Draw (hidden) shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           ## Draw (hidden) shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           ## Draw (hidden) unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           ## Draw (hidden) unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         ## Draw (hidden) fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL,
       panel.last = if(!hide) {  # draw points first
         if (shade) {
           ## Draw shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           ## Draw shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           ## Draw unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           ## Draw unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         ## Draw fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL, ...)
  
}

##' @rdname plotFit
##' @export
##' @method plotFit nls
plotFit.nls <- function(object, 
                        interval = c("none", "both", "confidence", "prediction"), 
                        level = 0.95, data,
                        adjust = c("none", "Bonferroni", "Scheffe"), k, ..., 
                        shade = FALSE, extend.range = FALSE, hide = TRUE,
                        col.conf = if (shade) grey(0.7) else "black", 
                        col.pred = if (shade) grey(0.9) else "black",  
                        border.conf = col.conf, border.pred = col.pred, 
                        col.fit = "black", lty.conf = if (shade) 1 else 2, 
                        lty.pred = if (shade) 1 else 3, lty.fit = 1, 
                        lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, 
                        xlab, ylab, xlim, ylim)
{
  
  ## Preliminary (extract variable information, etc.)
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  ## Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default plot domain
  xgrid <- if (extend.range) {  # set up plotting grid
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis
  
  ## Maximum and minimum of fitted values
  interval = match.arg(interval)
  if (interval == "none") {
    fitvals <- predict2(object, newdata = xgrid)$fit
    fit.ymin <- min(fitvals)
    fit.ymax <- max(fitvals)
  }
  
  ## Confidence interval for mean response
  adjust <- match.arg(adjust)
  if (interval == "confidence" || interval == "both") {
    conf <- predict2(object, newdata = xgrid, interval = "confidence", 
                     level = level, adjust = adjust, k = k)
    conf.lwr <- conf$lwr
    conf.upr <- conf$upr
    conf.ymin <- min(conf.lwr)
    conf.ymax <- max(conf.upr)
  }
  
  ## Prediction interval for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predict2(object, newdata = xgrid, interval = "prediction", 
                     level = level, adjust = adjust, k = k)
    pred.lwr <- pred$lwr
    pred.upr <- pred$upr
    pred.ymin <- min(pred.lwr)
    pred.ymax <- max(pred.upr)
  }
  
  ## Automatic limits for y-axis
  if (missing(ylim)) {
    if(interval == "prediction" || interval == "both") {
      ylim <- c(min(c(pred.ymin, yvals)), max(c(pred.ymax, yvals)))
    }
    if (interval == "confidence") {
      ylim <- c(min(c(conf.ymin, yvals)), max(c(conf.ymax, yvals)))
    } 
    if (interval == "none") {
      ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))
    }
  }
  
  ## Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {
         if (shade) {
           ## Draw (hidden) shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           ## Draw (hidden) shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           ## Draw (hidden) unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           ## Draw (hidden) unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         ## Draw (hidden) fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)
       } else NULL,
       panel.last = if(!hide) {
         if (shade) {
           ## Draw shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           ## Draw shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           ## Draw unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           ## Draw unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         ## Draw fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)
       } else NULL, ...)
  
}

##' @rdname plotFit
##' @export
##' @method plotFit lm
plotFit.glm <- function(object, link = FALSE,
                        data, ..., extend.range = FALSE, hide = TRUE, 
                        col.fit = "black", lty.fit = 1, lwd.fit = 1, n = 500, 
                        xlab, ylab, xlim, ylim) 
{
  
  ## Preliminary (extract data, variable names, etc.)
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  ## Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default plot domain
  xgrid <- if (extend.range) {  # the x values at which to evaluate
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis
  
  ## Maximum and minimum of fitted values
  fitvals <- predict(object, newdata = xgrid, 
                     type = if (link) "link" else "response")
  fit.ymin <- min(fitvals)
  fit.ymax <- max(fitvals)
    
  ## Automatic limits for y-axis
  if (missing(ylim)) {
      ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))
  }
  
  ## Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {  # draw points last
         ## Draw (hidden) fitted response curve
         lines(xgrid[[1]], fitvals, 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL,
       panel.last = if(!hide) {  # draw points first
         ## Draw fitted response curve
         lines(xgrid[[1]], fitvals, 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL, ...)
  
}