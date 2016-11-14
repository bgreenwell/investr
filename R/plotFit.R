# Expression to be evaluated after the plot axes are set up but before any 
# plotting takes place.
panel.first <- quote({
  if (hide) {  # draw points last
    if (shade) {
      # Draw (hidden) shaded prediction band
      if (interval == "prediction" || interval == "both") {
        graphics::polygon(c(xgrid[[1L]], rev(xgrid[[1L]])), 
                          c(pred[, "lwr"], rev(pred[, "upr"])), 
                          col = col.pred, border = border.pred, lty = lty.pred, 
                          lwd = lwd.pred)
      }
      # Draw (hidden) shaded confidence band
      if (interval == "confidence" || interval == "both") {
        graphics::polygon(c(xgrid[[1L]], rev(xgrid[[1L]])), 
                          c(conf[, "lwr"], rev(conf[, "upr"])), 
                          col = col.conf, border = border.conf, lty = lty.conf, 
                          lwd = lwd.conf)
      }
    } else {
      # Draw (hidden) unshaded prediction band
      if (interval == "prediction" || interval == "both") {
        graphics::lines(xgrid[[1L]], pred[, "lwr"], col = col.pred, 
                        lty = lty.pred, 
                        lwd = lwd.pred)
        graphics::lines(xgrid[[1L]], pred[, "upr"], col = col.pred, 
                        lty = lty.pred, 
                        lwd = lwd.pred)
      }
      # Draw (hidden) unshaded confidence band
      if (interval == "confidence" || interval == "both") {
        graphics::lines(xgrid[[1L]], conf[, "lwr"], col = col.conf, 
                        lty = lty.conf, 
                        lwd = lwd.conf)
        graphics::lines(xgrid[[1L]], conf[, "upr"], col = col.conf, 
                        lty = lty.conf, 
                        lwd = lwd.conf)
      }
    }
    # Draw (hidden) fitted response curve
    graphics::lines(xgrid[[1L]], suppressWarnings(
      stats::predict(object, newdata = xgrid, type = type)), 
      lty = lty.fit, lwd = lwd.fit, col = col.fit)  
  } else {
    NULL
  }
})


# Expression to be evaluated after plotting has taken place but before the axes, 
# title and box are added.
panel.last <- quote({
  if(!hide) {  # draw points first
    if (shade) {
      # Draw shaded prediction band
      if (interval == "prediction" || interval == "both") {
        graphics::polygon(c(xgrid[[1L]], rev(xgrid[[1L]])), 
                          c(pred[, "lwr"], rev(pred[, "upr"])), 
                          col = col.pred, border = border.pred, lty = lty.pred, 
                          lwd = lwd.pred)
      }
      # Draw shaded confidence band
      if (interval == "confidence" || interval == "both") {
        graphics::polygon(c(xgrid[[1L]], rev(xgrid[[1L]])), 
                          c(conf[, "lwr"], rev(conf[, "upr"])), 
                          col = col.conf, border = border.conf, lty = lty.conf, 
                          lwd = lwd.conf)
      }
    } else {
      # Draw unshaded prediction band
      if (interval == "prediction" || interval == "both") {
        graphics::lines(xgrid[[1L]], pred[, "lwr"], col = col.pred, 
                        lty = lty.pred, 
                        lwd = lwd.pred)
        graphics::lines(xgrid[[1L]], pred[, "upr"], col = col.pred, 
                        lty = lty.pred, 
                        lwd = lwd.pred)
      }
      # Draw unshaded confidence band
      if (interval == "confidence" || interval == "both") {
        graphics::lines(xgrid[[1L]], conf[, "lwr"], col = col.conf, 
                        lty = lty.conf, 
                        lwd = lwd.conf)
        graphics::lines(xgrid[[1L]], conf[, "upr"], col = col.conf, 
                        lty = lty.conf, 
                        lwd = lwd.conf)
      }
    }
    # Draw fitted response curve
    graphics::lines(xgrid[[1L]], suppressWarnings(
      stats::predict(object, newdata = xgrid, type = type)), 
      lty = lty.fit, lwd = lwd.fit, col = col.fit)  
  } else {
    NULL
  }
})


#' Plotting Fitted Models
#' 
#' Generic function for plotting predictions from various types of fitted 
#' models. \code{plotFit} currently supports objects of class \code{"lm"}, 
#' \code{"glm"}, and \code{"nls"}. A default method also exists which may be 
#' used for plotting the fitted mean response from other model fits 
#' (e.g., \code{"lqs"} and \code{"rlm"} from the \code{MASS} package.
#'
#' @param object A fitted model object. Typically, an object that inherits from 
#'   class \code{"lm"}, \code{"glm"}, or \code{"nls"}, but others may work too).
#' @param type The type of prediction required. The default is on the scale of 
#'   the response variable; the alternative \code{"link"} is on the scale of the 
#'   linear predictor. This option is only used when plotting \code{"glm"} 
#'   objects.
#' @param interval A character string indicating if a prediction band, 
#'   confidence band, both, or none should be plotted.
#' @param level The desired confidence level.
#' @param data An optional data frame containing the variables in the model. 
#' @param adjust A character string indicating the type of adjustment (if any) 
#'   to make to the confidence/prediction bands.
#' @param k An integer to be used in computing the critical value for the 
#'   confidence/prediction bands. Only needed when \code{adjust = "Bonferroni"},
#'   or when \code{adjust = "Scheffe"} and \code{interval = "prediction"}.
#' @param shade A logical value indicating if the band should be shaded.
#' @param extend.range A logical value indicating if the fitted regression line
#'   and bands (if any) should extend to the edges of the plot. Default is 
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
#'   (larger gives a smoother plot).
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param xlim The x limits (x1, x2) of the plot.
#' @param ylim The y limits (y1, y2) of the plot. 
#' @param hide A logical value indicating if the fitted model should be plotted 
#'   on top of the points (\code{FALSE}) or behind them (\code{TRUE}). Default 
#'   is \code{TRUE}.
#' @param ... Additional optional arguments passed on to \code{plot}.
#' 
#' @seealso \code{\link[nlstools]{plotfit}}.
#' @rdname plotFit
#' @export
#' 
#' @note
#' By default, the plotted intervals are unadjusted (i.e., pointwise) intervals.
#' For simultaneous intervals use \code{adjust = "Bonferroni"} or 
#' \code{adjust = "Scheffe"}. For the Bonferroni adjustment, you must specify a 
#' value for \code{k}, the number of intervals for which the coverage is to hold simultaneously. For the 
#' Scheffe adjustment, specifying a value for \code{k} is only required when
#' \code{interval = "prediction"}; if \code{interval = "confidence"}, \code{k} 
#' is set equal to \eqn{p}, the number of regression parameters. For example,
#' if \code{object} is a simple linear regression model, then calling 
#' \code{plotFit} with \code{interval = "confidence"} and 
#' \code{adjust = "Scheffe"} will plot the Working-Hotelling band.
#' 
#' Confidence/prediction bands for nonlinear regression (i.e., objects of class
#' \code{nls}) are based on the linear approximation described in Bates & Watts 
#' (2007). T
#' 
#' @references
#' Bates, D. M., and Watts, D. G. (2007)
#' \emph{Nonlinear Regression Analysis and its Applications}. Wiley.
#' 
#' Florent Baty, Christian Ritz, Sandrine Charles, Martin Brutsche, 
#' Jean-Pierre Flandrois, Marie-Laure Delignette-Muller (2015). 
#' A Toolbox for Nonlinear Regression in R: The Package nlstools. 
#' \emph{Journal of Statistical Software}, \bold{66}(5), 1-21.
#' 
#' @examples
#' # A nonlinear least squares example (see ?datasets::Puromycin and 
#' # ?investr::predFit)
#' data(Puromycin, package = "datasets")
#' Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
#' Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
#'                 start = c(Vm = 200, K = 0.05))
#' plotFit(Puro.nls, interval = "both", pch = 19, shade = TRUE, 
#'         col.conf = "skyblue4", col.pred = "lightskyblue2")  
plotFit <- function(object, ...) {
  UseMethod("plotFit")
} 


#' @rdname plotFit
#' @export
plotFit.default <- function(object, type = c("response", "link"), 
                            interval = c("none", "both", "confidence", "prediction"), 
                            level = 0.95, data,
                            adjust = c("none", "Bonferroni", "Scheffe"), k, ..., 
                            shade = FALSE, extend.range = FALSE, hide = TRUE,
                            col.conf = if (shade) grDevices::grey(0.7) else "black", 
                            col.pred = if (shade) grDevices::grey(0.9) else "black",  
                            border.conf = col.conf, border.pred = col.pred, 
                            col.fit = "black", lty.conf = if (shade) 1 else 2, 
                            lty.pred = if (shade) 1 else 3, lty.fit = 1, 
                            lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, 
                            xlab, ylab, xlim, ylim) {
  
  # Match arguments
  type <- match.arg(type)
  interval <- match.arg(interval)
  adjust <- match.arg(adjust)
  
  # Try to catch errors
  if (interval != "none") {
    if (!inherits(object, c("glm", "lm", "nls"))) {
      stop(paste('Confidence and prediction bands can only be plotted for',
                 '"glm", "lm", or "nls" objects.'))
    }
    if (inherits(object, "glm") && interval == "prediction") {
      stop(paste('Prediction bands can only be plotted for "lm", or "nls"',
                 'objects.'))    
    }
  }
  if (inherits(object, c("lm", "nls")) && type != "response") {
    warning('Option "type" is ignored for "lm" and "nls" objects.')
  }
  
  # Extract data for plotting
  if (missing(data)) {
    data <- eval(stats::getCall(object)$data)
  }
  if (is.null(data)) {  # throw error if no data are found
    stop(paste("Could not find data to plot."))
  }
  
  # Dependent variable
  yname <- all.vars(stats::formula(object)[[2L]])
  if (inherits(object, "glm")) {
    # For binomial and quasibinomial families the response can also be specified 
    # as a factor (when the first level denotes failure and all others success) 
    # or as a two-column matrix with the columns giving the numbers of successes 
    # and failures.
    if (stats::family(object)$family %in% c("binomial", "quasibinomial")) {
      if (length(yname) == 1) {
        yvals <- with(data, eval(stats::formula(object)[[2L]]))
      } else {
        ymat <- data[, yname]
        yvals <- ymat[, 1L] / ymat[, 2L]
      }
    } else {
      if (length(yname) != 1) {
        stop("Only one dependent variable allowed.")
      }
      yvals <- with(data, eval(stats::formula(object)[[2L]]))
    }
    if (type == "link") {
      yvals <- stats::family(object)$linkfun(yvals)
    }
  } else {
    if (length(yname) != 1) {
      stop("Only one dependent variable allowed.")
    }
    yvals <- with(data, eval(stats::formula(object)[[2L]]))
  }
  
  # Independent variable
  xname <- intersect(all.vars(stats::formula(object)[[3L]]), colnames(data)) 
  if (length(xname) != 1) {
    stop("Only one independent variable allowed.")
  }
  xvals <- data[[xname]]
  
  # Check if user requested a logged x-axis
  dots <- list(...)
  logx <- FALSE
  if ("log" %in% names(dots)) {
    if (grepl("x", dots$log)) {
      logx <- TRUE
    }
  } 
  
  # Determine x-axis limits
  if (missing(xlim)) {
    xlim <- range(xvals)  # default x-axis limits
  }
  ulim <- if (logx) {
    log(xlim) 
  } else {
    xlim
  }
  if (extend.range) {
    ulim <- grDevices::extendrange(ulim)
  }
  
  # Grid of predictor values
  xgrid <- if (logx) {
    list(exp(seq(from = ulim[1L], to = ulim[2L], length = n)))  # logspace
  } else {
    list(seq(from = ulim[1L], to = ulim[2L], length = n))
  }
  names(xgrid) <- xname
  
  # Axis labels
  if (missing(xlab)) {
    xlab <- xname  # default x-axis label
  }
  if (missing(ylab)) {
    ylab <- if (inherits(object, "glm")) {
      paste0(deparse(stats::formula(object)[[2L]]), "(", type, " scale)")
    } else {
      deparse(stats::formula(object)[[2L]])  # default y-axis label
    }
  }
  
  # Fitted mean response values
  if (interval == "none") {
    fit <- if (inherits(object, "glm")) {
      unname(stats::predict(object, newdata = xgrid, type = type))
    } else {
      predFit(object, newdata = xgrid)
    }
  }
  
  # Confidence intervals for mean response
  if (interval == "confidence" || interval == "both") {
    if (inherits(object, "glm")) {
      conf <- stats::predict(object, newdata = xgrid, type = "link", 
                             se.fit = TRUE)
      conf <- cbind("fit" = conf$fit, 
                    "lwr" = conf$fit - conf$se.fit * 
                      stats::qnorm((level+1) / 2), 
                    "upr" = conf$fit + conf$se.fit * 
                      stats::qnorm((level+1) / 2))
    } else {
      conf <- predFit(object, newdata = xgrid, interval = "confidence", 
                      level = level, adjust = adjust, k = k)
    }
    if (inherits(object, "glm") && type == "response") {
      conf <- apply(conf, MARGIN = 2, FUN = function(x) {
        stats::family(object)$linkinv(x)
      })
    }
  }
  
  # Prediction intervals for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predFit(object, newdata = xgrid, interval = "prediction", 
                    level = level, adjust = adjust, k = k)
  }
  
  # Determine y-axis limits
  if (missing(ylim)) {
    ylim <- if (interval == "none") {
      x <- fit
      c(min(c(min(fit), yvals)), max(c(max(fit), yvals)))
    } else if (interval == "confidence") {
      c(min(c(min(conf[, "lwr"]), yvals)), max(c(max(conf[, "upr"]), yvals)))
    } else {
      c(min(c(min(pred[, "lwr"]), yvals)), max(c(max(pred[, "upr"]), yvals)))
    }
  }
  
  # Plot data, mean response, etc.
  graphics::plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, 
                 ylim = ylim, panel.first = eval(panel.first),
                 panel.last = eval(panel.last), ...)
  
}