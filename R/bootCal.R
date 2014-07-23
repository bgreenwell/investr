source("/home/w108bmg/Desktop/Dropbox/devel/investr/R/utilities.R")
source("/home/w108bmg/Desktop/Dropbox/devel/investr/R/invest.R")
source("/home/w108bmg/Desktop/Dropbox/devel/investr/R/calibrate.R")
source("/home/w108bmg/Desktop/Dropbox/devel/investr/R/plotFit.R")
library(nlme)

bootCal <- function(object, nsim = 2) 
{ 
  data <- getData(object)
  t0 <- getVarInfo(object)$x[1]
  n <- length(res <- resid(object) - mean(resid(object)))  
  boot.data <- data.frame(data, res = res, fit = fitted(object))
  boot.fun <- function(.data, .ind) {
    d <- .data #.data[.ind, ]
    boot.object <- update(object, data = d)
    invest(boot.object, y0 = 600, interval = "none", lower = 0.01)
  }
  replicate(nsim, boot.fun(boot.data, sample(n, replace = TRUE)))
}


## Example
library(drc)
mod <- nls(weight ~ theta1/(1 + exp(theta2 + theta3*log(conc))),
           start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
           data = nasturtium)
bootCal(mod)

par(mfrow = c(1, 2))
plotFit(mod, interval = "both", shade = T, hide = T)
plotFit(mod, interval = "both", shade = T, hide = F)

invest(mod, y0 = 600)

f1 <- function(obj1) {
  f2 <- function(obj2) {
    plotFit(obj2)
  }
  f2(obj1)
}