investr
================================================================================

This page is still a work in progress...

investr stands for **inv**erse **est**imation in **R**. Inverse 
estimation, also referred to as the calibration problem, is in a sense the 
reverse of the prediction problem. More generally, consider two common uses
of a fitted regression model:

   * predict an individual response for a given value of the predictor
   (the prediction problem);

   * estimate the mean response for a given value of the predictor.

Loosely speaking, inverse estimation refers to the reverse of these. In other 
words, in regression, there is often a need to:

   * estimate the predictor value corresponding to an observed value of the 
   response (*calibration*);

   * estimate the predictor value corresponding to a specified value of the 
   mean response (*regulation*).

Calibration is a classical problem that has been thoroughly discussed in 
statistical literature, regulation, a related problem, is less well-known. The R
package `investr` was designed to estimate and make inferences on the unknown
value of the explanatory variable for both calibration- and regulation-type problems in both linear and nonlinear regression models. 

Installation
--------------------------------------------------------------------------------

```S
## Install from CRAN
install.packages("investr", dependencies = TRUE)
```

Usage
--------------------------------------------------------------------------------

```S
## Load package
library(investr)

## Crystal growth data from Graybill & Iyer (1994)
fit <- lm(weight ~ time, data = crystal) 
plotFit(fit, interval = "confidence", shade = T, col.conf = "lightblue", 
        pch = 19)
(res <- calibrate(fit, y0 = 8, interval = "inversion", mean.response = T))
abline(h = 8, v = c(res$lower, res$estimate, res$upper), lty = 2)

## Treatment group from Puromycin data frame
Puro.trt <- subset(Puromycin, state == "treated")
fit <- nls(rate ~ Vm * conc/(K + conc), data = Puro.trt, 
           start = c(Vm = 200, K = 0.05))
plotFit(fit, interval = "prediction", shade = T, pch = 19)
(res <- invest(fit, y0 = 150, interval = "inversion"))
abline(h = 150, v = c(res$lower, res$estimate, res$upper))
```

