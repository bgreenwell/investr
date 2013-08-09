investr
================================================================================

This page is still a work in progress...

investr is an acronym for **inv**erse **est**imation in **R**. Inverse 
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
predictor for both calibration- and regulation-type problems in both linear and
nonlinear regression models. To install:

```{r}
## Install from CRAN
install.packages("investr", dependencies = TRUE)

## Try it out
library(investr)

## Crystal growth data from Graybill & Iyer (1994)
fit <- lm(weight ~ time, data = crystal) 
plotFit(fit, interval = "both", shade = T, pch = 19)

## Treatment group for Puromycin data frame
Puro.trt <- subset(Puromycin, state == "treated")
fit <- nls(rate ~ Vm * conc/(K + conc), data = Puro.trt, 
           start = c(Vm = 200, K = 0.05))
plotFit(fit, interval = "both", shade = T, pch = 19)
```

Alternatively, the development version can be installed using ```install_github()```
from the ```devtools``` package, but it may be easier to download the *.zip 
file copy the source code directly.

