investr
================================================================================

`investr` stands for **inverse estimation in R**. Inverse 
estimation, also referred to as the calibration problem, is in a sense the 
reverse of the prediction problem. Inverse estimation is a classical and well-known problem in regression. In simple terms, it involves the use of an observed value of the response to make inference on the corresponding unknown value of the explanatory variable. To our knowledge, however, statistical software is somewhat lacking the capabilities for analyzing these types of problems. 

The `investr` package currently has three main functions:
  * `calibrate`;
  * `invest`;
  * `plotFit`.
  
`calibrate` operates on objects of class `lm` and can only be used for the simple linear regression model. For more complicated models (e.g., polynomial and nonlinear regression), use the `invest` function, which then calls `uniroot` from the `stats` package to do the computations numerically. The function `plotFit` produces a scatterplot of the data with fitted regression curve and the option to add confidence/prediction bands for the response (pointwise or adjusted). It can be used with single-predictor objects of class `lm` or `nls`; however, for objects of class `nls`, confidence/prediction bands are based on the linear approximation and can be misleading (Bates and Watts, 1988, pg. 65). The development version of `investr` can be found on GitHub at https://github.com/w108bmg/investr and can easily be installed using the the `devtools` package (Wickham and Chang, 2013):

```S
  ## Install development version from GitHub
  library(devtools)
  install_github(repo = "w108bmg/investr")
```

To report bugs or issues, contact the main author directly or submit them to https://github.com/w108bmg/investr/issues.

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

