investr
================================================================================

`investr` stands for **inverse estimation in R**. Inverse 
estimation, also referred to as the calibration problem, is a classical and well-known problem in regression. In simple terms, it involves the use of an observed value of the response (or specified value of the mean response) to make inference on the corresponding unknown value of the explanatory variable. 

The package is currently listed on CRAN and can easily be installed:
```S
  ## Install from CRAN
  install.packages("investr", dep = TRUE)
```
The package is also currently listed in the ChemPhys task view (http://cran.r-project.org/web/views/ChemPhys.html) --- a collection of R packages useful for analyzing data from chemistry and physics experiments. These packages can all be installed at once (including `investr`) using the `ctv` package (Zeileis, 2005):
```S
  ## Install the ChemPhys task view
  install.packages("ctv")
  library("ctv")
  install.views("ChemPhys")
```
The development version can easily be installed using the the `devtools` package (Wickham and Chang, 2013):
```S
  ## Install development version from GitHub
  install.packages("devtools")
  library(devtools)
  install_github(repo = "w108bmg/investr")
```
To report bugs or issues, contact the main author directly or submit them to https://github.com/w108bmg/investr/issues. 

There are currently three main funtions in the `investr` package:
 * `plotFit`;
 * `calibrate`;
 * `invest`;

`plotFit`
--------------------------------------------------------------------------------
The function `plotFit` produces a scatterplot of the data with fitted regression curve and the option to add confidence/prediction bands for the response (pointwise or adjusted). Currently, it can only be used with single-predictor objects of class `lm` or `nls`; however, for objects of class `nls`, the confidence/prediction bands are based on the linear approximation and can be misleading (Bates and Watts, 1988, pg. 65).

### Usage
```S
## Load package
library(investr)

## A linear regression example
data(cars, package = "datasets")
library(splines)
cars.lm1 <- lm(dist ~ speed, data = cars)
cars.lm2 <- lm(dist ~ speed + I(speed^2), data = cars)
cars.lm3 <- lm(dist ~ poly(speed, degree = 3), data = cars)
cars.lm4 <- lm(dist ~ ns(speed, df = 3), data = cars)
par(mfrow = c(2, 2))
plotFit(cars.lm1, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
        main = "linear", shade = T)
plotFit(cars.lm2, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
        main = "quadratic", shade = T)
plotFit(cars.lm3, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
        main = "cubic", shade = T)
plotFit(cars.lm4, interval = "both", xlim = c(-10, 40), ylim = c(-50, 150), 
        main = "cubic spline", shade = T)
        
## A nonlinear regression example
par(mfrow = c(1, 1))
data(Puromycin, package = "datasets")
Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
                start = c(Vm = 200, K = 0.05))
plotFit(Puro.nls, interval = "prediction", pch = 19, shade = T, 
        col.pred = rgb(0, 1, 1, 0.4))

```

`calibrate`
--------------------------------------------------------------------------------
`calibrate` only operates on objects of class `lm` and can only be used with the simple linear regression model.

### Usage
```S
## Crystal growth data from Graybill & Iyer (1994)
fit <- lm(weight ~ time, data = crystal) 
plotFit(fit, interval = "confidence", shade = T, col.conf = "lightblue", 
        pch = 19)
(res <- calibrate(fit, y0 = 8, interval = "inversion", mean.response = T))
abline(h = 8, v = c(res$lower, res$estimate, res$upper), lty = 2)
```

`invest`
--------------------------------------------------------------------------------
For more complicated models (e.g., polynomial and nonlinear regression), use the more genaral `invest` function which carries out the computations numerically.

### Usage
```S
## Treatment group from Puromycin data frame
plotFit(fit, interval = "prediction", shade = T, pch = 19)
(res <- invest(Puro.nls, y0 = 150, interval = "inversion"))
abline(h = 150, v = c(res$lower, res$estimate, res$upper))
```

