# investr: Inverse Estimation in R

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/investr)](http://cran.r-project.org/package=investr)
[![Build Status](https://travis-ci.org/bgreenwell/investr.svg?branch=master)](https://travis-ci.org/bgreenwell/investr)
[![Coverage Status](https://img.shields.io/codecov/c/github/bgreenwell/investr.svg)](https://codecov.io/github/bgreenwell/investr?branch=master)
[![Downloads](http://cranlogs.r-pkg.org/badges/investr)](http://cranlogs.r-pkg.org/badges/investr)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/investr)](http://cranlogs.r-pkg.org/badges/grand-total/investr)

Inverse estimation, also referred to as the calibration problem, is a classical and well-known problem in regression. In simple terms, it involves the use of an observed value of the response (or specified value of the mean response) to make inference on the corresponding unknown value of the explanatory variable. 

A detailed introduction to investr has been published in The R Journal: "investr: An R Package for Inverse Estimation", http://journal.r-project.org/archive/2014-1/greenwell-kabban.pdf. You can track development at https://github.com/bgreenwell/investr. To report bugs or issues, contact the main author directly or submit them to https://github.com/bgreenwell/investr/issues. 

As of right now, `investr` supports (univariate) inverse estimation with objects of class:

* `lm` - linear models (multiple predictor variables allowed)
* `glm` - generalized linear models (multiple predictor variables allowed)
* `nls` - nonlinear least-squares models
* `lme` - linear mixed-effects models (fit using the `nlme` package)

# Installation
The package is currently listed on CRAN and can easily be installed:
```r
  # Install from CRAN
  install.packages("investr", dep = TRUE)
```
The package is also part of the [ChemPhys task view](http://cran.r-project.org/web/views/ChemPhys.html), a collection of R packages useful for analyzing data from chemistry and physics experiments. These packages can all be installed at once (including `investr`) using the `ctv` package (Zeileis, 2005):
```r
  # Install the ChemPhys task view
  install.packages("ctv")
  ctv::install.views("ChemPhys")
```

# Examples

## Dobson's Beetle Data

In binomial regression, the estimated lethal dose corresponding to a specific probability _p_ of death is often referred to as _LDp_. `invest` obtains an estimate of _LDp_ by inverting the fitted mean response on the link scale. Similarly, a confidence interval for _LDp_ can be obtained by inverting a confidence interval for the mean response on the link scale.
```r
library(investr)

# Dobson's beetle data
head(beetle)

# Binomial regression
binom_fit <- glm(cbind(y, n-y) ~ ldose, data = beetle, 
                 family = binomial(link = "cloglog"))
plotFit(binom_fit, lwd.fit = 2, cex = 1.2, pch = 21, bg = "lightskyblue", 
        lwd = 2, xlab = "Log dose", ylab = "Probability")

# Inverse estimation
invest(binom_fit, y0 = 0.5)   # median lethal dose
invest(binom_fit, y0 = 0.9)   # 90% lethal dose
invest(binom_fit, y0 = 0.99)  # 99% lethal dose


# estimate    lower    upper 
#   1.7788   1.7702   1.7862
```
![Alt text](https://raw.githubusercontent.com/bgreenwell/investr/master/beetle_plotFit.png)

To obtain an estimate of the standard error, we can use the Wald method:
```r
invest(binom_fit, y0 = 0.5, interval = "Wald")

# estimate    lower    upper       se 
#   1.7788   1.7709   1.7866   0.0040

# The MASS package function dose.p works too 
MASS::dose.p(binom_fit, p = 0.5)

#              Dose         SE
# p = 0.5: 1.778753 0.00400654
```

### Including a factor variable

Multiple predictor variables are allowed for objects of class `lm` and `gls`. 
For instance, the example from `?MASS::dose.p` can be re-created as follows:
```r

# Load package, assuming it is already installed
library(MASS)

# Data
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20 - numdead)
budworm <- data.frame(ldose, numdead, sex, SF)

# Logistic regression
budworm.lg0 <- glm(SF ~ sex + ldose - 1, family = binomial, data = budworm)

# Using dose.p function from package MASS
dose.p(budworm.lg0, cf = c(1, 3), p = 1/4)

#               Dose        SE
# p = 0.25: 2.231265 0.2499089

# Using invest function from package investr
invest(budworm.lg0, y0 = 1/4, 
       interval = "Wald",
       x0.name = "ldose", 
       newdata = data.frame(sex = "F"))
       
# estimate    lower    upper       se 
#   2.2313   1.7415   2.7211   0.2499
```

## Bioassay on Nasturtium

The data here contain the actual concentrations of an agrochemical present in soil samples versus the weight of the plant after three weeks of growth. These data are stored in the data frame `nasturtium` and are loaded with the package. A simple
log-logistic model describes the data well:
```r
# Log-logistic model
log_fit <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
               start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
               data = nasturtium)
plotFit(log_fit, lwd.fit = 2)
```
![Alt text](https://raw.githubusercontent.com/bgreenwell/investr/master/nasturtium_plotFit.png)

Three new replicates of the response (309, 296, 419) at an unknown concentration of interest ($x_0$) are measured. It is desired to estimate $x_0$.
```r      
# Inversion method
invest(log_fit, y0 = c(309, 296, 419), interval = "inversion")

# estimate    lower    upper 
#   2.2639   1.7722   2.9694

# Wald method
invest(log_fit, y0 = c(309, 296, 419), interval = "Wald")  

# estimate    lower    upper       se 
#   2.2639   1.6889   2.8388   0.2847
```

The intervals both rely on large sample results and normality. In practice, the bootstrap may be more reliable:
```r
# Bootstrap calibration intervals (may take a few seconds)
boo <- invest(log_fit, y0 = c(309, 296, 419), interval = "percentile", 
              nsim = 9999, seed = 101, progress = TRUE)
boo  # print bootstrap summary

# estimate    lower    upper       se     bias 
#   2.2639   1.7890   2.9380   0.2947   0.0281

plot(boo)  # plot results
```
![Alt text](https://raw.githubusercontent.com/bgreenwell/investr/master/nasturtium_boot.png)
