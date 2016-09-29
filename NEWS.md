# NEWS for investr package

### Changes for version 1.4.1

* Better y-axis limits when x-axis is on log scale (e.g., using log = "x"). Fixes issue #27.
* Better default y-axis label when plotting a model with transformed response; for example, `log(y) ~ x`.

### Changes for version 1.4.0

* Multiple predictor variables are allowed for "lm" and "glm" objects.
* All non-base package functions are now imported.
* The generic function predFit is now exported. This function is used by investr 
  to obtain predictions, and hence, inverse predictions. For example, predFit 
  can be used to obtain prediction intervals for nonlinear least-squares fits 
  (i.e., models of class "nls").
* Improved tests and test coverage.
* plotFit gained methods for "rlm" and "lqs" objects from package MASS.

### Changes for version 1.3.0

* invest now accepts objects of class 'glm' (experimental).
* Functions calibrate and invest now return an object of class "invest".

### Changes for version 1.2.1

* Cleaned up documentation.
* Added AnyNA function for those using older versions of R.

### Changes for version 1.2.0

* Cleaned up examples.
* Added bootstrap option to invest.

### Changes for version 1.1.2

* Changed tests to satisfy CRAN check.

### Changes for version 1.1.1

* Updated citation file.
* Minor code changes.
* plotFit should now plot models with transformed responses correctly.
* Fixed error causing invest to fail because of a missing data argument.
* Added more tests.

### Changes for version 1.1.0

* invest now accepts objects of class 'lme' (experimental).
* A few minor bug fixes and code improvements.
* Added more tests.

### Changes for version 1.0.1

* A few minor bug fixes.
* Slightly better documentation.

### Changes for version 1.0

* Added functions.
* Fixed roxygen documentation.
