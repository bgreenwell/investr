#' investr: a package for inverse estimation in R
#' 
#' Inverse estimation, also referred to as the calibration problem, is a 
#' classical and well-known problem in regression. In simple terms, it involves 
#' the use of an observed value of the response (or specified value of the mean 
#' response) to make inference on the corresponding unknown value of the 
#' explanatory variable. 
#'
#' A detailed [introduction to investr](https://journal.r-project.org/archive/2014/RJ-2014-009/index.html) has been published in The R Journal: 
#' "investr: An R Package for Inverse Estimation." You 
#' can track development at <https://github.com/bgreenwell/investr>. To report 
#' bugs or issues, contact the main author directly or submit them to 
#' <https://github.com/bgreenwell/investr/issues>. 
#'
#' As of right now, `investr` supports (univariate) inverse estimation 
#' with objects of class:
#' \itemize{
#'   \item{`lm`} --- linear models (multiple predictor variables allowed)
#'   \item{`glm`} --- generalized linear models (multiple predictor variables allowed)
#'   \item{`nls`} --- nonlinear least-squares models
#'   \item{`lme`} --- linear mixed-effects models (fit using the 
#'     `nlme` package)
#' }
#' 
#' @keywords internal
"_PACKAGE"
