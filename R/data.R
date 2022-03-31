#' Crystal weight data
#'
#' The data give the growing time and final weight of crystals.
#'
#' \itemize{
#'   \item \code{time} Time taken to grow (hours).
#'   \item \code{weight} Final weight of the crystal (grams).
#' }
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @format A data frame with 14 rows and 2 columns.
#' 
#' @name crystal
#' 
#' @source
#' Graybill, F. A., and Iyer, H. K. (1994)
#' \emph{Regression analysis: Concepts and Applications}. Duxbury Press.
NULL


#' Concentrations of arsenic in water samples
#'
#' The data give the actual and measured concentrations of arsenic present in
#' water samples.
#'
#' \itemize{
#'   \item \code{actual} True amount of arsenic present.
#'   \item \code{measured} Measured amount of arsenic present.
#' }
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @format A data frame with 32 rows and 2 columns.
#' 
#' @name arsenic
#' 
#' @source
#' Graybill, F. A., and Iyer, H. K. (1994)
#' \emph{Regression analysis: Concepts and Applications}. Duxbury Press.
NULL


#' Bioassay on Nasturtium
#'
#' The data give the actual concentrations of an agrochemical present in soil
#' samples versus the weight of the plant after three weeks of growth.
#'
#' \itemize{
#'   \item \code{conc} True concentration of agrochemical (g/ha).
#'   \item \code{weight} Weight of plant (mg) after 3 weeks' growth.
#' }
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @format A data frame with 42 rows and 2 columns.
#' 
#' @name nasturtium
#' 
#' @source
#' Racine-Poon, A. (1988) A Bayesian Approach to Nonlinear Calibration
#' Problems, \emph{Journal of the American Statistical Association}, \bold{83},
#' 650--656.
#' 
#' @references
#' Huet, S., Bouvier, A., Poursat, M-A., and Jolivet, E. (2004)
#' \emph{Statistical Tools for Nonlinear Regression: A Practical Guide with
#' S-PLUS and R Examples}. Springer.
#'
NULL


#' Dobson's Beetle Data
#'
#' The data give the number of flour beetles killed after five hour exposure
#' to the insecticide carbon disulphide at eight different concentrations.
#'
#' \itemize{
#'   \item \code{ldose} Log dose of carbon disulphide.
#'   \item \code{y} Number of beetles subjected to insecticide.
#'   \item \code{n} Number of beetles killed.
#' }
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @format A data frame with 8 rows and 3 columns.
#' 
#' @name beetle
#' 
#' @source
#' Dobson, A. (2002) \emph{An Introduction to Generalized Linear Models}. 
#' Chapman & Hall/CRC.
#'
NULL


#' Bladder volume data
#'
#' A series of 23 women patients attending a urodynamic clinic were recruited 
#' for a study. After successful voiding of the bladder, sterile water was 
#' introduced in additions of 10, 15, and then 25 ml increments up to a final 
#' cumulative total of 175 ml. At each volume a measure of height (H) in mm and 
#' depth (D) in mm of largest ultrasound bladder images were taken. The product 
#' H × D was taken as a measure of liquid volume.
#'
#' \itemize{
#'   \item \code{subject} The subject ID.
#'   \item \code{HD} The product H × D (mm^2).
#'   \item \code{volume} The true volume of sterile water in the bladder (ml).
#' }
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @format A data frame with 184 rows and 3 columns.
#' 
#' @name bladder
#' 
#' @source
#' Brown, P. (1993) \emph{Measurement, Regression, and Calibration}. 
#' Oxford University Press.
NULL


#' Whisky data
#'
#' The data give the proof (measured as twice the percentage of alcohol by 
#' volume, denoted 2ABV) of whiskey stored in a charred oak barrel against time 
#' in years.
#'
#' \itemize{
#'   \item \code{age} The age of the whisky (years).
#'   \item \code{proof} The proof of the whisky (2ABV).
#' }
#' 
#' @docType data
#' 
#' @keywords datasets
#' 
#' @format A data frame with 10 rows and 2 columns.
#' 
#' @name whisky
#' 
#' @source
#' R. Schoeneman, R. Dyer, and E. Earl. Analytical profile of straight bourbon 
#' whiskies. Journal of the Association of Official Analytical Chemists, 
#' 54:1247--1261, 1971.
NULL