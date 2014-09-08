## Tests for invest.R

## TODO: Make sure invest produces the same results as calibrate
context("Numerical linear calibration/regulation")

test_that("invest and calibrate produce the same results", {
  
  ## Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  res1.cal <- calibrate(crystal.lm, y0 = 5)
  res2.cal <- calibrate(crystal.lm, y0 = 5, mean.response = TRUE)
  res3.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
  res4.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald", mean.response = TRUE)
  res1.inv <-    invest(crystal.lm, y0 = 5)
  res2.inv <-    invest(crystal.lm, y0 = 5, mean.response = TRUE)
  res3.inv <-    invest(crystal.lm, y0 = 5, interval = "Wald")
  res4.inv <-    invest(crystal.lm, y0 = 5, interval = "Wald", mean.response = TRUE)
  
  ## Expectations
  expect_true(all.equal(res1.cal, res1.inv, tol = 1e-05))
  expect_true(all.equal(res2.cal, res2.inv, tol = 1e-05))
  expect_true(all.equal(res3.cal, res3.inv, tol = 1e-05))
  expect_true(all.equal(res4.cal, res4.inv, tol = 1e-04))
  
})

## TODO: Test inversion/Wald interval on cars data frame?
context("Polynomial calibration/regulation - inversion interval")
context("Polynomial calibration/regulation - Wald interval")

## TODO: Test inversion/Wald interval on nasturtium data frame?
context("Nonlinear calibration/regulation - inversion interval")

context("Nonlinear calibration/regulation - Wald interval")

test_that("approximate standard error is correct", {
  
  ## Nasturtium data from the drc package
  nas <- data.frame(conc = rep(c(0.000, 0.025, 0.075, 0.250, 0.750, 2.000, 
                                 4.000), each = 6),
                    weight = c(920, 889, 866, 930, 992, 1017, 919, 878, 882, 
                               854, 851, 850, 870,  825, 953, 834, 810, 875, 
                               880, 834, 795,  837, 834, 810, 693, 690, 722, 
                               738, 563,  591, 429, 395, 435, 412, 273, 257, 
                               200,  244, 209, 225, 128, 221))
  nas.nls <- nls(weight ~ theta1/(1 + exp(theta2 + theta3*log(conc))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1), 
                 data = nas)
  
#   getData(nas.nls)
  
#   se1 <- invest(nas.nls, y0 = c(309, 296, 419), interval = "Wald")$se
  se2 <- invest(nas.nls, y0 = c(309, 296, 419), interval = "Wald", data = nas,
                tol = 1e-10)$se
#   expect_false(all.equal(se1, 0.2847019, tol = 1e-07))
  expect_true(all.equal(se2, 0.2847019, tol = 1e-07))
  
})