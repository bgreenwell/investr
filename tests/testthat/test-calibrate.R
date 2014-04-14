context("Simple linear calibration")

test_that("output matches answers to Graybill and Iyer (1996, chap. 6)", {
    
  ## Thermostat example from Graybill and Iyer (1996, p. 431)
  thermom <- data.frame(temp = seq(from = 96, to = 110, by = 2), 
                        read = c(95.71, 98.16, 99.52, 102.09, 103.79, 106.18, 
                                 108.14, 110.21))
  thermom.cal1 <- calibrate(thermom, y0 = 104)
  thermom.cal2 <- calibrate(thermom, y0 = 100, level = 0.9)
  expect_that(round(thermom.cal1$estimate, 3), equals(103.995))
  expect_that(round(thermom.cal1$lower, 1), equals(103.4))
  expect_that(round(thermom.cal1$upper, 1), equals(104.6))
  expect_that(round(thermom.cal2$estimate, 1), equals(100.1))
  expect_that(round(thermom.cal2$lower, 2), equals(99.63))
  expect_that(round(thermom.cal2$upper, 2), equals(100.59))
  
  ## Reaction chamber example from Graybill and Iyer (1996, p. 433)
  chamber <- data.frame(dial = seq(from = 0, to = 100, by = 10), 
                        temp = c(206.36, 225.52, 252.18, 289.33, 318.11, 349.49, 
                                 383.03, 410.70, 444.40, 469.14, 501.16))
  chamber.reg <- calibrate(chamber, y0 = 400, mean.response = TRUE, 
                           level = 0.99)
  expect_that(round(chamber.reg$estimate, 1), equals(66.5))
  expect_that(round(chamber.reg$lower, 2), equals(65.07))
  expect_that(round(chamber.reg$upper, 2), equals(68.03))
  
  ## Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.reg <- calibrate(crystal.lm, y0 = 5, mean.response = TRUE, 
                           level = 0.9)
  expect_that(round(crystal.reg$estimate, 2), equals(9.93))
  expect_that(round(crystal.reg$lower, 2), equals(8.65))
  expect_that(round(crystal.reg$upper, 2), equals(11.05))
  
})

test_that("standard error matches the one from car::deltaMethod", {
  
  ## Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
  crystal.reg <- calibrate(crystal.lm, y0 = 5, interval = "Wald", 
                           mean.response = TRUE)
  
  ## Calculate and compare standard error using invest and car::deltaMethod
  covmat.cal <- diag(3)
  covmat.cal[1:2, 1:2] <- vcov(crystal.lm)
  covmat.cal[3, 3] <- summary(crystal.lm)$sigma^2
  coefs <- unname(coef(crystal.lm))
  params <- c(b0 = coefs[1], b1 = coefs[2], y0 = 5)
  se.cal <- 2.211698 #car::deltaMethod(params, g = "(y0-b0)/b1", vcov. = covmat.cal)$SE
  se.reg <- 0.6658998 #car::deltaMethod(crystal.lm, g = "(5-b0)/b1", 
                             #parameterNames = c("b0", "b1"))$SE
  expect_that(round(crystal.cal$se, 5), equals(round(se.cal, 5))) # small diff
  expect_that(crystal.reg$se, equals(se.reg))
  
})

test_that("errors are handled appropriately", {
  
  ## Simulated data
  set.seed(101)
  x <- rep(seq(from = 0, to = 10, length = 10), 2)
  y <-  3 + 0.01*x + rnorm(length(x), sd = 0.5)
  d1 <- data.frame(x, y)
  d2 <- list(x = 1:11, y = 1:10 + rnorm(10, sd = 1))
#   fit <- lm(y ~ x, data = d)
#   plotFit(fit, interval = "both", xlim = c(-10, 25))
#   abline(h = c(2, 3), col = "red")
  expect_that(calibrate(d1, y0 = 3), gives_warning())
  expect_that(calibrate(d1, y0 = 3, mean.response = TRUE), gives_warning())
  expect_that(calibrate(d1, y0 = 2), throws_error())
  expect_that(calibrate(d1, y0 = 2.5, mean.response = TRUE), throws_error())
  expect_that(calibrate(d2, y0 = 2.5), throws_error())
  expect_that(calibrate(y ~ x + I(x^2), y0 = 2.5), throws_error())

})
