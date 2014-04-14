context("Inverse estimation")

test_that("output matches answers from calibrate", {
  
  ## Thermostat example from Graybill and Iyer (1996, p. 431)
  read <- c(95.71, 98.16, 99.52, 102.09, 103.79, 106.18, 108.14, 110.21)
  thermom <- data.frame(read, temp = seq(from = 96, to = 110, by = 2))
  thermom.lm <- lm(read ~ temp, data = thermom)
  thermom.cal1 <- calibrate(thermom.lm, y0 = 104)
  thermom.cal2 <- calibrate(thermom.lm, y0 = 100, level = 0.9)
  thermom.cal3 <- invest(thermom.lm, y0 = 104, tol = 1e-10)
  thermom.cal4 <- invest(thermom.lm, y0 = 100, level = 0.9, tol = 1e-10)
  expect_that(thermom.cal3$estimate, equals(thermom.cal1$estimate))
  expect_that(thermom.cal3$lower, equals(thermom.cal1$lower))
  expect_that(thermom.cal3$upper, equals(thermom.cal1$upper))
  expect_that(thermom.cal4$estimate, equals(thermom.cal2$estimate))
  expect_that(thermom.cal4$lower, equals(thermom.cal2$lower))
  expect_that(thermom.cal4$upper, equals(thermom.cal2$upper))
  
  ## Reaction chamber example from Graybill and Iyer (1996, p. 433)
  dial <- seq(from = 0, to = 100, by = 10)
  temp <- c(206.36, 225.52, 252.18, 289.33, 318.11, 349.49, 383.03, 410.70, 
            444.40, 469.14, 501.16)
  chamber <- data.frame(dial, temp)
  chamber.reg <- calibrate(chamber.lm, y0 = 400, mean.response = TRUE, 
                           level = 0.99)
  chamber.reg2 <- invest(chamber.lm, y0 = 400, mean.response = TRUE, 
                         level = 0.99, tol = 1e-10)
  expect_that(round(chamber.reg2$estimate, 1), equals(66.5))
  expect_that(round(chamber.reg2$lower, 2), equals(65.07))
  expect_that(round(chamber.reg2$upper, 2), equals(68.03))
  expect_that(chamber.reg2$estimate, equals(chamber.reg$estimate))
  expect_that(chamber.reg2$lower, equals(chamber.reg$lower))
  expect_that(chamber.reg2$upper, equals(chamber.reg$upper))
  
  ## Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.reg <- calibrate(crystal.lm, y0 = 5, mean.response = TRUE, 
                           level = 0.9)
  expect_that(round(crystal.reg$estimate, 2), equals(9.93))
  expect_that(round(crystal.reg$lower, 2), equals(8.65))
  expect_that(round(crystal.reg$upper, 2), equals(11.05))
  
})

test_that("standard error matches the one from calibrate", {
  
  ## Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal.lm <- lm(weight ~ time, data = crystal)
  crystal.cal1 <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
  crystal.reg1 <- calibrate(crystal.lm, y0 = 5, interval = "Wald", 
                            mean.response = TRUE)
  crystal.cal2 <- invest(crystal.lm, y0 = 5, interval = "Wald", tol = 1e-10)
  crystal.reg2 <- invest(crystal.lm, y0 = 5, interval = "Wald", 
                         mean.response = TRUE, tol = 1e-10)
  expect_that(round(crystal.cal2$se, 4), equals(round(se.cal, 4))) # small diff
  expect_that(round(crystal.reg2$se, 4), equals(round(se.reg, 4)))
  
})

test_that("standard error matches the one from car::deltaMethod", {
  
  ## Nasturtium data from drc package
  weight <- c(920, 889, 866, 930, 992, 1017, 919, 878, 882, 854, 851, 850, 870,
              825, 953, 834, 810, 875,  880, 834, 795, 837, 834, 810, 693, 690,
              722, 738, 563, 591, 429,  395, 435, 412, 273, 257, 200, 244, 209,
              225, 128, 221)
  conc <- rep(c(0, 0.025, 0.075, 0.250, 0.750, 2.000, 4.000), each = 6)
  nas <- data.frame(conc, weight)
  nas.nls <- nls(weight ~ ifelse(conc == 0, theta1, 
                                 theta1/(1 + exp(theta2 + theta3*log(conc)))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
                 data = nas)
  nas.cal <- invest(nas.nls, y0 = c(309, 296, 419), interval = "Wald", 
                    tol = 1e-10)
  
  ## Calculate and compare standard error using packages investr and car
#   covmat <- diag(4)
#   covmat[1:3, 1:3] <- vcov(nas.nls)
#   covmat[4, 4] <- summary(nas.nls)$sigma^2/3
#   dm <- deltaMethod(c(coef(nas.nls), eta = mean(c(309, 296, 419))), 
#               g = "exp((log(theta1/eta - 1) - theta2) / theta3)", 
#               vcov. = covmat)
  expect_that(round(nas.cal$estimate, 6), equals(round(2.263852, 6)))
  expect_that(round(nas.cal$se, 7), equals(round(0.2847019, 7)))
  
})

test_that("errors are handled appropriately", {
  
  ## Nasturtium data from drc package
  weight <- c(920, 889, 866, 930, 992, 1017, 919, 878, 882, 854, 851, 850, 870,
              825, 953, 834, 810, 875,  880, 834, 795, 837, 834, 810, 693, 690,
              722, 738, 563, 591, 429,  395, 435, 412, 273, 257, 200, 244, 209,
              225, 128, 221)
  conc <- rep(c(0, 0.025, 0.075, 0.250, 0.750, 2.000, 4.000), each = 6)
  nas <- data.frame(conc, weight)
  nas.nls <- nls(weight ~ ifelse(conc == 0, theta1, 
                                 theta1/(1 + exp(theta2 + theta3*log(conc)))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
                 data = nas)
  expect_that(invest(nas.nls, y0 = 175), throws_error())
  expect_that(invest(nas.nls, y0 = 175, upper = 10), throws_error())

})
