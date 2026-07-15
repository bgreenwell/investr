# Inverse estimation in the SLR model

# output matches answers to Graybill and Iyer (1996, chap. 6) ----------------

# Thermostat example from Graybill and Iyer (1996, p. 431)
thermom <- data.frame(temp = seq(from = 96, to = 110, by = 2),
                      read = c(95.71, 98.16, 99.52, 102.09, 103.79, 106.18,
                               108.14, 110.21))
thermom.cal1 <- calibrate(thermom, y0 = 104)
thermom.cal2 <- calibrate(thermom, y0 = 100, level = 0.9)
expect_equal(round(thermom.cal1$estimate, 3), 103.995)
expect_equal(round(thermom.cal1$lower, 1), 103.4)
expect_equal(round(thermom.cal1$upper, 1), 104.6)
expect_equal(round(thermom.cal2$estimate, 1), 100.1)
expect_equal(round(thermom.cal2$lower, 2), 99.63)
expect_equal(round(thermom.cal2$upper, 2), 100.59)

# Reaction chamber example from Graybill and Iyer (1996, p. 433)
chamber <- data.frame(dial = seq(from = 0, to = 100, by = 10),
                      temp = c(206.36, 225.52, 252.18, 289.33, 318.11, 349.49,
                               383.03, 410.70, 444.40, 469.14, 501.16))
chamber.reg <- calibrate(chamber, y0 = 400, mean.response = TRUE,
                         level = 0.99)
expect_equal(round(chamber.reg$estimate, 1), 66.5)
expect_equal(round(chamber.reg$lower, 2), 65.07)
expect_equal(round(chamber.reg$upper, 2), 68.03)

# Crystal weight example from Graybill and Iyer (1996, p. 434)
crystal.lm <- lm(weight ~ time, data = crystal)
crystal.reg <- calibrate(crystal.lm, y0 = 5, mean.response = TRUE,
                         level = 0.9)
expect_equal(round(crystal.reg$estimate, 2), 9.93)
expect_equal(round(crystal.reg$lower, 2), 8.65)
expect_equal(round(crystal.reg$upper, 2), 11.05)

# errors are handled appropriately --------------------------------------------

set.seed(101)
x <- rep(seq(from = 0, to = 10, length = 10), 2)
y <-  3 + 0.01*x + rnorm(length(x), sd = 0.5)
d1 <- data.frame(x, y)
d2 <- list(x = 1:11, y = 1:10 + rnorm(10, sd = 1))

expect_warning(calibrate(d1, y0 = 3))
expect_warning(calibrate(d1, y0 = 3, mean.response = TRUE))
expect_error(calibrate(d1, y0 = 2))
expect_error(calibrate(d1, y0 = 2.5, mean.response = TRUE))
expect_error(calibrate(d2, y0 = 2.5))
expect_error(calibrate(y ~ x + I(x^2), y0 = 2.5))

# approximate standard error is correct ---------------------------------------

crystal.lm <- lm(weight ~ time, data = crystal)
crystal.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
crystal.reg <- calibrate(crystal.lm, y0 = 5, interval = "Wald",
                         mean.response = TRUE)

# Calculated using car::deltaMethod (see the R Journal article)
se.cal <- 2.211698
se.reg <- 0.6658998

expect_equal(crystal.cal$se, se.cal, tol = 1e-04)  # small diff
expect_equal(crystal.reg$se, se.reg, tol = 1e-04)

# all methods produce equivalent results --------------------------------------

set.seed(101)
x <- rep(1:10, each = 3)
y <- 2 + 3 * x + rnorm(length(x), sd = 1)
d <- data.frame(x = x, y = y)

# Inversion interval, across all supported input methods
cal1 <- calibrate(cbind(x, y), y0 = 15, mean.response = FALSE)
cal2 <- calibrate(data.frame(x, y), y0 = 15, mean.response = FALSE)
cal3 <- calibrate(list(x, y), y0 = 15, mean.response = FALSE)
cal4 <- calibrate(y ~ x, y0 = 15, mean.response = FALSE)
cal5 <- calibrate(y ~ x, data = d, y0 = 15, mean.response = FALSE)
cal6 <- calibrate(lm(y ~ x), y0 = 15, mean.response = FALSE)
cal7 <- calibrate(exp(log(y)) ~ sqrt(x^2), y0 = 15, mean.response = FALSE)

# Wald interval, across all supported input methods
cal8 <- calibrate(cbind(x, y), y0 = 15, mean.response = FALSE,
                  interval = "Wald")
cal9 <- calibrate(data.frame(x, y), y0 = 15, mean.response = FALSE,
                  interval = "Wald")
cal10 <- calibrate(list(x, y), y0 = 15, mean.response = FALSE,
                   interval = "Wald")
cal11 <- calibrate(y ~ x, y0 = 15, mean.response = FALSE, interval = "Wald")
cal12 <- calibrate(y ~ x, data = d, y0 = 15, mean.response = FALSE,
                  interval = "Wald")
cal13 <- calibrate(lm(y ~ x), y0 = 15, mean.response = FALSE,
                   interval = "Wald")
cal14 <- calibrate(exp(log(y)) ~ sqrt(x^2), y0 = 15, mean.response = FALSE,
                   interval = "Wald")

expect_identical(cal1, cal2)
expect_identical(cal1, cal3)
expect_identical(cal1, cal4)
expect_identical(cal1, cal5)
expect_identical(cal1, cal6)
expect_equal(cal1, cal7)  # not identical due to the formula transformations

expect_identical(cal8, cal9)
expect_identical(cal8, cal10)
expect_identical(cal8, cal11)
expect_identical(cal8, cal12)
expect_identical(cal8, cal13)
expect_equal(cal8, cal14)  # not identical due to the formula transformations

# errors get handled appropriately, including poly()-transformed predictors --
# (regression test for issue #48 / PR #49: calibrate.lm() must read x/y from
# the model frame, not model.matrix(), so poly() terms aren't expanded)

nls.fit <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
               start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
               data = nasturtium)
mlr.fit1 <- lm(weight ~ time + I(time ^ 2), data = crystal)
mlr.fit2 <- lm(cbind(weight, weight ^ 2) ~ time, data = crystal)
mlr.poly <- lm(weight ~ poly(time, 2), data = crystal)

expect_error(calibrate(nls.fit, y0 = c(309, 296, 419)))
expect_error(calibrate(mlr.fit1, y0 = c(309, 296, 419)))
expect_error(calibrate(mlr.fit2, y0 = c(309, 296, 419)))
expect_error(calibrate(mlr.poly, y0 = c(309, 296, 419)))

# multiple inference procedures work ------------------------------------------

crystal.lm <- lm(weight ~ time, data = crystal)
crystal.cal <- calibrate(crystal.lm, y0 = 5, interval = "Wald")
crystal.cal.multi <- calibrate(crystal.lm, y0 = 5, interval = "Wald",
                               adjust = "Scheffe", k = 1)

expect_equal(crystal.cal, crystal.cal.multi)
