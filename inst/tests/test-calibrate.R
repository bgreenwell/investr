context("Calibration")

test_that("correct output for crystal data to four decimal places", {
    
  ## Setup
  crystal.lm <- lm(weight ~ time, data = crystal)
  cal.inv <- calibrate(crystal.lm, y0 = 8, interval = "inversion")
  cal.del <- calibrate(crystal.lm, y0 = 8, interval = "Wald")
  reg.inv <- calibrate(crystal.lm, y0 = 8, interval = "inversion", 
                       mean.response = TRUE)
  reg.del <- calibrate(crystal.lm, y0 = 8, interval = "Wald", 
                       mean.response = TRUE)
  cal.inv.est <- round(cal.inv$estimate, digits = 4)
  cal.del.est <- round(cal.del$estimate, digits = 4)
  cal.inv.lwr <- round(cal.inv$lower, digits = 4)
  cal.del.lwr <- round(cal.del$lower, digits = 4)
  cal.inv.upr <- round(cal.inv$upper, digits = 4)
  cal.del.upr <- round(cal.del$upper, digits = 4)
  cal.del.se <- round(cal.del$se, digits = 4)
  reg.inv.est <- round(reg.inv$estimate, digits = 4)
  reg.del.est <- round(reg.del$estimate, digits = 4)
  reg.inv.lwr <- round(reg.inv$lower, digits = 4)
  reg.del.lwr <- round(reg.del$lower, digits = 4)
  reg.inv.upr <- round(reg.inv$upper, digits = 4)
  reg.del.upr <- round(reg.del$upper, digits = 4)
  reg.del.se <- round(reg.del$se, digits = 4)
  
  ## Expectations
  expect_that(cal.inv.est, equals(15.8882))
  expect_that(cal.del.est, equals(15.8882))
  expect_that(cal.inv.lwr, equals(11.0946))
  expect_that(cal.del.lwr, equals(11.1297))
  expect_that(cal.inv.upr, equals(20.7240))
  expect_that(cal.del.upr, equals(20.6467))
  expect_that(cal.del.se, equals(2.1840))
  expect_that(reg.inv.est, equals(15.8882))
  expect_that(reg.del.est, equals(15.8882))
  expect_that(reg.inv.lwr, equals(14.6590))
  expect_that(reg.del.lwr, equals(14.6526 ))
  expect_that(reg.inv.upr, equals(17.1596))
  expect_that(reg.del.upr, equals(17.1238))
  expect_that(reg.del.se, equals(0.5671))
  
})

test_that("correct output for arsenic data to four decimal places", {
    
  ## Setup
  arsenic.lm <- lm(measured ~ actual, data = arsenic)
  cal.inv <- calibrate(arsenic.lm, y0 = 3, interval = "inversion")
  cal.del <- calibrate(arsenic.lm, y0 = 3, interval = "Wald")
  reg.inv <- calibrate(arsenic.lm, y0 = 3, interval = "inversion", 
                       mean.response = TRUE)
  reg.del <- calibrate(arsenic.lm, y0 = 3, interval = "Wald", 
                       mean.response = TRUE)
  cal.inv.est <- round(cal.inv$estimate, digits = 4)
  cal.del.est <- round(cal.del$estimate, digits = 4)
  cal.inv.lwr <- round(cal.inv$lower, digits = 4)
  cal.del.lwr <- round(cal.del$lower, digits = 4)
  cal.inv.upr <- round(cal.inv$upper, digits = 4)
  cal.del.upr <- round(cal.del$upper, digits = 4)
  cal.del.se <- round(cal.del$se, digits = 4)
  reg.inv.est <- round(reg.inv$estimate, digits = 4)
  reg.del.est <- round(reg.del$estimate, digits = 4)
  reg.inv.lwr <- round(reg.inv$lower, digits = 4)
  reg.del.lwr <- round(reg.del$lower, digits = 4)
  reg.inv.upr <- round(reg.inv$upper, digits = 4)
  reg.del.upr <- round(reg.del$upper, digits = 4)
  reg.del.se <- round(reg.del$se, digits = 4)
  
  ## Expectations
  expect_that(cal.inv.est, equals(2.9314))
  expect_that(cal.del.est, equals(2.9314))
  expect_that(cal.inv.lwr, equals(2.5367))
  expect_that(cal.del.lwr, equals(2.5374))
  expect_that(cal.inv.upr, equals(3.3251))
  expect_that(cal.del.upr, equals(3.3255))
  expect_that(cal.del.se, equals(0.1929))
  expect_that(reg.inv.est, equals(2.9314))
  expect_that(reg.del.est, equals(2.9314))
  expect_that(reg.inv.lwr, equals(2.8603))
  expect_that(reg.del.lwr, equals(2.8608))
  expect_that(reg.inv.upr, equals(3.0016))
  expect_that(reg.del.upr, equals(3.0021))
  expect_that(reg.del.se, equals(0.0346))
  
})

