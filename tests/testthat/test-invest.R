context("Inverse estimation")

test_that("correct output for crystal data to four decimal places", {
  
  ## Setup
  crystal.lm <- lm(weight ~ time, data = crystal)
  cal.inv <- invest(crystal.lm, y0 = 8, interval = "inversion")
  cal.del <- invest(crystal.lm, y0 = 8, interval = "Wald")
  reg.inv <- invest(crystal.lm, y0 = 8, interval = "inversion", 
                       mean.response = TRUE)
  reg.del <- invest(crystal.lm, y0 = 8, interval = "Wald", 
                       mean.response = TRUE)
  cal.inv.est <- round(cal.inv$estimate, digits = 3)
  cal.del.est <- round(cal.del$estimate, digits = 3)
  cal.inv.lwr <- round(cal.inv$lower, digits = 3)
  cal.del.lwr <- round(cal.del$lower, digits = 3)
  cal.inv.upr <- round(cal.inv$upper, digits = 3)
  cal.del.upr <- round(cal.del$upper, digits = 3)
  cal.del.se <- round(cal.del$se, digits = 3)
  reg.inv.est <- round(reg.inv$estimate, digits = 3)
  reg.del.est <- round(reg.del$estimate, digits = 3)
  reg.inv.lwr <- round(reg.inv$lower, digits = 3)
  reg.del.lwr <- round(reg.del$lower, digits = 3)
  reg.inv.upr <- round(reg.inv$upper, digits = 3)
  reg.del.upr <- round(reg.del$upper, digits = 3)
  reg.del.se <- round(reg.del$se, digits = 3)
  
  ## Expectations
  expect_that(cal.inv.est, equals(15.888))
  expect_that(cal.del.est, equals(15.888))
  expect_that(cal.inv.lwr, equals(11.095))
  expect_that(cal.del.lwr, equals(11.130))
  expect_that(cal.inv.upr, equals(20.724))
  expect_that(cal.del.upr, equals(20.647))
  expect_that(cal.del.se, equals(2.184))
  expect_that(reg.inv.est, equals(15.888))
  expect_that(reg.del.est, equals(15.888))
  expect_that(reg.inv.lwr, equals(14.659))
  expect_that(reg.del.lwr, equals(14.653 ))
  expect_that(reg.inv.upr, equals(17.160))
  expect_that(reg.del.upr, equals(17.124))
  expect_that(reg.del.se, equals(0.567))
  
})

test_that("correct output for arsenic data to four decimal places", {
  
  ## Setup
  arsenic.lm <- lm(measured ~ actual, data = arsenic)
  cal.inv <- invest(arsenic.lm, y0 = 3, interval = "inversion")
  cal.del <- invest(arsenic.lm, y0 = 3, interval = "Wald")
  reg.inv <- invest(arsenic.lm, y0 = 3, interval = "inversion", 
                       mean.response = TRUE)
  reg.del <- invest(arsenic.lm, y0 = 3, interval = "Wald", 
                       mean.response = TRUE)
  cal.inv.est <- round(cal.inv$estimate, digits = 3)
  cal.del.est <- round(cal.del$estimate, digits = 3)
  cal.inv.lwr <- round(cal.inv$lower, digits = 3)
  cal.del.lwr <- round(cal.del$lower, digits = 3)
  cal.inv.upr <- round(cal.inv$upper, digits = 3)
  cal.del.upr <- round(cal.del$upper, digits = 3)
  cal.del.se <- round(cal.del$se, digits = 3)
  reg.inv.est <- round(reg.inv$estimate, digits = 3)
  reg.del.est <- round(reg.del$estimate, digits = 3)
  reg.inv.lwr <- round(reg.inv$lower, digits = 3)
  reg.del.lwr <- round(reg.del$lower, digits = 3)
  reg.inv.upr <- round(reg.inv$upper, digits = 3)
  reg.del.upr <- round(reg.del$upper, digits = 3)
  reg.del.se <- round(reg.del$se, digits = 3)
  
  ## Expectations
  expect_that(cal.inv.est, equals(2.931))
  expect_that(cal.del.est, equals(2.931))
  expect_that(cal.inv.lwr, equals(2.537))
  expect_that(cal.del.lwr, equals(2.537))
  expect_that(cal.inv.upr, equals(3.325))
  expect_that(cal.del.upr, equals(3.325))
  expect_that(cal.del.se, equals(0.193))
  expect_that(reg.inv.est, equals(2.931))
  expect_that(reg.del.est, equals(2.931))
  expect_that(reg.inv.lwr, equals(2.860))
  expect_that(reg.del.lwr, equals(2.861))
  expect_that(reg.inv.upr, equals(3.002))
  expect_that(reg.del.upr, equals(3.002))
  expect_that(reg.del.se, equals(0.035))
  
})
