# The following tests are for calibration with linear regression models fit 
# using the lm() function.
context("Inverse estimation with linear models")

test_that("invest() and calibrate() produce the same results", {
  
  # Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal_lm <- lm(weight ~ time, data = crystal)
  res1.cal <- calibrate(crystal_lm, y0 = 5)
  res2.cal <- calibrate(crystal_lm, y0 = 5, mean.response = TRUE)
  res3.cal <- calibrate(crystal_lm, y0 = 5, interval = "Wald")
  res4.cal <- calibrate(crystal_lm, y0 = 5, interval = "Wald", mean.response = TRUE)
  res1.inv <-    invest(crystal_lm, y0 = 5)
  res2.inv <-    invest(crystal_lm, y0 = 5, mean.response = TRUE)
  res3.inv <-    invest(crystal_lm, y0 = 5, interval = "Wald")
  res4.inv <-    invest(crystal_lm, y0 = 5, interval = "Wald", mean.response = TRUE)
  
  # Expectations
  expect_true(all.equal(res1.cal, res1.inv, tol = 1e-05))
  expect_true(all.equal(res2.cal, res2.inv, tol = 1e-05))
  expect_true(all.equal(res3.cal, res3.inv, tol = 1e-05))
  expect_true(all.equal(res4.cal, res4.inv, tol = 1e-04))
  
})

# The following tests are for calibration with nonlinear regression models fit
# using the nls() function.
context("Inverse estimation with nonlinear models")

test_that("approximate standard error is correct", {
  
  # Nasturtium data from the drc package
  nas <- data.frame(conc = rep(c(0.000, 0.025, 0.075, 0.250, 0.750, 2.000, 
                                 4.000), each = 6),
                    weight = c(920, 889, 866, 930, 992, 1017, 919, 878, 882, 
                               854, 851, 850, 870,  825, 953, 834, 810, 875, 
                               880, 834, 795,  837, 834, 810, 693, 690, 722, 
                               738, 563,  591, 429, 395, 435, 412, 273, 257, 
                               200,  244, 209, 225, 128, 221))
  
  # Log-logistic model
  nas.nls <- nls(weight ~ theta1/(1 + exp(theta2 + theta3*log(conc))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1), 
                 data = nas)
  
  # Calculate standard errors using default and user-specified precision. The
  # estimate based on the deltaMethod function from car is 0.2847019. This was
  # calculated in the R Journal article.
  se1 <- invest(nas.nls, y0 = c(309, 296, 419), interval = "Wald")$se
  se2 <- invest(nas.nls, y0 = c(309, 296, 419), interval = "Wald", data = nas,
                tol = 1e-10)$se
  
  # Expectations
  expect_true(all.equal(se1, 0.2847019, tol = 1e-05))  # less precise
  expect_true(all.equal(se2, 0.2847019, tol = 1e-05))  # more precise
  
})

# The following tests are for linear calibration with generalized linear models
# (GzLMs) fit using the glm() function from the stats package.
context("Inverse estimation with generalized linear models")

test_that("inversion and Wald methods work", {
  
  # Dobson's beetle data
  beetle <- data.frame(
    x = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839),
    n = c(59, 60, 62, 56, 63, 59, 62, 60),
    y = c(6, 13, 18, 28, 52, 53, 61, 60)
  )
  beetle_glm <- glm(cbind(y, n-y) ~ x, data = beetle, family = "binomial")
  
  #
  # Page 207 from Categorical Data Analysis, 2nd Edition, by Alan Agresti.
  #
  
  # Inversion interval
  res <- invest(beetle_glm, y0 = 0.5, interval = "inversion", tol = 1e-10)
  a <- unname(coef(beetle_glm)[1])
  b <- unname(coef(beetle_glm)[2])
  var_a <- vcov(beetle_glm)[1, 1]
  var_b <- vcov(beetle_glm)[2, 2]
  cov_ab <- vcov(beetle_glm)[1, 2]
  fun <- function(x, p = 0.5) {
    abs(a + b*x - qlogis(p)) / sqrt(var_a + x^2*var_b + 2*x*cov_ab) - 
      qnorm(0.975)
  }
  lwr <- uniroot(fun, lower = 1.76, upper = 1.77, tol = 1e-10)$root
  upr <- uniroot(fun, lower = 1.77, upper = 1.79, tol = 1e-10)$root
  expect_true(all.equal(res$lower, lwr))
  expect_true(all.equal(res$upper, upr))
  
  #
  # Check Taylor series approximation of standard error using MASS::dose.p
  #
  
  # Wald interval
  #   mass_se <- MASS::dose.p(beetle_glm, p = 0.5)
  wald_se <- invest(beetle_glm, y0 = 0.5, interval = "Wald")$se
  expect_that(wald_se, equals(0.003858052, tol = 1e-05)) 

})

test_that("invest.glm with Gaussian family matches invest.lm", {
  
  # Using glm
  gauss_glm <- glm(weight ~ time, data = crystal, family = gaussian)
  gauss_glm_inversion <- invest(gauss_glm, y0 = 5, interval = "inversion")
  gauss_glm_wald <- invest(gauss_glm, y0 = 5, interval = "Wald")
  
  # Using lm
  gauss_lm <- glm(weight ~ time, data = crystal)
  gauss_lm_inversion <- invest(gauss_lm, y0 = 5, interval = "inversion")
  gauss_lm_wald <- invest(gauss_lm, y0 = 5, interval = "Wald")
  
  # Results should match
  expect_true(all.equal(gauss_glm_inversion$lower, gauss_lm_inversion$lower))
  expect_true(all.equal(gauss_glm_inversion$upper, gauss_lm_inversion$upper))
  expect_true(all.equal(gauss_glm_wald$upper, gauss_lm_wald$upper))
  expect_true(all.equal(gauss_glm_wald$upper, gauss_lm_wald$upper))
  expect_true(all.equal(gauss_glm_wald$se, gauss_lm_wald$se))
  
})

# The following tests are for linear calibration with random coefficient models
# fit using the lme() function from the nlme package.
context("Inverse estimation with linear mixed-effects models")

test_that("inversion method works", {

  # Load nlme package
  require(nlme)
  
  # Bladder data
  subject <- rep(1:23, times = 8)
  volume <- rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23) / 10
  HD <- c(13.2, 11.1, 10.3, NA, 4.8, 7.7, NA, 5.9, 1.9, 6.5, 19.8, 
          14.6, NA, NA, 9.7, 17.2, 10.6, 19.3, 8.5, 6.9, 8.1, 14.8, 13.7, 
          27.4, 27.5, 15, 10, 18.6, 12.6, 24, 28.4, 12.5, 16.7, 29.6, 
          27.1, 14, 18.7, 20.3, 35.8, 23.6, 37.4, 31.3, 23.7, 22, 34.3, 
          28.5, 41.6, 58.1, 34.2, 28.8, 29.9, 31.4, 46.9, 44.4, 26.8, 
          30.6, 51.7, 49.8, 19.1, 35.8, 38.9, 41.4, 49.9, 58.6, 54.8, 44, 
          39.1, 58.5, 41.5, 60.1, 78.8, 49.4, 46.4, 39.4, 45.3, 50.4, 
          70.7, 54.4, 41.8, 72.2, 67.5, 39.2, 49.6, 65.1, 69.7, 67.7, 
          73.7, 78.3, 65.7, 44.7, 72.1, 59.8, 73.9, 91.5, 71.3, 54.8, NA, 
          48, 67.8, 89.4, 63.1, 49.6, 81.9, 79.1, 48.7, 65.6, 65.1, 81.9,
          87.7, 79.4, 93, 80.3, 68.9, 90.9, 77.5, 85.5, 98.3, 81.3, 69.4, 
          NA, 66.6, 81, 105.8, 83.5, 60.8, 95.1, 95.1, 67, 85.3, 86.9, 
          96.6, 89.3, 102.6, NA, 93.6, 93.3, 105, 92.9, 95.6, 111.4, 94, 
          73.9, NA, NA, 91.2, 113.5, 114.5, 80.1, 115.4, 109.8, 72.7, 
          90.4, 98.6, 115, 108, 110.9, NA, 99.2, 102.4, 117.5, 99.4, 
          107.4, 121, 104.3, NA, NA, NA, 99.8, 127.3, 124, 87.1, NA, NA, 
          NA, NA, 107.2, 117, 114.8, 122.4, NA, 112.2, 104.7, 124.2, 113)
  bladder <- na.omit(data.frame(subject = subject, HD = HD, volume = volume))
  
  # Random intercept and slope model
  fitLME <- lme(HD^(3/2) ~ volume, random = list(subject = pdDiag(~volume)), 
                data = bladder)
  
  # Inversion method using default precision
  res <- invest(fitLME, y0 = 500, interval = "inversion")
  
  # Expectations
  expect_true(all.equal(res$estimate, 8.015521, tol = 1e-05)) 
  expect_true(all.equal(res$lower, 4.227962, tol = 1e-05))  # 4.227965
  expect_true(all.equal(res$upper, 11.91918, tol = 1e-05)) 
  
})

test_that("Wald method works", {
  
  # Load required packages
  require(nlme)
  
  # Bladder data
  subject <- rep(1:23, times = 8)
  volume <- rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23) / 10
  HD <- c(13.2, 11.1, 10.3, NA, 4.8, 7.7, NA, 5.9, 1.9, 6.5, 19.8, 
          14.6, NA, NA, 9.7, 17.2, 10.6, 19.3, 8.5, 6.9, 8.1, 14.8, 13.7, 
          27.4, 27.5, 15, 10, 18.6, 12.6, 24, 28.4, 12.5, 16.7, 29.6, 
          27.1, 14, 18.7, 20.3, 35.8, 23.6, 37.4, 31.3, 23.7, 22, 34.3, 
          28.5, 41.6, 58.1, 34.2, 28.8, 29.9, 31.4, 46.9, 44.4, 26.8, 
          30.6, 51.7, 49.8, 19.1, 35.8, 38.9, 41.4, 49.9, 58.6, 54.8, 44, 
          39.1, 58.5, 41.5, 60.1, 78.8, 49.4, 46.4, 39.4, 45.3, 50.4, 
          70.7, 54.4, 41.8, 72.2, 67.5, 39.2, 49.6, 65.1, 69.7, 67.7, 
          73.7, 78.3, 65.7, 44.7, 72.1, 59.8, 73.9, 91.5, 71.3, 54.8, NA, 
          48, 67.8, 89.4, 63.1, 49.6, 81.9, 79.1, 48.7, 65.6, 65.1, 81.9,
          87.7, 79.4, 93, 80.3, 68.9, 90.9, 77.5, 85.5, 98.3, 81.3, 69.4, 
          NA, 66.6, 81, 105.8, 83.5, 60.8, 95.1, 95.1, 67, 85.3, 86.9, 
          96.6, 89.3, 102.6, NA, 93.6, 93.3, 105, 92.9, 95.6, 111.4, 94, 
          73.9, NA, NA, 91.2, 113.5, 114.5, 80.1, 115.4, 109.8, 72.7, 
          90.4, 98.6, 115, 108, 110.9, NA, 99.2, 102.4, 117.5, 99.4, 
          107.4, 121, 104.3, NA, NA, NA, 99.8, 127.3, 124, 87.1, NA, NA, 
          NA, NA, 107.2, 117, 114.8, 122.4, NA, 112.2, 104.7, 124.2, 113)
  bladder <- na.omit(data.frame(subject = subject, HD = HD, volume = volume))

  # Random intercept and slope model
  fitLME <- lme(HD^(3/2) ~ volume, random = list(subject = pdDiag(~volume)), 
                data = bladder)
  
  # Wald method using default precision
  tvals <- qt(c(0.025, 0.975), length(resid(fitLME)) - 1)
  res.norm <- invest(fitLME, y0 = 500, interval = "Wald")
  res.t <- invest(fitLME, y0 = 500, interval = "Wald", q1 = tvals[1],
                  q2 = tvals[2])
  
  # True values calulated by hand
  ci.norm <- 8.015521 + qnorm(c(0.025, 0.975))*1.954191
  ci.t <- 8.015521 + tvals*1.954191
  
  # Expectations
  expect_true(all.equal(res.norm$estimate, 8.015521, tol = 1e-05))  # estimate
  expect_true(all.equal(res.t$estimate, 8.015521, tol = 1e-05))     # estimate
  expect_true(all.equal(res.norm$se, 1.954191, tol = 1e-05))  # SE 
  expect_true(all.equal(res.t$se, 1.954191, tol = 1e-05))     # SE 
  expect_true(all.equal(res.norm$lower, ci.norm[1], tol = 1e-05))  # lower limit
  expect_true(all.equal(res.t$lower, ci.t[1], tol = 1e-05))        # lower limit
  expect_true(all.equal(res.norm$upper, ci.norm[2], tol = 1e-05))  # upper limit
  expect_true(all.equal(res.t$upper, ci.t[2], tol = 1e-05))        # upper limit
  
})