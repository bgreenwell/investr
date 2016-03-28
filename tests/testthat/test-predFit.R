## The following tests are for testing the prediction function and estimated
## standard errors of the fitted values
context("Prediction")

test_that("results from predFit match results from PROC NLIN in SAS", {
  
  ## DNase data from the dataframes package
  DNase1 <- data.frame(conc = c(0.04882812, 0.04882812, 0.19531250, 0.19531250, 
                                0.39062500, 0.39062500, 0.78125000, 0.78125000, 
                                1.56250000, 1.56250000, 3.12500000, 3.12500000, 
                                6.25000000, 6.25000000, 12.50000000, 12.50000000),
                       density = c(0.017, 0.018, 0.121, 0.124, 0.206, 0.215, 
                                   0.377, 0.374, 0.614, 0.609, 1.019, 1.001, 
                                   1.334, 1.364, 1.730, 1.710))
  DNase1.new <- data.frame(conc = c(8.0, 12.0, 1.0, 10.0, 5.5))
  
  DNase1.nls <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
                    data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
  DNase1.conf <- predFit(DNase1.nls, interval = "confidence")
  DNase1.conf2 <- predFit(DNase1.nls, newdata = DNase1.new, 
                           interval = "confidence")
  DNase1.pred <- predFit(DNase1.nls, interval = "prediction")
  DNase1.pred2 <- predFit(DNase1.nls, newdata = DNase1.new, 
                           interval = "prediction")
  
  ## Fitted value standard errors from PROC NLIN in SAS/STATS
  SAS.se.fit <- c(0.002899, 0.002899, 0.006079, 0.006079, 0.007461, 0.007461, 
                  0.007704, 0.007704, 0.007581, 0.007581, 0.009277, 0.009277,
                  0.009038, 0.009038, 0.012925, 0.012925)
  SAS.se.fit.new <- c(0.008857, 0.012272, 0.007548, 0.010014, 0.009314)
  
  ## Confidence limits from PROC NLIN in SAS/STATS
  SAS.conf.lwr <- c(0.02442, 0.02442, 0.09892, 0.09892, 0.19246, 0.19246,
                    0.35768, 0.35768, 0.61640, 0.61640, 0.96082, 0.96082,
                    1.34799, 1.34799, 1.68706, 1.68706)
  SAS.conf.upr <- c(0.03694, 0.03694, 0.12518, 0.12518, 0.22470, 0.22470,
                    0.39097, 0.39097, 0.64916, 0.64916, 1.00090, 1.00090,
                    1.38704, 1.38704, 1.74291, 1.74291)
  SAS.conf.lwr.new <- c(1.48029, 1.67025, 0.43872, 1.58988, 1.27678)
  SAS.conf.upr.new <- c(1.51856, 1.72327, 0.47133, 1.63315, 1.31703)
  
  ## Prediction limits from PROC NLIN in SAS/STATS
  SAS.pred.lwr <- c(-0.01126, -0.01126,  0.06855,  0.06855,  0.16409,  0.16409,  
                    0.32965,  0.32965,  0.58819,  0.58819,  0.93481,  0.93481,
                    1.32168,  1.32168,  1.66500,  1.66500)
  SAS.pred.upr <- c(0.07262, 0.07262, 0.15555, 0.15555, 0.25307, 0.25307,
                    0.41901, 0.41901, 0.67736, 0.67736, 1.02692, 1.02692,
                    1.41335, 1.41335, 1.76498, 1.76498)
  SAS.pred.lwr.new <- c(1.45376, 1.64754, 0.41047, 1.56474, 1.25081)
  SAS.pred.upr.new <- c(1.54510, 1.74598, 0.49958, 1.65829, 1.34300)
  
  ## Expectations for original data
  expect_true(all.equal(round(DNase1.conf[, "se.fit"], 6), SAS.se.fit))
  expect_true(all.equal(round(DNase1.pred[, "se.fit"], 6), SAS.se.fit))
  expect_true(all.equal(round(DNase1.conf[, "lwr"], 5), SAS.conf.lwr, tol = 1e-05))
  expect_true(all.equal(round(DNase1.conf[, "upr"], 5), SAS.conf.upr, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred[, "lwr"], 5), SAS.pred.lwr, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred[, "upr"], 5), SAS.pred.upr, tol = 1e-05))
  
  ## Expectations for new data
  expect_true(all.equal(round(DNase1.conf2[, "se.fit"], 6), SAS.se.fit.new))
  expect_true(all.equal(round(DNase1.pred2[, "se.fit"], 6), SAS.se.fit.new))
  expect_true(all.equal(round(DNase1.conf2[, "lwr"], 5), SAS.conf.lwr.new, tol = 1e-05))
  expect_true(all.equal(round(DNase1.conf2[, "upr"], 5), SAS.conf.upr.new, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred2[, "lwr"], 5), SAS.pred.lwr.new, tol = 1e-05))
  expect_true(all.equal(round(DNase1.pred2[, "upr"], 5), SAS.pred.upr.new, tol = 1e-05))
  
})


test_that("", {

  # Simulate some data
  set.seed(101)  # for reproducibilty
  x <- rep(1:10, each = 3)
  y <- 1 + 2 * x + rnorm(length(x), sd = 1)
  d <- data.frame("x" = x, "y" = y)

  # Fit some linear models
  lm1 <- lm(y ~ x, data = d)
  lm2 <- lm(y ~ x)

  # Predictions only
  pred.investr <- predFit(lm1)
  pred.stats <- predict(lm1)

  # Predictions and confidence intervals
  pred.investr.se.conf <- predFit(lm1, interval = "confidence")
  pred.stats.se.conf <- predict(lm1, interval = "confidence")

  # Predictions and prediction intervals
  pred.investr.se.pred <- predFit(lm1, interval = "prediction")
  pred.stats.se.pred <- predict(lm1, interval = "prediction")

  # Predictions and standard errors
  pred.investr.se <- predFit(lm1, se.fit = TRUE)
  pred.stats.se <- predict(lm1, se.fit = TRUE)

  # Predictions, confidence intervals, and standard errors
  pred.investr.se.conf <- predFit(lm1, se.fit = TRUE, interval = "confidence")
  pred.stats.se.conf <- predict(lm1, se.fit = TRUE, interval = "confidence")

  # Predictions, prediction intervals, and standard errors
  pred.investr.se.pred <- predFit(lm1, se.fit = TRUE, interval = "prediction")
  pred.stats.se.pred <- predict(lm1, se.fit = TRUE, interval = "prediction")

  # Expectations
  expect_equal(pred.investr, pred.stats)
  expect_equal(pred.investr.conf, pred.stats.conf)
  expect_equal(pred.investr.pred, pred.stats.pred)
  expect_equal(pred.investr.se$se.fit, pred.stats.se$se.fit)
  expect_equal(pred.investr.se.conf$se.fit, pred.stats.se$se.fit)
  expect_equal(pred.investr.se.pred$se.fit, pred.stats.se$se.fit)
  expect_equal(pred.investr.se.conf$fit, pred.stats.se.conf$fit)
  expect_equal(pred.investr.se.pred$fit, pred.stats.se.pred$fit)

  # Using predFit on an object with no data component should cause an error if no
  # data frame is supplied via the newdata argument.
  expect_error(predFit(lm2))
  
})
