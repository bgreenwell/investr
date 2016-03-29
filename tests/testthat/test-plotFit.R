context("Basic functionality")

library(MASS)

test_that("plotFit works", {
  
  # Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal_lm <- lm(weight ~ time, data = crystal)
  crystal_rlm <- rlm(weight ~ time, data = crystal)  # requires MASS
  crystal_lqs <- lqs(weight ~ time, data = crystal)  # requires MASS

  # Dobson's beetle example
  beetle_glm <- glm(cbind(y, n-y) ~ ldose, data = beetle, 
                    family = binomial(link = "cloglog"))
  
  # Nasturtium example
  nas_nls <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
                 data = nasturtium)
  
  # Simulated data
  x <- rnorm(10)
  y <- rnorm(10)
  
  # Expectations
  expect_error(plotFit(lm(y ~ x)))
  expect_silent(plotFit(crystal_rlm))  # requires MASS
  expect_silent(plotFit(crystal_lqs))  # requires MASS
  expect_silent(plotFit(crystal_lm))
  expect_warning(plotFit(crystal_lm, interval = "both"))
  expect_warning(plotFit(crystal_lm, interval = "both", extend.range = TRUE, shade = TRUE))
  expect_silent(plotFit(crystal_lm, interval = "both", hide = FALSE))
  expect_silent(plotFit(crystal_lm, interval = "both", extend.range = TRUE, shade = TRUE, hide = FALSE))
  
  expect_silent(plotFit(beetle_glm))
  expect_silent(plotFit(beetle_glm, interval = "confidence"))
  expect_silent(plotFit(beetle_glm, interval = "confidence", hide = FALSE))
  expect_silent(plotFit(beetle_glm, interval = "confidence", extend.range = TRUE, shade = TRUE, hide = FALSE))
  
  expect_silent(plotFit(nas_nls))
  expect_silent(plotFit(nas_nls, interval = "both"))
  expect_silent(plotFit(nas_nls, interval = "both", 
                        extend.range = TRUE, shade = TRUE, xlim = c(1, 4)))
  expect_silent(plotFit(nas_nls, interval = "both", hide = FALSE))
  expect_silent(plotFit(nas_nls, interval = "both", extend.range = TRUE, 
                        shade = TRUE, hide = FALSE, xlim = c(1, 4)))
  
})
