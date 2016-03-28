context("Basic functionality")


test_that("plotFit works", {
  
  # Crystal weight example from Graybill and Iyer (1996, p. 434)
  crystal_lm <- lm(weight ~ time, data = crystal)
  
  # Dobson's beetle example
  beetle_glm <- glm(cbind(y, n-y) ~ ldose, data = beetle, 
                    family = binomial(link = "cloglog"))

  # Nasturtium example
  nas_nls <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
                 start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
                 data = nasturtium)
  
  # Expectations
  expect_silent(plotFit(crystal_lm))
  expect_warning(plotFit(crystal_lm, interval = "both"))
  expect_warning(plotFit(crystal_lm, interval = "both", extend.range = TRUE, shade = TRUE))
  expect_silent(plotFit(crystal_lm, interval = "confidence", hide = FALSE))
  expect_silent(plotFit(crystal_lm, interval = "confidence", extend.range = TRUE, shade = TRUE, hide = FALSE))

  expect_silent(plotFit(beetle_glm))
  expect_silent(plotFit(beetle_glm, interval = "confidence"))
  expect_silent(plotFit(beetle_glm, interval = "confidence", hide = FALSE))
  expect_silent(plotFit(beetle_glm, interval = "confidence", extend.range = TRUE, shade = TRUE, hide = FALSE))

  expect_silent(plotFit(nas_nls))
  expect_warning(plotFit(nas_nls, interval = "both"))
  expect_warning(plotFit(nas_nls, interval = "both", extend.range = TRUE, shade = TRUE))
  expect_silent(plotFit(nas_nls, interval = "confidence", hide = FALSE))
  expect_silent(plotFit(nas_nls, interval = "confidence", extend.range = TRUE, shade = TRUE, hide = FALSE))
  
})
