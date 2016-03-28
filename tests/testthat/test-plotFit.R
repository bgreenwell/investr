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
  expect_silent(plotFit(beetle_glm, pch = 19, cex = 1.2, lwd = 2, 
                        xlab = "Log dose of carbon disulphide",
                        interval = "confidence", shade = TRUE, 
                        col.conf = "lightskyblue"))
  expect_silent(plotFit(nas_nls, lwd.fit = 2))
  
})
