# plotFit works

exit_if_not(requireNamespace("MASS", quietly = TRUE))
library(MASS)

# Crystal weight example from Graybill and Iyer (1996, p. 434)
crystal_lm <- lm(weight ~ time, data = crystal)
crystal_rlm <- rlm(weight ~ time, data = crystal)  # requires MASS
crystal_lqs <- lqs(weight ~ time, data = crystal)  # requires MASS

# Nasturtium example
nas_nls <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
               start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
               data = nasturtium)

# Simulated data
set.seed(101)
x <- rnorm(10)
y <- rnorm(10)

expect_error(plotFit(lm(y ~ x)))
expect_error(plotFit(rlm(crystal$weight ~ crystal$time)))  # requires MASS
expect_error(plotFit(lqs(crystal$weight ~ crystal$time)))  # requires MASS
expect_silent(plotFit(crystal_rlm))  # requires MASS
expect_silent(plotFit(crystal_lqs))  # requires MASS
expect_silent(plotFit(crystal_lm))

expect_silent(plotFit(nas_nls))
expect_silent(plotFit(nas_nls, interval = "both"))
expect_silent(plotFit(nas_nls, interval = "both",
                      extend.range = TRUE, shade = TRUE, xlim = c(1, 4)))
expect_silent(plotFit(nas_nls, interval = "both", hide = FALSE))
expect_silent(plotFit(nas_nls, interval = "both", extend.range = TRUE,
                      shade = TRUE, hide = FALSE, xlim = c(1, 4)))
