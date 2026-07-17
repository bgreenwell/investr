# simulate.nls works as expected

mod <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
           start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
           data = nasturtium)

sim1 <- investr:::simulate.nls(mod, nsim = 10L, seed = 101)
sim2 <- investr:::simulate.nls(mod, nsim = 10L, seed = 101)

expect_inherits(sim1, "data.frame")
expect_equal(dim(sim1), c(42, 10))
expect_identical(sim1, sim2)
