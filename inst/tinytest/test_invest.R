# Setup ------------------------------------------------------------------------

exit_if_not(requireNamespace("nlme", quietly = TRUE))
library(nlme)

# Simulate some data for multiple linear regression
set.seed(101)
x1 <- runif(50, min = 0, max = 10)
x2 <- runif(50, min = 0, max = 10)
x3 <- gl(2, 50, labels = c("Control", "Treat"))
e <- rnorm(50, sd = 3)
y <- ifelse(x3 == "Control", 3 + 2*x1 - 3*x2, 6 + 7*x1 - 3*x2) + e
d <- data.frame("x1" = x1, "x2" = x2, "x3" = x3, "y" = y)

# DNase data from the datasets package
data(DNase, package = "datasets")
DNase1 <- subset(DNase, Run == 1)

# Bladder volume data
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

# Crystal weight example from Graybill and Iyer (1996, p. 434)
crystal_lm <- lm(weight ~ time, data = crystal)
crystal_glm <- glm(weight ~ time, data = crystal, family = gaussian)

# Simulated multiple linear regression model
multi_lm <- lm(y ~ x1 + x2 + x3 + x1 * x3, data = d)

# Log-logistic model for nasturtium data
nas_nls <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
               start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
               data = nasturtium)

# Dobson's beetle data
beetle_glm <- glm(cbind(y, n-y) ~ ldose, data = beetle, family = "binomial")

# Random intercept and slope model
bladder_lme <- lme(HD^(3/2) ~ volume,
                   random = list(subject = pdDiag(~volume)), data = bladder)

# Inverse estimation with linear models -----------------------------------------

# invest and calibrate produce the same results
res1.cal <- calibrate(crystal_lm, y0 = 5)
res2.cal <- calibrate(crystal_lm, y0 = 5, mean.response = TRUE)
res3.cal <- calibrate(crystal_lm, y0 = 5, interval = "Wald")
res4.cal <- calibrate(crystal_lm, y0 = 5, interval = "Wald", mean.response = TRUE)
res1.inv <- invest(crystal_lm, y0 = 5)
res2.inv <- invest(crystal_lm, y0 = 5, mean.response = TRUE)
res3.inv <- invest(crystal_lm, y0 = 5, interval = "Wald")
res4.inv <- invest(crystal_lm, y0 = 5, interval = "Wald", mean.response = TRUE)

expect_equal(calibrate(crystal_lm, y0 = 5, interval = "none"),
             unname(invest(crystal_lm, y0 = 5, interval = "none")),
             tol = 1e-05)
expect_equal(res1.cal, res1.inv, tol = 1e-05)
expect_equal(res2.cal, res2.inv, tol = 1e-05)
expect_equal(res3.cal, res3.inv, tol = 1e-05)
expect_equal(res4.cal, res4.inv, tol = 1e-04)

# invest should throw an error when extrapolating beyond the range of the data
expect_silent(calibrate(crystal_lm, y0 = 20, interval = "none"))  # calibrate should still work!
expect_error(invest(crystal_lm, y0 = 2, interval = "inversion"))  # lwr limit should fail
expect_error(invest(crystal_lm, y0 = 14, interval = "inversion"))  # upr limit should fail

# bootstrap method produces reasonable results
boot.npar <- invest(crystal_lm, y0 = 5, data = crystal, interval = "percentile",
                    boot.type = "nonparametric", nsim = 9, seed = 101)
boot.par <- invest(crystal_lm, y0 = 5, data = crystal, interval = "percentile",
                   boot.type = "parametric", nsim = 9, seed = 101)

expect_inherits(boot.npar, c("invest", "bootCal"))
expect_inherits(boot.par, c("invest", "bootCal"))
expect_silent(plot(boot.npar))

# Percentile-interval results are also valid "boot" objects, usable with
# the boot package (regression tests for issue #32)
expect_inherits(boot.npar, "boot")
expect_inherits(boot.par, "boot")
expect_identical(boot.npar$t0, boot.npar$estimate)
expect_identical(boot.npar$t, matrix(boot.npar$bootreps, ncol = 1L))
expect_identical(boot.npar$R, length(boot.npar$bootreps))
if (requireNamespace("boot", quietly = TRUE)) {
  boo32 <- invest(crystal_lm, y0 = 5, interval = "percentile", nsim = 199,
                  seed = 101)
  ci32 <- boot::boot.ci(boo32, type = c("norm", "basic", "perc"))
  expect_inherits(ci32, "bootci")
  # boot.ci's percentile interval uses the (R+1) order-statistic convention
  # rather than stats::quantile()'s type-7 default, so allow slack
  expect_equal(ci32$percent[4:5], c(boo32$lower, boo32$upper), tol = 0.1)
}

# multiple predictor results match output from JMP (v11)
contr <- data.frame("x1" = mean(x1), "x3" = "Control")
treat <- data.frame("x1" = mean(x1), "x3" = "Treat")
inv.contr <- invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = contr, lower = -10)
inv.treat <- invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = treat)

b <- unname(coef(multi_lm))

expect_equal(inv.contr$estimate,
             (20 - b[1L] - b[2L]*mean(x1)) / b[3L])
expect_equal(inv.treat$estimate,
             (20 - b[1L] - (b[5L] + b[2L])*mean(x1) - b[4L]) / b[3L])
expect_equal(inv.contr$lower, -4.28662, tol = 1e-04)
expect_equal(inv.contr$upper, -0.35969, tol = 1e-04)
expect_equal(inv.treat$lower, 5.13027, tol = 1e-04)
expect_equal(inv.treat$upper, 8.94294, tol = 1e-04)

expect_error(invest(multi_lm, y0 = 20, interval = "inversion",
                    newdata = contr, lower = -10))
expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = data.matrix(treat), lower = -10))
expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = d, lower = -10))
expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = d[, 1, drop = FALSE], lower = -10))
expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = setNames(treat, c("y1", "y2")), lower = -10))

# newdata columns must have classes compatible with the fitted model's data
# (regression tests for issue #36)

expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = data.frame(x1 = mean(x1), x3 = "Bogus"),
                    lower = -10))
expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = data.frame(x1 = mean(x1), x3 = 1),
                    lower = -10))
expect_error(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                    newdata = data.frame(x1 = "five", x3 = "Control"),
                    lower = -10))
expect_silent(invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
                     newdata = contr, lower = -10))

# The "Use plotFit for guidance" hint on a failed confidence-limit search
# must not appear for multi-predictor models, since plotFit() doesn't
# support more than one predictor variable (regression test for issue #35)
err_msg <- tryCatch(
  invest(multi_lm, y0 = 20, x0.name = "x2", interval = "inversion",
         newdata = contr, lower = inv.contr$estimate - 0.01,
         upper = inv.contr$estimate + 0.01),
  error = function(e) conditionMessage(e)
)
expect_true(grepl("confidence limit not found", err_msg))
expect_false(grepl("plotFit", err_msg))

# Inverse estimation with nonlinear models --------------------------------------

# approximate standard error is correct
# The estimate based on the deltaMethod function from car is 0.2847019; this
# was calculated in the R Journal article.
se1 <- invest(nas_nls, y0 = c(309, 296, 419), interval = "Wald")$se
se2 <- invest(nas_nls, y0 = c(309, 296, 419), interval = "Wald",
              data = nasturtium, tol = 1e-10)$se

expect_equal(se1, 0.2847019, tol = 1e-05)  # less precise
expect_equal(se2, 0.2847019, tol = 1e-05)  # more precise

# Wald and inversion methods produce the same point estimate
res.inv <- invest(nas_nls, y0 = c(309, 296, 419), interval = "inversion")
res.wald <- invest(nas_nls, y0 = c(309, 296, 419), interval = "Wald")

expect_identical(res.inv$estimate, res.wald$estimate)

# bootstrap produces reasonable results
expect_silent(invest(nas_nls, y0 = c(309, 296, 419), interval = "percentile",
                     nsim = 9, seed = 101))

# invest works properly on 'special' nls fits
DNase1_nls_1 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                  data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))

# Using conditional linearity
DNase1_nls_2 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
                    data = DNase1,
                    start = list(xmid = 0, scal = 1),
                    algorithm = "plinear")

# Using selfStart
DNase1_nls_3 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal),
                    data = DNase1)

# Using Port's nl2sol algorithm
DNase1_nls_4 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                    data = DNase1,
                    start = list(Asym = 3, xmid = 0, scal = 1),
                    algorithm = "port")

res_1 <- invest(DNase1_nls_1, y0 = 0.5)
res_3 <- invest(DNase1_nls_3, y0 = 0.5)
res_4 <- invest(DNase1_nls_4, y0 = 0.5)

expect_error(invest(DNase1_nls_2, y0 = 0.5))
expect_equal(res_1$estimate, res_3$estimate, tol = 1e-05)
expect_equal(res_1$lower, res_3$lower, tol = 1e-05)
expect_equal(res_1$upper, res_3$upper, tol = 1e-05)
expect_equal(res_1$estimate, res_4$estimate, tol = 1e-05)
expect_equal(res_1$lower, res_4$lower, tol = 1e-05)
expect_equal(res_1$upper, res_4$upper, tol = 1e-05)

# Inverse estimation with generalized linear models ------------------------------

# inversion and Wald methods work
# Based on methodology from Page 207 from Categorical Data Analysis, 2nd
# Edition, by Alan Agresti.
res <- invest(beetle_glm, y0 = 0.5, interval = "inversion", tol = 1e-10)
a <- unname(coef(beetle_glm)[1])
b2 <- unname(coef(beetle_glm)[2])
var_a <- vcov(beetle_glm)[1, 1]
var_b <- vcov(beetle_glm)[2, 2]
cov_ab <- vcov(beetle_glm)[1, 2]
fun <- function(x, p = 0.5) {
  abs(a + b2*x - qlogis(p)) / sqrt(var_a + x^2*var_b + 2*x*cov_ab) -
    qnorm(0.975)
}
lwr <- uniroot(fun, lower = 1.76, upper = 1.77, tol = 1e-10)$root
upr <- uniroot(fun, lower = 1.77, upper = 1.79, tol = 1e-10)$root
expect_equal(res$lower, lwr)
expect_equal(res$upper, upr)

# Check Taylor series approximation of standard error using MASS::dose.p
wald_se <- invest(beetle_glm, y0 = 0.5, interval = "Wald")$se
expect_equal(wald_se, 0.003858052, tol = 1e-05)

# ?MASS::dose.p
#
#                 Dose        SE
#   p = 0.25: 2.231265 0.2499089
#   p = 0.50: 3.263587 0.2297539
#   p = 0.75: 4.295910 0.2746874
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20 - numdead)
budworm <- data.frame(SF, sex, ldose)
budworm_glm <- glm(SF ~ sex + ldose - 1, family = binomial, data = budworm)
p.25 <- invest(budworm_glm, y0 = 1/4, x0.name = "ldose",
               newdata = data.frame(sex = "F"), interval = "Wald")
p.50 <- invest(budworm_glm, y0 = 1/2, x0.name = "ldose",
               newdata = data.frame(sex = "F"), interval = "Wald")
p.75 <- invest(budworm_glm, y0 = 3/4, x0.name = "ldose",
               newdata = data.frame(sex = "F"), interval = "Wald")

expect_error(invest(budworm_glm, y0 = 1/2, newdata = data.frame(sex = "F")))  # missing x0.name
expect_error(invest(budworm_glm, y0 = 1/2, x0.name = "ldose"))  # missing newdata
expect_error(invest(budworm_glm, y0 = 1/2, x0.name = "ldose",
             newdata = data.matrix(data.frame(sex = "F"))))  # newdata must be a data frame
expect_equal(p.25$estimate, 2.231265, tol = 1e-05)
expect_equal(p.50$estimate, 3.263587, tol = 1e-05)
expect_equal(p.75$estimate, 4.295910, tol = 1e-05)
expect_equal(p.25$se, 0.2499089, tol = 1e-05)
expect_equal(p.50$se, 0.2297539, tol = 1e-05)
expect_equal(p.75$se, 0.2746874, tol = 1e-05)
expect_error(invest(beetle_glm, y0 = 0.5, lower = 1.770, upper = 1.800))  # lwr should fail
expect_error(invest(beetle_glm, y0 = 0.5, lower = 1.700, upper = 1.772))  # upr should fail
expect_error(invest(beetle_glm, y0 = 0.5, interval = "percentile"))  # bootstrap should always fail

# Inverse estimation for Gaussian GLM (closely) matches LM results
glm_inversion <- invest(crystal_glm, y0 = 5, interval = "inversion")
glm_wald <- invest(crystal_glm, y0 = 5, interval = "Wald")

lm_inversion <- invest(crystal_lm, y0 = 5, interval = "inversion", mean.response = TRUE)
lm_wald <- invest(crystal_lm, y0 = 5, interval = "Wald", mean.response = TRUE)

# The GLM method uses a critical value based on the standard normal
# distribution; hence, these intervals will be too optimistic. However, the
# point estimate should still be the same.
expect_equal(glm_inversion$estimate, lm_inversion$estimate)
expect_equal(glm_wald$estimate, lm_wald$estimate)
expect_true(glm_inversion$lower > lm_inversion$lower)
expect_true(glm_inversion$upper < lm_inversion$upper)
expect_true(glm_wald$lower > lm_wald$lower)
expect_true(glm_wald$upper < lm_wald$upper)

# Inverse estimation with linear mixed-effects models ----------------------------

# inversion method works
res <- invest(bladder_lme, y0 = 500, interval = "inversion")

expect_equal(res$estimate, 8.015521, tol = 1e-05)
expect_equal(res$lower, 4.227962, tol = 1e-05)  # 4.227965
expect_equal(res$upper, 11.91918, tol = 1e-05)

# Wald method works
tvals <- qt(c(0.025, 0.975), length(resid(bladder_lme)) - 1)
res.norm <- invest(bladder_lme, y0 = 500, interval = "Wald")
res.t <- invest(bladder_lme, y0 = 500, interval = "Wald", q1 = tvals[1],
                q2 = tvals[2])

# True values calculated by hand
ci.norm <- 8.015521 + qnorm(c(0.025, 0.975))*1.954191
ci.t <- 8.015521 + tvals*1.954191

expect_equal(res.norm$estimate, 8.015521, tol = 1e-05)  # estimate
expect_equal(res.t$estimate, 8.015521, tol = 1e-05)     # estimate
expect_equal(res.norm$se, 1.954191, tol = 1e-05)        # SE
expect_equal(res.t$se, 1.954191, tol = 1e-05)           # SE
expect_equal(res.norm$lower, ci.norm[1], tol = 1e-05)   # lower limit
expect_equal(res.t$lower, ci.t[1], tol = 1e-05)         # lower limit
expect_equal(res.norm$upper, ci.norm[2], tol = 1e-05)   # upper limit
expect_equal(res.t$upper, ci.t[2], tol = 1e-05)         # upper limit

# Models fit inside a function must not require their training data to be
# reachable from invest()'s caller (regression tests for issues #41, #42,
# #45: .data was reconstructed via eval(object$call$data, parent.frame()),
# which fails once the fitting function's local environment is gone)

# #42: glm fit with a data= argument that's a local variable
fit_glm_in_fun <- function() {
  dat <- beetle
  glm(cbind(y, n-y) ~ ldose, data = dat, family = binomial(link = "cloglog"))
}
mod42 <- fit_glm_in_fun()
mod42_ref <- glm(cbind(y, n-y) ~ ldose, data = beetle, family = binomial(link = "cloglog"))
expect_silent(res42 <- invest(mod42, y0 = 0.5))
expect_equal(res42$estimate, invest(mod42_ref, y0 = 0.5)$estimate, tol = 1e-05)

# #41: glm fit with no data= argument at all (variables taken from the
# environment directly)
fit_glm_no_data <- function() {
  x1_ <- rep(1:10, each = 3)
  y1_ <- rbinom(30, 1, plogis(x1_ - 5))
  glm(y1_ ~ x1_, family = binomial)
}
mod41 <- fit_glm_no_data()
expect_silent(res41 <- invest(mod41, y0 = 0.5, interval = "none"))
expect_true(is.numeric(res41))

# #45: nls fit with a data= argument that's a local variable
fit_nls_in_fun <- function() {
  dat <- nasturtium
  nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
      start = list(theta1 = 1000, theta2 = -1, theta3 = 1), data = dat)
}
mod45 <- fit_nls_in_fun()
expect_silent(res45 <- invest(mod45, y0 = c(309, 296, 419), interval = "inversion"))
expect_equal(res45$estimate, invest(nas_nls, y0 = c(309, 296, 419), interval = "inversion")$estimate,
             tol = 1e-05)

# lm fit with a data= argument that's a local variable
fit_lm_in_fun <- function() {
  dat <- crystal
  lm(weight ~ time, data = dat)
}
mod_lm_local <- fit_lm_in_fun()
expect_silent(res_lm_local <- invest(mod_lm_local, y0 = 5))
expect_equal(res_lm_local$estimate, invest(crystal_lm, y0 = 5)$estimate, tol = 1e-05)
