# Inverse estimation for Kaplan-Meier survival curves (issue #30). Note:
# survival is deliberately not attached via library(), since it would mask
# investr's bladder data set.

exit_if_not(requireNamespace("survival", quietly = TRUE))

km <- survival::survfit(survival::Surv(time, status) ~ 1,
                        data = survival::aml)

# Median survival time matches survival's own quantile.survfit inversion
res <- invest(km)
q <- stats::quantile(km, probs = 0.5, conf.int = TRUE)
expect_inherits(res, "invest")
expect_identical(res$interval, "inversion")
expect_equal(res$estimate, unname(q$quantile))
expect_equal(res$lower, unname(q$lower))
expect_equal(res$upper, unname(q$upper))

# Known values for the aml data (matches print.survfit's median and 0.95 CI)
expect_equal(res$estimate, 27)
expect_equal(res$lower, 18)
expect_equal(res$upper, 48)

# Other survival probabilities invert consistently
res75 <- invest(km, y0 = 0.75)
q75 <- stats::quantile(km, probs = 0.25, conf.int = TRUE)
expect_equal(res75$estimate, unname(q75$quantile))

# The Surv method fits the curve internally and matches the survfit method
s <- survival::Surv(survival::aml$time, survival::aml$status)
expect_equal(invest(s), res)

# The Surv method honors level directly (narrower interval at 90%)
res90 <- invest(s, level = 0.90)
expect_equal(res90$estimate, res$estimate)
expect_true(res90$upper <= res$upper)

# A survfit fit at a different conf.int than level must error informatively
expect_error(invest(km, level = 0.90), pattern = "Refit")

# Multiple strata are not supported (invert each stratum separately)
km2 <- survival::survfit(survival::Surv(time, status) ~ x,
                         data = survival::aml)
expect_error(invest(km2), pattern = "single survival curve")
res_s1 <- invest(km2[1])
expect_inherits(res_s1, "invest")
expect_true(is.numeric(res_s1$estimate) && !is.na(res_s1$estimate))

# y0 must be a single probability in (0, 1)
expect_error(invest(km, y0 = 1.5))
expect_error(invest(km, y0 = 0))
expect_error(invest(km, y0 = c(0.25, 0.5)))

# A curve that never drops below y0 warns and returns NA
expect_warning(invest(km, y0 = 0.01))
res_na <- suppressWarnings(invest(km, y0 = 0.01))
expect_true(is.na(res_na$estimate))
