.investrEnv <- new.env()  # FIXME: Is this needed anywhere?

##' Nonparametric bootstrap for (controlled) calibration
##' 
##' Perform nonparametric bootstrap for controlled calibration.
bootCal <- function(object, ...) {
  UseMethod("bootCal")
} 

##' @rdname bootCal
##' @export
##' @method bootCal lm
bootCal.lm <- function(object, y0, nsim = 99, mean.response = FALSE, 
                       seed = NULL, lower, upper, tol = 1e-10, maxiter = 1000) 
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (!is.null(seed)) set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  
  vnms <- getVars(object, return.data = TRUE)
  xname <- vnms$x.names  # covariate name
  yname <- vnms$y.names  # response name
  data <- vnms$data
  if (missing(lower)) assign("lower", min(vnms$x), envir = .investrEnv)
  if (missing(upper)) assign("upper", max(vnms$x), envir = .investrEnv)
  
  ## Classical estimator function
  xFun <- function(object, y) {
    uniroot(function(x) {
      predict(object, newdata = setNames(list(x), xname)) - y
    }, tol = tol, maxiter = maxiter,
    lower = get("lower", envir = .investrEnv), 
    upper = get("upper", envir = .investrEnv))$root
  }
  t0 <- invest(object, y0 = y0, interval = "none")  # original estimate
  
  ## Bootstrap setup
  n <- length(res <- resid(object))  # sample size and residuals
  res <- res / sqrt(1 - hatvalues(object))  # scaled residuals
  res <- res - mean(res)  # cetered residuals
  boot.data <- data.frame(data, res = res, fit = fitted(object))
  boot.fun <- function(.data, .ind) {
    .data[, yname] <- .data$fit + .data$res[.ind]
    boot.object <- update(object, data = .data)
    ## Simulate the correct variance
    Y0 <- y0 + sample(.data$res, size = length(y0), replace = TRUE)
    ## Make sure the original estimate also gets returned
    if (all(.ind == 1:n)) t0 else xFun(boot.object, y = Y0)
  }
  
  ## Bootstrap replicates
  t.star <- replicate(nsim, {
    boot.fun(boot.data, sample(n, replace = TRUE))
  })
  
  ## Bootstrap object
  mc <- match.call()
  mc[[1L]] <- "boot"
  structure(list(t0 = t0, 
                 t = as.matrix(t.star), 
                 R = nsim, 
                 data = boot.data, 
                 seed = .Random.seed, 
                 statistic = boot.fun, 
                 sim = "ordinary", 
                 stype = "i",
                 strata = NULL,#rep(1, n),
                 weights = NULL,
                 call = mc), class = "boot")
}

crystal.lm <- lm(weight ~ time, data = crystal)
invest(crystal.lm, y0 = 8, interval = "Wald")
invest(crystal.lm, y0 = 8, interval = "inversion")

res <- bootCal(crystal.lm, y0 = 8, nsim = 99)
boot.ci(res)
hist(res$t, freq=F, 50)
quantile(res$t, c(0.025, 0.975))
sd(out)


