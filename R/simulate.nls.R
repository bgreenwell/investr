##' @keywords internal
simulate.nls <- function (object, nsim = 1, seed = NULL, ...) {
  
  ## Set seed value
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  ## Mostly code from simulate
  ftd <- fitted(object)
  nm <- names(ftd)
  n <- length(ftd)
  ntot <- n * nsim
  vars <- deviance(object)/df.residual(object)
  if (!is.null(object$weights)) vars <- vars/object$weights
  val <- ftd + rnorm(ntot, sd = sqrt(vars))
  if (!is.list(val)) {
    dim(val) <- c(n, nsim)
    val <- as.data.frame(val)
  } else {
    class(val) <- "data.frame"
  }
  names(val) <- paste("sim", seq_len(nsim), sep = "_")
  if (!is.null(nm)) 
    row.names(val) <- nm
  attr(val, "seed") <- RNGstate
  val
  
}
