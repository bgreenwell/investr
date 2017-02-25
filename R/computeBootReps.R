#' @keywords internal
computeBootReps <- function(object, nsim, seed, boot.type, yname,
                            mean.response, y0, x0.name, lower, upper, 
                            extendInt, tol, maxiter, progress) {
  UseMethod("computeBootReps")
}


#' @keywords internal
computeBootReps.lm <- function(object, nsim, seed, boot.type, yname,
                               mean.response, y0, x0.name, lower, upper, 
                               extendInt, tol, maxiter, progress) {
  
  # Sanity check
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  
  # Set up progress bar (if requested)
  if (progress) { 
    pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)
  }
  
  # Initialize random number generator
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    stats::runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  # Simulate new response vectors
  ftd <- stats::fitted(object)  # fitted values
  res <- stats::residuals(object) # redisuals
  if (boot.type == "parametric") {  
    ss <- stats::simulate(object, nsim = nsim)
  } else {
    ss <- replicate(nsim, ftd + sample(res, replace = TRUE), simplify = FALSE)
  }
  
  # Bootstrap function
  x0Fun <- function(i) {
    
    # Update model using simulated response data
    boot.data <- eval(object$call$data)  # copy data
    boot.data[, yname] <- ss[[i]]  # simulated response vector
    boot.object <- tryCatch(stats::update(object, data = boot.data),
                            error = function(e) NULL)
    
    # If updating the model fails, then return value is NA
    if (is.null(boot.object)) {
      ret <- NA
    } else {
      
      # Simulate new response (different from simulated response vector)
      if (mean.response) {  # regulation
        y0.star <- y0  # hold constant in bootstrap replications
      } else {  # calibration
        if (boot.type == "parametric") {
          y0.star <- y0 + stats::rnorm(length(y0), sd = stats::sigma(object))
        } else {
          y0.star <- y0 + sample(res, size = length(y0), replace = TRUE)
        }
      }
      
      # Calculate point estimate
      ret <- tryCatch(stats::uniroot(function(x) {
        stats::predict(boot.object, newdata = makeData(x, x0.name)) - mean(y0.star)
      }, interval = c(lower, upper), extendInt = extendInt, 
      tol = tol, maxiter = maxiter)$root, 
      error = function(e) NA)
    }
    
    # Update progress bar
    if (progress) { 
      utils::setTxtProgressBar(pb, i) 
    }
    
    # Return estimate
    ret
    
  }
  
  # Calculate bootstrap replicates
  x0.star <- sapply(seq_len(nsim), x0Fun)
  
  # Check for errors and return the runs that did not fail
  if (anyNA(x0.star)) {
    num.fail <- sum(is.na(x0.star))
    warning("some bootstrap runs failed (", num.fail, "/", nsim, 
            ")")
    x0.star <- stats::na.omit(x0.star)  # remove runs that failed
    attr(x0.star, "bootFail") <- num.fail  # remove attributes
  } else {
    num.fail <- NULL
  }

  # Return bostrap replicates
  x0.star
  
}


#' @keywords internal
computeBootReps.nls <- function(object, nsim, seed, boot.type, yname,
                                mean.response, y0, x0.name, lower, upper, 
                                extendInt, tol, maxiter, progress) {
  
  # Sanity check
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  
  # Set up progress bar (if requested)
  if (progress) { 
    pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)
  }
  
  # Initialize random number generator
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    stats::runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  # Simulate new response vectors
  ftd <- stats::fitted(object)  # fitted values
  res <- stats::residuals(object) # redisuals
  if (boot.type == "parametric") {  
    ss <- stats::simulate(object, nsim = nsim)
  } else {
    ss <- replicate(nsim, ftd + sample(res, replace = TRUE), simplify = FALSE)
  }
  
  # Bootstrap function
  x0Fun <- function(i) {
    
    # Update model using simulated response data
    boot.data <- eval(object$call$data)  # copy data
    boot.data[, yname] <- ss[[i]]  # simulated response vector
    boot.object <- tryCatch(stats::update(object, data = boot.data),
                            error = function(e) NULL)
    
    # If updating the model fails, then return value is NA
    if (is.null(boot.object)) {
      ret <- NA
    } else {
      
      # Simulate new response (different from simulated response vector)
      if (mean.response) {  # regulation
        y0.star <- y0  # hold constant in bootstrap replications
      } else {  # calibration
        if (boot.type == "parametric") {
          y0.star <- y0 + stats::rnorm(length(y0), sd = stats::sigma(object))
        } else {
          y0.star <- y0 + sample(res, size = length(y0), replace = TRUE)
        }
      }
      
      # Calculate point estimate
      ret <- tryCatch(stats::uniroot(function(x) {
        stats::predict(boot.object, newdata = makeData(x, x0.name)) - mean(y0.star)
      }, interval = c(lower, upper), extendInt = extendInt, 
      tol = tol, maxiter = maxiter)$root, 
      error = function(e) NA)
    }
    
    # Update progress bar
    if (progress) { 
      utils::setTxtProgressBar(pb, i) 
    }
    
    # Return estimate
    ret
    
  }
  
  # Calculate bootstrap replicates
  x0.star <- sapply(seq_len(nsim), x0Fun)
  
  # Check for errors and return the runs that did not fail
  if (anyNA(x0.star)) {
    num.fail <- sum(is.na(x0.star))
    warning("some bootstrap runs failed (", num.fail, "/", nsim, 
            ")")
    x0.star <- stats::na.omit(x0.star)  # remove runs that failed
    attr(x0.star, "bootFail") <- num.fail  # remove attributes
  } else {
    num.fail <- NULL
  }
  
  # Return bostrap replicates
  x0.star
  
}
