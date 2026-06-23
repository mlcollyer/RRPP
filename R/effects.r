
#' Obtain P-value from a vector of statistics
#'
#' A function to find the probability of values greater or lesser than target,
#' from a vector of statistics presumably obtained in random permutations.
#'
#' @param s The sampling distribution vector to use.
#' @param target The value to target in the distribution.  (If null, the first value
#' in the vector is used.).  If the target exists outside the range of s,
#' a probability of 0 or 1 is certain.
#' @param greater Logical value for whether the probability should be "greater than 
#' or equal to".  Change to greater = FALSE for "less than or equal to".
#' @export
#' @author Michael Collyer
#' @keywords utilities
pval <- function(s, target = NULL, greater = TRUE){
  if(is.null(target)) target <- s[1]
  p <- length(s)
  pv <- if(greater) 
    length(which(s >= target)) else 
      length(which(s <= target))
  pv / p
}


box.cox.true <- function(y, eps = 0.001){
  
  n <- length(y) - 1
  y.obs <- y[1]
  y <- y[-1]
  
  # Note the code below (to get yy) is short-cut for Jacobian adjustment 
  yy <- y / exp(mean(log(y)))
  logy <- log(yy)
  
  lambda <- seq(-2, 2, 0.001)
  m <- length(lambda)
  
  loglik <- sapply(1:m, function(j){ # same as MASS::boxcox loglik 
    la <- lambda[j]
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  })
  
  lambda.opt <- lambda[which.max(loglik)][[1]]
  
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  y <- c(y.obs, y)
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  
  list(opt.lambda = lambda.opt, transformed = res, lambda = lambda, 
       loglik = loglik)
  
}

box.cox.spline <- function(y, eps = 0.001) {
  
  n <- length(y) - 1
  y.obs <- y[1]
  y <- y[-1]
  
  yy <- y / exp(mean(log(y)))
  logy <- log(yy)
  m <- 20
  lambda <- seq(-1, 1.5, length.out = m)
  
  loglik <- sapply(1:m, function(j){ # same as MASS::boxcox loglik 
    la <- lambda[j]
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  })
  
  lambda.opt <- lambda[which.max(loglik)][[1]]
  
  sp <- spline(lambda, loglik, n = 300)
  lambda.opt <- sp$x[which.max(sp$y)]
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  y <- c(y.obs, y)
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  
  list(opt.lambda = lambda.opt, transformed = res, lambda = sp$x, loglik = sp$y)
  
}

box.cox.interp <- function(y, eps = 0.001) {
  bc <- box.cox.spline(y, eps = eps)
  if(bc$opt.lambda == -1 || bc$opt.lambda == 1.5) 
    bc <- box.cox.true(y, eps = eps)
  return(bc)
}

box.cox.fast <- function(y, eps = 0.001) {
  
  n <- length(y) - 1
  y.obs <- y[1]
  y <- y[-1]

  yy <- y / exp(mean(log(y)))
  logy <- log(yy)
  
  logLik <- function(lambda) {
    la <- lambda
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  }
  
  result <- optimise(logLik, lower = -2, upper = 2, maximum = TRUE)
  lambda.opt <- result$maximum
  
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  y <- c(y.obs, y)
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  list(opt.lambda = lambda.opt, transformed = res, lambda = NULL, loglik = NULL)
}



#' Power Transformation
#'
#' A function that transforms a vector in data to be 
#' approximately normally distributed.
#' 
#' This function uses a Box-Cox power transformation for a vector of values,
#' (presumably statistics from RRPP) with the following conditions:
#' 
#' (1) If any values are 0 or negative, a Yeo-Johnson transformation is performed
#' to optimize lambda with only positive values.
#' 
#' (2) The lambda parameter is optimized over a range of -2 to +2,
#' as is typical for most applications.
#' 
#' (3) The first (or targeted) value (assumed to be the observed value in a 
#' distribution of random values) is removed from the vector to 
#' perform optimization of lambda, but is returned to the vector for
#' transformation.  This assures that large, aberrant values, such
#' as might be found via RRPP for large effects, do not have 
#' leverage for the optimization.
#' 
#' A Yeo-Johnson transformation can be forced, which for positive values
#' merely shifts the distribution by adding 1.  For most RRPP statistics
#' (SS. MS, Rsq, F, d), this is not needed for effect sizes, as values are positive.  
#' However, for log-likelihoods, this might be a good idea for comparisons, 
#' because values might be strongly positive or strongly negative.
#' 
#' Using standard deviates (useStDev) is a good idea for statistics with
#' aberrant ranges (like 0 to 1, or 1,000 to 2,000).  If used, a Yeo-Johnson
#' transformation is implicit, since standardizing introduces negative values.
#' 
#'
#' @param y A vector of values to use.
#' @param eps Tolerance for lambda = 0.
#' @param target The value to target in the distribution.  (If null, the first value
#' in the vector is used.)  This will be moved to the first value in results.
#' @param interp Logical value for whether to spline interpolate lambda 
#' (ideal for vectors of small length).  This will default if the number
#' of values is less than 100, as in \link[MASS]{boxcox}.
#' @param forceYJ A logical value whether to force a Yeo-Johnson transformation for non-negative
#' values.  Negative values trigger
#' the transformation, but positive values can be transformed similarly.
#' @param useStDev A logical value whether to initially use standard deviates of y.
#' This might be useful if the range of y is restricted, say between 0 and 1.
#' @return A list containing the following:
#' \item{opt.lambda}{The optimized lambda parameter.}
#' \item{transformed}{The power-transformed values.}
#' \item{lambda}{If spline interpolation is used, the values of lambda used.}
#' \item{loglik}{If spline interpolation is used, the log-liklihoods calculated for lambda.}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @references Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations 
#' (with discussion). Journal of the Royal Statistical Society B, 26, 211–252.
#' @references Yeo, I. and Johnson, R. (2000) A new family of power transformations 
#' to improve normality or symmetry. Biometrika, 87, 954-959.
powerTrans <- function(y, 
                       eps = 0.001, 
                       target = NULL,
                       interp = FALSE,
                       forcYJ = FALSE,
                       useStDev = FALSE) {
  
  whichNegs <- NULL
  which0 <- NULL
  useYJ <- FALSE
  
  if(!is.null(target))
    y <- c(target, y)
  
  if(forcYJ && !useStDev) {
    useYJ <- TRUE
    y <- abs(y) + 1
  }
    
  if(any(y <= 0)) {
    useYJ <- TRUE
    whichNegs <- which(y < 0)
    which0 <- which(abs(y) <= eps)
    y <- abs(y) + 1
  }
  
  if(useStDev){
    useYJ <- TRUE
    n <- length(y) - 1
    yc <- y - mean(y[-1])
    sy <- sqrt(sum(center(y[-1])^2)/ n)
    y <- yc / sy
    whichNegs <- which(y < 0)
    which0 <- which(abs(y) <= eps)
    y <- abs(y) + 1
  }
  
  n <- length(y)
  if(n < 100) interp <- TRUE
  result <- if(interp) 
    try(box.cox.interp(y, eps = eps), silent = TRUE) else
    box.cox.fast(y, eps = eps)
  if(inherits(result, "try-error"))
    result <- suppressWarnings(
      box.cox.fast(y, eps = eps))
  
  if(useYJ){
    tf <- result$transformed
    lambda <- result$opt.lambda
    whichPos <- setdiff(1:n, union(whichNegs, which0))
    tf[whichPos] <- (y[whichPos]^lambda - 1)/lambda
    tf[which0] <- log(y[which0])
    if(lambda == 2){
      tf[whichNegs] <- -log(tf[whichNegs])
    } else tf[whichNegs] <- 
      -(y[whichNegs]^(2 - lambda) -1)/(2 - lambda)
  }
  
  return(result)
}


# effect.size
# Effect sizes (standard deviates) form random outcomes
# any analytical function

#' Obtain Effect-size from a vector of values
#'
#' A function to find the effect size (Z-score) of a target,
#' from a vector of values presumably obtained in random permutations.
#'
#' @param x The vector of data to use.
#' @param center Logical value for whether to center x.
#' @param ... Arguments passed to \link{powerTrans}.
#' @export
#' @author Michael Collyer
#' @keywords utilities
effect.size <- function(x, center = TRUE, ...) {
  es.args <- c(list(y = x), list(...))
  if(length(unique(x)) == 1) {
    sdx <- 1
    x <- 0
  } else {
    bc <- do.call(powerTrans, es.args)
    x <- bc$transformed
    n <- length(x)
    if(center) x <- center(x)
    sdx <- sqrt((sum(x^2)/n))
  }
  
  (x[1]- mean(x)) / sdx
}



