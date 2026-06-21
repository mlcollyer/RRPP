
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
  
  y <- scale(y)
  y <- y - min(y) + 1
  y.obs <- y[1]
  y <- y[-1]
  
  # Note the code below (to get yy) is short-cut for Jacobian adjustment 
  n <- length(y)
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
  
  y <- scale(y)
  y <- y - min(y) + 1
  y.obs <- y[1]
  y <- y[-1]
  
  n <- length(y)
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
  
  y <- scale(y)
  y <- y - min(y) + 1
  y.obs <- y[1]
  y <- y[-1]
  
  n <- length(y)
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

#' Box-Cox Power Transformation
#'
#' A function that transforms a vector in data to be 
#' approximately normally distributed.
#' 
#' This function uses a Box-Cox power transformation for a vector of values,
#' (presumably statistics from RRPP) with the following conditions:
#' 
#' (1) Values are first transformed as z - min(z) + 1, where z 
#' is a vector of standard deviates of the original values.  
#' This assures values are positive
#' with a minimum value of 1, such that a lambda of 0 returns a 
#' minimum transformed value of 0.  This also allows transformation of
#' negative values.
#' 
#' 
#' (2) The lambda parameter is optimized over a range of -2 to +2,
#' as is typical for most applications.
#' 
#' (3) The first value (assumed to be the observed value in a 
#' distribution of random values) is removed from the vector to 
#' perform optimization of lambda, but is returned to the vector for
#' transformation.  This assures that large, aberrant values, such
#' as might be found via RRPP for large effects, do not have 
#' leverage for the optimization.
#' 
#' This function is similar to \link[MASS]{boxcox} in most regards
#' but uses a consistent initial scaling and shifting of values for comparative
#' purposes, when calculating effect sizes with \link{effect.size}.  In essence,
#' it is the same as using the \link[MASS]{boxcox} function, after
#' first scaling values and adding 1.  However, it also assumes
#' the first value is an observed value in a vector of random statistics.  It
#' is intended for RRPP functions, only, but could be used with other results.  
#' However, one must be aware that the first value will be removed for optimization.
#'
#' @param y A vector of values to use.
#' @param eps Tolerance for lambda = 0.
#' @param interp Logical value for whether to spline interpolate lambda 
#' (ideal for vectors of small length).  This will default if the number
#' of values is less than 100, as in \link[MASS]{boxcox}.
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

box.cox <- function(y, eps = 0.001, interp = FALSE) {
  if(length(y) < 100) interp <- TRUE
  result <- if(interp) box.cox.interp(y, eps = eps) else
    box.cox.fast(y, eps = eps)
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
#' @param target The value to target in the distribution.  (If null, the first value
#' in the vector is used.).  If the target exists outside the range of x,
#' very small or very large z-scores are possible.  Additionally, if the target
#' is excessively outside of the range of x, it could affect the Box-Cox transformation 
#' used to transform x.
#' @export
#' @author Michael Collyer
#' @keywords utilities
effect.size <- function(x, center = TRUE, target = NULL) {
  if(is.null(target)) target <- x[1] else x <- c(target, x)
  if(length(unique(x)) == 1) {
    sdx <- 1
    x <- 0
  } else {
    x <- box.cox(x)$transformed
    n <- length(x)
    if(center) x <- center(x)
    sdx <- sqrt((sum(x^2)/n))
  }
  
  (x[1]- mean(x)) / sdx
}
