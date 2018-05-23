#' Model Comparisons, in terms of the log-likelihood or covariance trace
#'
#' Function calculates either log-likelihoods or traces of covariance matrices for comparison with 
#' respect to parameter penalties.
#' 
#'
#' The function calculates either log-likelihoods or traces of (residual) covariance matrices, plus parameter 
#' penalties, to assist in comparative model evaluation or selection.  Because high-dimensional data often
#' produce singular or ill-conditioned residual covariance matrices, this function does one of two things: 1) uses 
#' the trace of a covariance matrix rather than its determinant; or 2) provides a ridge-regularization (Warton, 2008)
#' of the covariance matrix, only if it is determined that it is ill-conditioned.  Regardless of implementation,
#' covariance matrices are projected onto a principal component manifold of appropriate dimensions.
#' 
#' The parameter penalty is based on that proposed by Bedrick and Tsai (1994), equal to 2(pk + p(p + 1)/2), where 
#' p is the appropriate dimension (not number of variables) of the covariance matrix.  The parameter, k,
#'  is the rank of the model design matrix.
#'  
#' In the case that "logLik" is chosen for the argument, type, AIC scores are calculated.  These scores
#' may not perfectly match other packages or software that calculate AIC for multivariate data, if ridge regualrization
#' was used (and if other packages let p = the number of data variables).  Users can construct their own tables 
#' from the results but this function does not attempt to summarize results, as interpreting results requires 
#' some arbitrary decisions.  The \code{\link{anova}} function explicitly tests multiple models and can be used for nested 
#' model comparisons.
#' 
#' Results can also be plotted using the generic \code{\link{plot}} function.
#' 
#' Caution: For models with GLS estimation, determinants of large Kronecker products are required for 
#' log-likelihood.  This requires reducing the data space dimensions via projection onto
#' principal components (about 50% reduction).  It might be better to rely on the trace statistic.  (- Infinite log-likelihoods
#' are possible and no further attempt to avoid them is provided.)
#' 
#' @param ... Any number of lm.rrpp class objects for model fits to be compared.
#' @param type An argument to choose between log-likelihood or covariance trace results
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{model.comparison} is a data frame with either log-likelihoods
#' or covariance traces, plus parameter penalties.  AIC scores might be include, if applicable
#' @references Bedrick, E.J., and C.L. Tsai. 1994. Model selection for multivariate regression in small samples. 
#' Biometrics, 226-231.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric
#' datasets. Evolution. 70:2623-2631.
#' @examples Warton, D.I., 2008. Penalized normal likelihood and ridge regularization of correlation and covariance matrices. 
#' Journal of the American Statistical Association. 103: 340-349.
#' 
#' data(Pupfish)
#' Pupfish$logSize <- log(Pupfish$CS)
#' fit1 <- lm.rrpp(coords ~ logSize, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit2 <- lm.rrpp(coords ~ Pop, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit3 <- lm.rrpp(coords ~ Sex, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit4 <- lm.rrpp(coords ~ logSize + Sex, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit5 <- lm.rrpp(coords ~ logSize + Pop, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit6 <- lm.rrpp(coords ~ logSize + Sex * Pop, data = Pupfish, iter = 0, print.progress = FALSE)
#' 
#' modComp1 <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6, type = "cov.trace")
#' modComp2 <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6, type = "logLik")
#' 
#' summary(modComp1)
#' summary(modComp2)
#' 
#' par(mfcol = c(1,2))
#' plot(modComp1)
#' plot(modComp2)
#' 
model.comparison<- function(..., type = c("cov.trace", "logLik")) {
  
  dots <- list(...)
  check <- lapply(dots, class)
  if(any(check != "lm.rrpp")) stop("\nObjects must be lm.rrpp fits\n.")
  dot.names <- lapply(dots, function(x) x$LM$Terms[[3]])
  type = match.arg(type)
  if(type == "logLik") res <- sapply(dots, logL) else
    res <- sapply(dots, cov.trace) 
  
  par.pen <- function(f){
    p <- f$LM$p.prime
    k <- f$LM$QR$rank
    2 * (p * k + 0.5 * p* (p + 1))
  }
  pp <- sapply(dots, par.pen)
  
  if(type == "logLik") {
    out <- cbind(res, pp, -2*res + pp)
    colnames(out) <- c("logLik", "penalty", "AIC")
  } else {
    out <- cbind(res, pp)
    colnames(out) <- c("cov.trace", "penalty")
  }
  
  out <- list(table = out, names = dot.names)
  attr(out, "class") = "model.comparison"
  out
  
}