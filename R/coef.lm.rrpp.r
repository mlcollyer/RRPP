#' coef for lm.rrpp model fits
#'
#' @description Computes ordinary or generalized least squares coefficients
#' over the permutations of an \code{\link{lm.rrpp}} model fit with predefined random permutations.
#' For each coefficient vector, the Euclidean distance is calculated as an estimate of
#' the amount of change in Y, the n x p matrix of dependent variables; larger distances mean more change 
#' in location in the data space associated with a one unit change in the model design, for the parameter
#' described.  Random coefficients are based on either RRPP or FRPP, as defined by the 
#' \code{\link{lm.rrpp}} model fit.  If RRPP is used, all distributions of coefficient vector distances are 
#' based on appropriate null models as defined by SS type.
#'
#' @param object Object from \code{\link{lm.rrpp}}
#' @param confidence The desired confidence interval level to print with a table of summary statistics
#' @param ... Other arguments (currently none)
#' @export
#' @author Michael Collyer
#' @keywords utilities
coef.lm.rrpp <- function(object, confidence = 0.95, ...) {
  x <- beta.lm.rrpp(object)
  coef.obs = x$coef.obs
  rc = x$random.coef
  rd = x$random.distances
  n = x$n; p = x$p; k = x$k.terms
  model.terms = x$model.terms
  perms = x$nperms
  SS.type <- x$SS.type
  RRPP <- x$RRPP
  gls <- x$gls
  PV <- apply(rd, 1, pval)
  Z <- apply(rd, 1, effect.size)
  alpha = 1 - confidence
  if(alpha < 0) stop("Confidence level should be between 0 and 1")
  lcl <- apply(rd, 1, function(x) quantile(x, alpha/2))
  ucl <- apply(rd, 1, function(x) quantile(x, (1 - alpha/2)))
  stat.tab <- data.frame(d.obs = rd[,1], lcl=lcl, ucl=ucl, Z=Z, P=PV)
  colnames(stat.tab)[2] <- paste("LCL (", alpha/2*100,"%)", sep = "")
  colnames(stat.tab)[3] <- paste("UCL (", (1 -alpha/2)*100,"%)", sep = "")
  colnames(stat.tab)[5] <- "Pr(>d)"
  out <- list(coef.obs = coef.obs,
              random.coef = rc,
              random.distances = rd,
              n = n, p=p, k.terms = k, confidence = confidence,
              model.terms = model.terms, nperms = perms,
              RRPP = RRPP, gls=gls, SS.type = SS.type,
              stat.table = stat.tab)
  class(out) <- "coef.lm.rrpp"
  out
}