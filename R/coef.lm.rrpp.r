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
#' @examples 
#' # See examples for lm.rrpp to see how anova.lm.rrpp works in conjunction
#' # with other functions
#' 
#' data(Pupfish)
#' names(Pupfish)
#' Pupfish$logSize <- log(Pupfish$CS) # better to not have functions in formulas
#'
#' fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", data = Pupfish) 
#' 
#' coef(fit)
#' coef(fit, confidence = 0.99)
coef.lm.rrpp <- function(object, confidence = 0.95, ...) {
  x <- object
  rc <- x$LM$random.coef
  rd <- x$LM$random.coef.distances
  n <- x$LM$n; p <- x$LM$p
  model.terms <- x$LM$Terms
  k <- length(x$LM$term.labels)
  coef.obs <- rc[[k]][[1]]
  perms <- x$PermInfo$perms
  SS.type <- x$ANOVA$SS.type
  RRPP <- x$PermInfo$perm.method
  gls <- x$L$gls
  alpha = 1 - confidence
  if(alpha < 0) stop("Confidence level should be between 0 and 1")
  if(k > 0) {
    PV <- apply(rd, 1, pval)
    Z <- apply(rd, 1, effect.size)
    lcl <- apply(rd, 1, function(x) quantile(x, alpha/2))
    ucl <- apply(rd, 1, function(x) quantile(x, (1 - alpha/2)))
    stat.tab <- data.frame(d.obs = rd[,1], lcl=lcl, ucl=ucl, Z=Z, P=PV)
  } else {
    PV <- pval(rd)
    Z <- pval(rd)
    lcl <- quantile(rd, alpha/2)
    ucl <- quantile(rd, (1 - alpha/2))
    stat.tab <- data.frame(d.obs = rd[1], lcl=lcl, ucl=ucl, Z=Z, P=PV)
  }
  
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