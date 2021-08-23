#' coef for lm.rrpp model fits
#'
#' @description Computes ordinary or generalized least squares coefficients
#' over the permutations of an \code{\link{lm.rrpp}} model fit with predefined 
#' random permutations.
#' For each coefficient vector, the Euclidean distance is calculated as an 
#' estimate of
#' the amount of change in Y, the n x p matrix of dependent variables; larger 
#' distances mean more change 
#' in location in the data space associated with a one unit change in the model 
#' design, for the parameter
#' described.  Random coefficients are based on either RRPP or FRPP, as defined 
#' by the 
#' \code{\link{lm.rrpp}} model fit.  If RRPP is used, all distributions of 
#' coefficient vector distances are 
#' based on appropriate null models as defined by SS type.
#' 
#' This function can be used to test the specific coefficients of an 
#' lm.rrpp fit.  The test
#' statistics are the distances (d), which are also standardized (Z-scores).  
#' The Z-scores might be easier to compare,
#' as the expected values for random distances can vary among coefficient 
#' vectors (Adams and Collyer 2016).
#'
#' @param object Object from \code{\link{lm.rrpp}}
#' @param test Logical argument that if TRUE, performs hypothesis tests 
#' (Null hypothesis is vector distance = 0)
#' for the observed coefficients.  If FALSE, only the observed coefficients 
#' are returned.
#' @param confidence The desired confidence limit to print with a table of 
#' summary statistics,
#' if test = TRUE.  Because distances are directionless, confidence limits 
#' are one-tailed.
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
#' Pupfish$logSize <- log(Pupfish$CS)
#'
#' fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", data = Pupfish) 
#' 
#' coef(fit)
#' coef(fit, test = TRUE, confidence = 0.99)
coef.lm.rrpp <- function(object, test = FALSE, confidence = 0.95, ...) {
  x <- object
  k <- length(x$LM$term.labels)
  rc <- x$LM$random.coef
  rd <- x$LM$random.coef.distances
  coef.obs <- if(x$LM$gls) x$LM$gls.coefficients else x$LM$coefficients
  n <- x$LM$n; p <- x$LM$p; p.prime = x$LM$p.prime
  model.terms <- x$LM$Terms
  
  perms <- x$PermInfo$perms
  SS.type <- x$ANOVA$SS.type
  RRPP <- x$PermInfo$perm.method
  gls <- x$LM$gls
  
  if(test && x$turbo) {
    cat("\nCoefficients test not available because you turbo-charged your model fit.\n")
    cat("Go back to lm.rrpp and choose turbo = FALSE ")
    cat("if you wish to also test coefficients.\n\n")
    test = FALSE
  }
  if(test){
    if(confidence < 0) stop("Confidence level should be between 0 and 1")
    if(confidence > 1) stop("Confidence level should be between 0 and 1")
    if(k > 0) {
      PV <- apply(rd, 1, pval)
      Z <- apply(rd, 1, effect.size)
      ucl <- apply(rd, 1, function(x) quantile(x, confidence))
      stat.tab <- data.frame(d.obs = rd[,1], ucl=ucl, Zd=Z, P=PV)
      colnames(stat.tab)[2] <- paste("UCL (", confidence*100,"%)", sep = "")
      colnames(stat.tab)[4] <- "Pr(>d)"
    } else {
      PV <- Z <- ucl <- stat.tab <- NULL
    }
    
    out <- list(coef.obs = coef.obs,
                random.coef = rc,
                random.distances = rd,
                n = n, p=p, p.prime=p.prime, k.terms = k, 
                confidence = confidence,
                model.terms = model.terms, nperms = perms,
                RRPP = RRPP, gls=gls, SS.type = SS.type,
                stat.table = stat.tab, test = test)
    class(out) <- "coef.lm.rrpp"
  } else out <- coef.obs
  
  out
}