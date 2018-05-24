
#' Model Comparison tests, using z-tests
#'
#' Function find the pairwise z-scores between model fits, following ANOVA, as a means
#' to test differences between model deviances, given the same null model
#' 
#' The protocol for this function is as follows:
#' 1. Perform a multi-model ANOVA, with a specific null model (generally containing
#' only an intercept, but an alternative null model can be used).
#' 2. Obtain the pairwise z-scores following the method of Adams and Collyer (2016).
#' 3. Calculate P-values for the pairwise z-scores from the sampling distributions of
#' pairwise differences in statistics (e.g., SS), based on RRPP.
#' 4. Summarize the pariwise model comparisons.
#' 5. Perform a principal component analysis (multi-dimensional scaling) of the pairwise 
#' Z-score matrix.
#' 
#' The statistics used for the test is the one chosen for effect type in ANOVA.  The SS
#' statistics is recommended, as the z-score would most resemble the deviance between
#' models.
#' 
#' This function is experimental, and is currently the source for conceptual/theoretical 
#' work and should not be cited or used, haphazardly until confirmation of its appropriateness.  
#' 
#' @param A Object of class \code{\link{anova.lm.rrpp}}, performed on multiple models.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{z.model.test} with a list of Z-scores, P-values,
#' model names, and PCoA results, using the function \code{\link{cmdscale}}
#' @examples 
#' 
#' data(Pupfish)
#' Pupfish$logSize <- log(Pupfish$CS)
#' fit0 <- lm.rrpp(coords ~ 1, data = Pupfish, iter = 499, print.progress = FALSE)
#' fit1 <- lm.rrpp(coords ~ logSize, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit2 <- lm.rrpp(coords ~ Pop, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit3 <- lm.rrpp(coords ~ Sex, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit4 <- lm.rrpp(coords ~ logSize + Sex, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit5 <- lm.rrpp(coords ~ logSize + Pop, data = Pupfish, iter = 0, print.progress = FALSE)
#' fit6 <- lm.rrpp(coords ~ logSize + Sex * Pop, data = Pupfish, iter = 0, print.progress = FALSE)
#' 
#' aovPup <- anova(fit0, fit1, fit2, fit3, fit4, fit5, fit6, 
#' effect,type = "SS", print.progress = FALSE)
#' 
#' summary(aovPup)
#' 
#' zPup <- z.model.test(aovPup)
#' 
#' summary(zPup, names = "model")
#' plot(zPup, label = "form", pch=19)
#' 
z.model.test <- function(A){
  if(is.null(A$F)) stop("\nCannot perform a z-test on only one model.\n")
  out <- z.test(A)
  out$anova <- A
  class(out) <- "z.model.test"
  out
}