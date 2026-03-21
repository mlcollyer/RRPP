#' Likelihood ratio test for a linear model, based on RRPP
#' 
#' @description 
#' Function performs likelihood ratio tests on an lm.rrpp fit, using 
#' RRPP or FRPP.  Likelihood ratio statistics are calculated for every random 
#' permutation, and the effect size is estimated from the distribution of 
#' random statistics.  The likelihood ratio tests has some resemblance to
#' MANOVA, especially using Wilks' lambda.  Sums of squares and cross-products 
#' (SSCP) matrices are calculated over the random permutations of 
#' a \code{\link{lm.rrpp}} fit.  SSCP matrices are 
#' computed, as are the inverse of R times H (invR.H), where R is a SSCP 
#' for the residuals or random effects and H is
#' the difference between SSCP matrices of full and reduced models 
#' (see \code{\link{manova.update}}).   From invR.H, Wilks lambda is first estimated, and the 
#' likelihood ratio stat is then estimated as -n * log(Wilks).
#' 
#' This function does one of two things.  It either performs an update using
#' \code{\link{manova.update}}, using Wilks' lambda as the test statistic, converting 
#' Wilks' lambda to likelihood ratio statistics or it uses the results from
#' a previously performed update to calculate new statistics.  
#' 
#' 
#' @references Adams, D. C., and M. L. Collyer. 2024. Extended phylogenetic 
#' regression models for comparing within-species patterns across the 
#' tree of life. Methods in Ecology and Evolution. In review.
#' @param fit Linear model fit from \code{\link{lm.rrpp}} or a fit
#' that has already been updated with \code{\link{manova.update}}.
#' @param verbose Logical value for whether to include all random 
#' Wilks' lambda and likelihood ratio statistics from random permutations.
#' @param ... Arguments passed onto \code{\link{manova.update}}.

#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @examples 
#'    
#' # Body Shape Analysis (Multivariate) ----------------
#' 
#' \dontrun{
#' data(Pupfish)
#' 
#' # Although not recommended as a practice, this example will use only
#' # three principal components of body shape for demonstration.  
#' # A larger number of random permutations should also be used.
#' 
#' Pupfish$shape <- ordinate(Pupfish$coords)$x[, 1:3]
#' 
#' 
#' fit <- lm.rrpp(shape ~ log(CS) + Sex, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 499) 
#' summary(fit, formula = FALSE)
#' anova(fit) # ANOVA table
#' 
#' # MANOVA
#' 
#' fit.m <- manova.update(fit, print.progress = FALSE, tol = 0.001)
#' summary(fit.m, test = "Roy")
#' summary(fit.m, test = "Wilks")
#' 
#' 
#' # Likelihood Ratio Test
#' 
#' LRT <- lr_test(fit.m)
#' summary(LRT)
#' 
#' }


lr_test <- function(fit, verbose = FALSE, ...){
  
  manova.args <- list(...)
  need.invR.H <- is.null(fit$MANOVA$invR.H)
  n <- fit$LM$n
  
  if(need.invR.H) {
    if(inherits(fit, "manova.lm.rrpp")) {
      man.loc <- which(class(fit) == "manova.lm.rrpp")
      class(fit) <- class(fit)[-man.loc]
    }
    
    manova.args$fit <- fit
    fit <- NULL
    manova.args$error <- NULL
    if(is.null(manova.args$tol)) manova.args$tol <- 1e-07
    if(is.null(manova.args$print.progess)) manova.args$print.progress <- TRUE
    manova.args$verbose <- TRUE
    fit <- do.call(manova.update, manova.args)
  }
  
  tol <- min(c(1e-07, manova.args$tol))
  eigs <- fit$MANOVA$eigs
  fit <- NULL
  
  Wilks <- sapply(eigs, function(x){
    apply(x, 1, wilks)
  })
  
  eigs <- NULL
  LR <- -n * log(Wilks + tol)
  
  tab <- data.frame(
    Wilks  = Wilks[,1],
    LR = LR[,1],
    Z = apply(LR, 1, effect.size),
    pval = apply(LR, 1, pval)
  )
  
  out <- list(tab = tab)
  
  if(verbose){
    out$wilks <- Wilks
    out$lr <- LR
  }
  
  class(out) <- "lr_test"
  out
    
}
