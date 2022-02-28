#' Reveal model designs used in lm.rrpp fit
#'
#' Function returns every full and reduced model for model terms used in 
#' lm.rrpp fits.  This function is useful for revealing 
#' the null and full model that would be used in the pairwise function, 
#' if a specific null model is not declared as an argument
#' (fit.null in the \code{\link{pairwise}} function).
#' It also helps to demonstrate how sums of squares and cross-products 
#' (SSCP) are calculated in lm.rrpp permutations (iterations),
#' from the difference between fitted values for null and full designs.
#' 
#' @param fit A linear model fit from \code{\link{lm.rrpp}}.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @examples
#'
#'data(Pupfish)
#'fit1 <- lm.rrpp(coords~ Pop*Sex, data = Pupfish, 
#'SS.type = "I", print.progress = FALSE, iter = 0)
#'fit2 <- lm.rrpp(coords~ Pop*Sex, data = Pupfish, 
#'SS.type = "II", print.progress = FALSE, iter = 0)
#'fit3 <- lm.rrpp(coords~ Pop*Sex, data = Pupfish, 
#'SS.type = "III", print.progress = FALSE, iter = 0)
#'
#'reveal.model.designs(fit1)
#'reveal.model.designs(fit2)
#'reveal.model.designs(fit3)
#' 
reveal.model.designs <- function(fit) {
  model.sets <- fit$Models
  terms.f <- lapply(model.sets$full, function(x) x$terms)
  terms.r <- lapply(model.sets$reduced, function(x) x$terms)
  forms.r <- lapply(terms.r, function(x) formula(x)[[3]])
  forms.f <- lapply(terms.f, function(x) formula(x)[[3]])
  reduced <- lapply(forms.r, function(x) Reduce(paste, deparse(x)))
  full<- lapply(forms.f, function(x) Reduce(paste, deparse(x)))
  term.labels <- names(forms.f)
  k <- length(term.labels)
  blank <- target <- rep("", k)
  target[k] <- "<- Null/Full inherent in pairwise"
  df <- as.data.frame(cbind(blank, reduced = reduced, blank, full = full, target))
  names(df) <- c("", "Reduced", "", "Full", "")
  return(df)
}