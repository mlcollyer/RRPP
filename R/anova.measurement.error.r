#' ANOVA for lm.rrpp model fits used in measurement.error
#'
#' @description Computes an analysis of variance (ANOVA) table using 
#' distributions of random statistics from \code{\link{lm.rrpp}}.  
#' This function is the same as \code{\link{anova.lm.rrpp}} but includes statistics
#' specific to \code{\link{measurement.error}}, and with restrictions on
#' how P-values and effect.sizes are calculated.
#'
#' @param object Object from \code{\link{lm.rrpp}}
#' @param ... Additional lm.rrpp model fits or other arguments passed to anova. 
#' @export anova.measurement.error
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' # See measurement.error help file examples for use.
#' 
anova.measurement.error <- function(object, ...) {
  out <- aov.me(object)
  out
}
  



