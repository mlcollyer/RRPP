#' Diagnostic cross-validation tool for ordination based on fitted values
#'
#' Function performs a leave-one-out cross-validation estimate of ordination
#' scores, which is helpful for determining if apparent "group differences"
#' in ordination plots arise merely from data dimensionality.
#' 
#' The function uses the strategy of Thioulouse et al. (2021) to perform N
#' ordinations for N observations, in which each of the N observations are left
#' out of the estimation of linear model coefficients, but the vector of data for
#' the left-out observation is projected on the eigenvectors of the fitted values 
#' obtained from the leave-one-out cross-validation (jackknife) strategy.
#' The purpose of this diagnostic tool is to determine whether apparent "group differences"
#' in an ordination plot (using the function, \code{\link{ordinate}}) are
#' because of high-dimensional data (number of variables exceed number of observations)
#' rather than real differences.  An apparent group difference is common for high-dimensional
#' data, when variables are far greater in number than observations (Cardini et al., 2019).
#' However, leave-one-out cross-validation can help elucidate whether an observed visual 
#' difference is spurious.
#' 
#' This function differs from the strategy of Thioulouse et al. (2021) in two important 
#' ways.  First, this function uses the linear model design from a \code{\link{lm.rrpp}}
#' fit, and can contain any number of independent variables, rather than a single factor 
#' for groups.  Second, after obtaining leave-one-out cross-validated scores, a Procrustes
#'  alignment between cross-validated scores and "observed" (real) scores is performed, which 
#' minimizes summed squared distances between the alternative ordinations.  This latter
#' step assures comparisons are appropriate.
#' 
#' The type = "PC" plot from \code{\link{plot.lm.rrpp}} has the same scores as obtained
#' from ordinate(Y, A = H), using the \code{\link{ordinate}} function, where H is a hat
#' matrix (that can be calculated from \code{\link{plot.lm.rrpp}} output), and Y is a matrix 
#' of data.  This function updates H for every possible case that one row of Y is left out
#' (meaning the rotation matrix from \code{\link{ordinate}} is updated N times).  If
#' the H matrix is robust in spite of dropped data and design matrix parameters, the result
#' will be similar to the original ordination.  If apparent group differences are spurious,
#' H will tend to change, as will data projections.
#' 
#' The functions \code{\link{summary.looCV}} and \code{\link{plot.looCV}} are essential for 
#' evaluating results.  These support functions compare eigenvalues and 
#' projected scores, between observed and cross-validated cases.  
#' 
#' This function should be viewed as a diagnostic tool and not as a data transformation tool!
#' The cross-validated scores will not retain Euclidean distances among observations.  This
#' could cause problems in analyses that substitute cross-validated scores as data.
#' 
#' 
#' @param fit A \code{\link{lm.rrpp}} fit.
#' @param ... Arguments passed to \code{\link{ordinate}}
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{looCV} is a list containing 
#' the following
#' \item{x}{Aligned component scores for all observations}

#' @references Thioulouse, J., Renaud, S., Dufour, A. B., & Dray, S. (2021). 
#' Overcoming the Spurious Groups Problem in Between-Group PCA. 
#' Evolutionary Biology, In press.
#' @references Cardini, A., O’Higgins, P., & Rohlf, F. J. (2019). Seeing distinct groups 
#' where there are none: spurious patterns from between-group PCA. 
#' Evolutionary Biology, 46(4), 303-316.
#' 
#' @seealso \code{\link{summary.looCV}}, \code{\link{plot.looCV}}
#' @examples
#' 
#' # Example with real group differences
#' 
#' data(Pupfish)
#' fit <- lm.rrpp(coords ~ Pop*Sex, data = Pupfish, iter = 0)
#' CV1 <- looCV(fit)
#' summary(CV1)
#' group <- interaction(Pupfish$Pop, Pupfish$Sex)
#' plot(CV1, flip = 1, pch = 19, col = group)
#' 
#' # Example with apparent but not real group differences
#' 
#' n <- NROW(Pupfish$coords)
#' p <- NCOL(Pupfish$coords)
#' set.seed(1001)
#' Yr <- matrix(rnorm(n * p), n, p) # random noise
#' 
#' fit2 <-lm.rrpp(Yr ~ Pop*Sex, data = Pupfish, iter = 0)
#' CV2 <- looCV(fit2)
#' summary(CV2)
#' group <- interaction(Pupfish$Pop, Pupfish$Sex)
#' plot(CV2, pch = 19, col = group) 
#' 
looCV <- function(fit, ...){
  res <- looPCAll(fit, ...)
  d <- list(obs = res$raw$d, cv = res$cv$d)
  scores <- list(obs = res$raw$x, cv = res$cv$x)
  out <- list(d = d, scores = scores)
  class(out) <- "looCV"
  out
}

