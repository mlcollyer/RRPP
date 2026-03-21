#' Calculate the pairwise Mahalanobis distances between observations
#'
#' This function emulates the \code{\link[stats]{dist}} function
#' but allows a covariance matrix (Cov) to be included for standardizing
#' distances.  It is assumed that the Covariance matrix makes sense with
#' respect to the data, and that the number of variables match between data and covariance matrix.
#' 
#' No tests are performed on distances but could be performed with the 
#' \code{\link{pairwise}} function.  Distances are only calculated if
#' the covariance matrix is not singular.
#'
#' @param x A numeric matrix of data frame.
#' @param Cov A covariance matrix with the same number of variables as the data.
#' @param ... Other arguments passed to \code{\link[stats]{dist}}.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class "dist".
#' 
#' @examples 
#' 
#' # Using the Pupfish data (see lm.rrpp help for more detail)
#' 
#' data(Pupfish)
#' Pupfish$Y <- ordinate(Pupfish$coords)$x[, 1:3]
#' fit <- lm.rrpp(Y ~ Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 0)
#' X <- unique(as.matrix(model.matrix(fit)))
#' means <- X %*% coef(fit)
#' rownames(means) <- unique(interaction(Pupfish$Sex, Pupfish$Pop))
#' means
#' S <- getResCov(fit)
#' dist(means)
#' mahal_dist(means, S)

mahal_dist <- function(x, Cov, ...){
  x <- as.matrix(x)
  if(NCOL(x) != NCOL(Cov))
    stop("\nThe number of data variables do not match the number of covariance matrix variables.\n",
         call. = FALSE)
  Pcov <- Cov.proj(Cov)
  TX <- tcrossprod(x, Pcov)
  dist(TX, ...)
}