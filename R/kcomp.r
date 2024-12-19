#' K-component analysis
#'
#' Function performs a K-comonent analysis, which is an eigen decomposition
#' of a ratio of ordinary least-squares (OLS) to generalized least squares (GLS)
#' covariance matrices.
#' 
#' The function performs a K-component analysis (KCA), as described in Mitteroecker et al. (2025).  
#' The analysis is similar to many that can be performed with  the \code{\link{ordinate}} function,
#' but whereas that function finds ordinate scores for a rotation or shear of a data space, KCA performs 
#' a relative eigen decomposition of a matrix product that represents a ratio between two covariance
#' matrices (Mitteroecker and Bookstein, 2014).  The eigenvalues are helpful for understanding the importance
#' of the GLS estimation of a covariance matrix.  For example, if the GLS estimation is made with
#' respect to a matrix of phylogenetic covariances, the eigenvalues express the distribution of 
#' phylogenetic signal across components of the data space.  
#' 
#' The matrix product (K matrix) is the inverse of the OLS covariance matrix times 
#' the GLS covariance matrix.  
#' Redundancies between these matrices will mean the K matrix is likely not full rank.  Therefore,
#' an algorithm is used to compute scores in appropriate dimensions.  First, data are aligned to the 
#' square-root matrix of the GLS
#' covariance matrix (Collyer and Adams, 2021), using the \code{\link{ordinate}} function.  Second, relative eigen decomposition is 
#' performed, and the number of real eigenvectors that can be computed are retained.  Lastly, the same number of 
#' aligned components are used for projection of data onto relative eigenvectors (the K components).
#' 
#' 
#' @param Y An n x p data matrix.
#' @param Cov An required n x n covariance matrix to describe 
#' the non-independence among
#' observations in Y, and provide a GLS-centering of data. 
#' @param transform. An optional argument to transform GLS-centered residuals, 
#' if TRUE.  If FALSE, only GLS-centering is performed.  
#' @param scale. A logical value indicating whether the variables 
#' should be scaled to have unit variance before the analysis 
#' takes place. The default is FALSE.
#' @param tol A value indicating the magnitude below which 
#' components should be omitted. (Components are omitted if their 
#' standard deviations are less than or equal to tol times the 
#' standard deviation of the first component.) 
#' @param rank. Optionally, a number specifying the maximal rank, 
#' i.e., maximal number of K components to be used. 
#' This argument can be set as alternative or in addition to tol, 
#' useful notably when the desired rank is considerably 
#' smaller than the dimensions of the matrix of the K matrix.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{kcomp} is a list containing 
#' the following
#' \item{values}{The eigenvalues of the K matrix,}
#' \item{vectors}{The eigenvectors of the K matrix.}
#' \item{scores}{The projected scores of data onto the eigenvectors.}

#' @references Collyer, M.L. and D.C. Adams.  2021. 
#' Phylogenetically-aligned Component Analysis. Methods 
#' in Ecology and evolution. In press.
#' @references Bookstein, F. L., & Mitteroecker, P. (2014). Comparing 
#' covariance matrices by relative eigenanalysis, with applications to 
#' organismal biology. Evolutionary biology, 41, 336-350.
#' @references Mitteroecker, P., Collyer, M. L., & Adams, D. C. (2024). 
#' Exploring Phylogenetic Signal in Multivariate Phenotypes by Maximizing 
#' Blomberg’s K. Systematic Biology, syae035.
#' 
#' @seealso \code{\link{plot.kcomp}}, \code{\link{ordinate}}, 
#' \code{\link{plot.default}}, \code{\link{rrpp.data.frame}}
#' @examples
#'  # TBD
#'  
kcomp <- function(Y, Cov, 
                  transform. = TRUE,
                     scale. = FALSE, 
                     tol = NULL, rank. = NULL) {
  
  if(is.null(rownames(Y)))
    stop("Y must contain row names that are the same as in Cov.",
         call. = FALSE)
  
  Pcov <- Cov.proj(Cov, id = rownames(Y), 
                   symmetric = TRUE)
  Covsqrt <- fast.solve(Pcov)
  
  PACA <- ordinate(Y, Cov = Covsqrt, transform. = transform., 
                   scale. = scale., tol = tol, rank. = rank.)
  
  z <- PACA$x
  
  SSCPo <- crossprod(z)
  SSCPg <- crossprod(Pcov %*% z)
  
  K <- fast.solve(SSCPg) %*% SSCPo
  eig <- eigen(K)
  scores <- z %*% eig$vectors
  out <- list(values = eig$values, 
              vectors = eig$vectors, 
              scores = scores)
  
  class(out) <- "kcomp"
  out
  
}
