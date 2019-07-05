#' Ordination tool for data aligned to another matrix
#'
#' Function performs a singular value decomposition of ordinary least squares (OLS) or 
#' generalized least squares (GLS) residuals, aligned to an alternative matrix, plus projection of 
#' data onto vectors obtained.  
#' 
#' The function performs a singular value decomosition, \bold{A'Z} = \bold{UDV'}, where
#' \bold{Z} is a matrix of residuals (obtained from \bold{Y} - see below) and \bold{A}
#' is an alignment matrix with the same number of rows as \bold{Z}.  
#' (\bold{'} indicates matrix transposition.)  \bold{U}  and \bold{V}
#' are the matrices of left and right singular vectors, and \bold{D} is a diagonal matrix of singular
#' values. \bold{V} are the vectors that describe maximized covariation between \bold{Y} and \bold{A}. 
#' If \bold{A} = \bold{I}, an n x n identity matrix, \bold{V} are the
#' eigen vectors (principal components) of \bold{Y}.
#' 
#' \bold{Z} represents a centered and potentially standardized form of \bold{Y}.  This
#' function can center data via OLS or GLS means (the latter if a covariance matrix  to describe
#' the non-independence among observations is provided).  If standardizing variables is preferred,
#' then \bold{Z} both centers and scales the vectors of \bold{Y} by their standard deviations.
#' 
#' Data are prokected onto aligned vectors, \bold{ZV}, which in the case of OLS residuals is
#' an orthogonal projection and in the case of GLS is an oblique projection.
#'
#' The versatility of using an alignment approach is that alternive data space rotations are possible.
#' Principal components are thus the vectors that maximize variance with respect to the data, themselves,
#' but "components" of (co)variation can be described for any inter-matrix relationship, including
#' phylogenetic signal, ecological signal, ontogenetic signal, size allometry, etc.
#' More details are provided in Collyer and Adams (in review).
#' 
#' Much of this function is consisten with the \code{\link{prcomp}} function, except that centering data
#' is not an option (it is required).
#' 
#' @param Y An n x p data matrix.
#' @param A An optional n x n symmetric matrix or an n x k data matrix, where k is the number of variables that could
#' be associated with the p variables of Y.  If NULL, an n x n identity matrix will be used.
#' @param A An optional n x n covariance matrix to describe the non-independence among
#' observations in Y, and provide a GLS-centering of data.  Note that Cov and A can be the same, if one
#' wishes to align GLS residuals to the same matrix used to obtain them.
#' @param scale. a logical value indicating whether the variables should be scaled to have unit variance before the analysis 
#' takes place. The default is FALSE.
#' @param tol A value indicating the magnitude below which components should be omitted. (Components are omitted if their 
#' standard deviations are less than or equal to tol times the standard deviation of the first component.) 
#' With the default null setting, no components are omitted (unless rank. is specified less than min(dim(x)).). 
#' Other settings for tol could be tol = 0 or tol = sqrt(.Machine$double.eps), which would omit essentially constant components.
#' This argument is exactly the same as in \code{\link{prcomp}}
#' @param rank. Aptionally, a number specifying the maximal rank, i.e., maximal number of aligned components to be used. 
#' This argument can be set as alternative or in addition to tol, useful notably when the desired rank is considerably 
#' smaller than the dimensions of the matrix.  This argument is exactly the same as in \code{\link{prcomp}}
#' @param newdata An optional data frame of values for the same variables of Y to be projected onto 
#' aligned components.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{ordinate} is a list containing the following
#' \item{x}{Aigned component scores for all observations}
#' \item{xn}{Optional projection of new data onto components.}
#' \item{d}{The portion of the singular values attributed to the aligned components.}
#' \item{sdev}{Standard deviations of d; i.e., the scale of the components.}
#' \item{rot}{The matrix of variable loadings, i.e. the singular vectors, \bold{V}.}
#' \item{center}{The OLS or GLS means vector used for centering.}
#' \item{scale}{The sclaing used, or FALSE.}
#' @references Collyer and Adams, in prep.
#' @seealso \code{\link{prcomp}} and \code{gm.prcomp} within \code{geomorph}
#' @examples
#' 
#' # Examples use geometric morphometric data
#' # See the package, geomorph, for details about obtaining such data
#'
#' data("PlethMorph")
#' R <- lm.rrpp(cbind (TailLength, HeadLength, Snout.eye, BodyWidth, 
#' Forelimb, Hindlimb) ~ SVL, iter = 0, 
#' data = PlethMorph, print.progress = FALSE)$LM$residuals
#' 
#' PCA.ols <- ordinate(R, scale. = TRUE)
#' PCA.ols$rot
#' prcomp(R, scale. = TRUE)$rotation # should be the same
#' 
#' PCA.gls <- ordinate(R, scale. = TRUE, Cov = PlethMorph$PhyCov)
#' 
#' # Align to phylogenetic signal
#' 
#' PaCA.ols <- ordinate(R, A = PlethMorph$PhyCov, scale. = TRUE)
#' PaCA.gls <- ordinate(R, A = PlethMorph$PhyCov, scale. = TRUE,
#' Cov = PlethMorph$PhyCov)
#' 
#' par(mfrow = c(2,2))
#' plot(PCA.ols, main = "PCA OLS")
#' plot(PCA.gls, main = "PCA GLS")
#' plot(PaCA.ols, main = "PaCA OLS")
#' plot(PaCA.gls, main = "PaCA GLS")
#' par(mfrow = c(1,1))
#' 
ordinate <- function(Y, A = NULL, Cov = NULL, scale. = FALSE, 
                     tol = NULL, rank. = NULL, newdata = NULL) {
  Y <- try(as.matrix(Y), silent = TRUE)
  if(inherits(Y, "try-error"))
    stop("Y must be a matrix or data frame.\n", call. = FALSE)

  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  I <- diag(n)
  if(is.null(A)) A <- I
  if(!is.matrix(A))
    stop("A must be a matrix with the same number of rows as Y\n", 
         call. = FALSE)
  if(NROW(A) != n)
    stop("A must be a matrix with the same number of rows as Y\n", 
         call. = FALSE)
  if(is.null(rownames(A))) rownames(A) <- rownames(Y) else {
    nnames <- length(intersect(rownames(Y), rownames(A)))
    if(nnames > n)
      stop("The row names of A are not the same as the row names of Y\n",
           call. = FALSE)
    if(isSymmetric(A)) A <- A[rownames(Y), rownames(Y)] else
      A <- A[rownames(Y),]
  }
  X <- matrix(1, n)
  rownames(X) <- rownames(Y)
  if(!is.null(Cov)) Pcov <- Cov.proj(Cov, rownames(Y))
  if(is.null(Cov)) H <- tcrossprod(X)/n else 
    H <- X %*% solve(crossprod(crossprod(Pcov, X))) %*% crossprod(crossprod(Pcov, X), Pcov)
  
  Z <- scale(as.matrix((I - H) %*% Y), center = FALSE, scale = scale.)
  sc <- attr(Z, "scaled:scale")
  if (any(sc == 0)) 
    stop("cannot rescale a constant/zero column to unit variance")
  
  k <- if (!is.null(rank.)) {
    stopifnot(length(rank.) == 1, is.finite(rank.), as.integer(rank.) > 
                0)
    min(as.integer(rank.), n, p)
  } else min(n, p)
  
  s <-  if(!is.null(Cov)) svd(crossprod(A, Pcov %*% Z), 
                              nu = 0, nv = k) else svd(crossprod(A, Z), 
                                                       nu = 0, nv = k)
   
  j <- seq_len(k)
  s$v <- s$v[,j]
  x <- Z %*% s$v
  s$d <- apply(x, 2, sd)
  if (!is.null(tol)) {
    rank <- sum(s$d > (s$d[1L] * tol))
    if (rank < k) {
      j <- seq_len(k <- rank)
      s$v <- s$v[, j, drop = FALSE]
      s$d <- s$d[j]
      x <- x[j]
    }
  }
  dimnames(s$v) <- list(colnames(Z), paste0("PC", j))
  
  r <- list(d = s$d^2, sdev = s$d, rot = s$v, 
            center = colMeans(H %*% Y), 
            scale = if (is.null(sc)) FALSE else sc,
            x = x)

  if(!is.null(newdata)){
    xn <- as.matrix(newdata)
    xn <- scale(xn, center = r$center, scale = r$scale)
    if(NCOL(Z) != NCOL(xn)) stop("Different number of variables in newdata\n", call. = FALSE) 
    r$xn <- xn %*% s$v
  }

  class(r) <- "ordinate"
  r
  
}