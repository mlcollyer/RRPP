#' Ordination tool for data aligned to another matrix
#'
#' Function performs a singular value decomposition of ordinary 
#' least squares (OLS) or 
#' generalized least squares (GLS) residuals, aligned to an alternative 
#' matrix, plus projection of 
#' data onto vectors obtained.  
#' 
#' The function performs a singular value decomposition, 
#' \bold{t(A)Z} = \bold{UDt(V)}, where
#' \bold{Z} is a matrix of residuals (obtained from \bold{Y}; 
#' see below) and \bold{A}
#' is an alignment matrix with the same number of rows as \bold{Z}.  
#' (\bold{t} indicates matrix transposition.)  \bold{U}  and \bold{V}
#' are the matrices of left and right singular vectors, and \bold{D} 
#' is a diagonal matrix of singular
#' values. \bold{V} are the vectors that describe maximized 
#' covariation between \bold{Y} and \bold{A}. 
#' If \bold{A} = \bold{I}, an n x n identity matrix, \bold{V} are the
#' eigen vectors (principal components) of \bold{Y}.
#' 
#' \bold{Z} represents a centered and potentially standardized 
#' form of \bold{Y}.  This
#' function can center data via OLS or GLS means (the latter if 
#' a covariance matrix to describe
#' the non-independence among observations is provided).  If 
#' standardizing variables is preferred,
#' then \bold{Z} both centers and scales the vectors of \bold{Y} 
#' by their standard deviations.  
#' 
#' Data are projected onto aligned vectors, \bold{ZV}.  If a 
#' GLS computation is made, the option to transform centered values 
#' (residuals) before projection is 
#' available.  This is required for orthogonal projection, but 
#' from a transformed data space.  Not transforming
#' residuals maintains the Euclidean distances among observations 
#' and the OLS multivariate variance, but
#' the projection is oblique (scores can be correlated).
#'
#' The versatility of using an alignment approach is that 
#' alternative data space rotations are possible.
#' Principal components are thus the vectors that maximize 
#' variance with respect to the data, themselves,
#' but "components" of (co)variation can be described for any 
#' inter-matrix relationship, including
#' phylogenetic signal, ecological signal, ontogenetic signal, 
#' size allometry, etc.
#' More details are provided in Collyer and Adams (2021).
#' 
#' Much of this function is consistent with the \code{\link{prcomp}} 
#' function, except that centering data
#' is not an option (it is required).
#' 
#' SUMMARY STATISTICS: For principal component plots, the traditional 
#' statistics to summarize the analysis include
#' eigenvalues (variance by component), proportion of variance by 
#' component, and cumulative proportion of variance. 
#' When data are aligned to an alternative matrix, the statistics 
#' are less straightforward.  A summary of
#' of such an analysis (performed with \code{\link{summary.ordinate}}) 
#' will produce these additional statistics:
#'
#' \itemize{
#' \item{\bold{Singular Value}}{  Rather than eigenvalues, the 
#' singular values from singular value decomposition of the 
#' cross-product of the scaled alignment matrix and the data.}
#' \item{\bold{Proportion of Covariance}}{  Each component's 
#' singular value divided by the sum of singular values.  The cumulative
#' proportion is also returned.  Note that these values do not 
#' explain the amount of covariance between the alignment matrix 
#' and data, but
#' explain the distribution of the covariance.  Large proportions 
#' can be misleading.}
#' \item{\bold{RV by Component}}{  The partial RV statistic by 
#' component.  Cumulative values are also returned.  The sum of partial
#' RVs is Escoffier's RV statistic, which measures the amount of 
#' covariation between the alignment matrix and data.  Caution should
#' be used in interpreting these values, which can vary with the 
#' number of observations and number of variables.  However,
#' the RV is more reliable than proportion of singular value for 
#' interpretation of the strength of linear association for 
#' aligned components.  (It is most analogous to proportion of 
#' variance for principal components.)}

#' }
#' 
#' @param Y An n x p data matrix.
#' @param A An optional n x n symmetric matrix or an n x k data 
#' matrix, where k is the number of variables that could
#' be associated with the p variables of Y.  If NULL, an n x n 
#' identity matrix will be used.
#' @param Cov An optional n x n covariance matrix to describe 
#' the non-independence among
#' observations in Y, and provide a GLS-centering of data.  
#' Note that Cov and A can be the same, if one
#' wishes to align GLS residuals to the same matrix used to 
#' obtain them.  Note also that no explicit GLS-centering
#' is performed on A.  If this is desired, A should be 
#' GLS-centered beforehand.
#' @param transform. An optional argument if a covariance matrix 
#' is provided to transform GLS-centered residuals, if TRUE.  If FALSE, 
#' only GLS-centering is performed.  Only if transform = TRUE 
#' (the default) can one expect the variances of ordinate scores 
#' in a principal component analysis to match eigenvalues.
#' @param scale. a logical value indicating whether the variables 
#' should be scaled to have unit variance before the analysis 
#' takes place. The default is FALSE.
#' @param tol A value indicating the magnitude below which 
#' components should be omitted. (Components are omitted if their 
#' standard deviations are less than or equal to tol times the 
#' standard deviation of the first component.) 
#' With the default null setting, no components are omitted 
#' (unless rank. is provided). 
#' Other settings for tol could be tol = sqrt(.Machine$double.eps), 
#' which would omit essentially constant components, or tol = 0,
#' to retain all components, even if redundant.
#' This argument is exactly the same as in \code{\link{prcomp}}
#' @param rank. Optionally, a number specifying the maximal rank, 
#' i.e., maximal number of aligned components to be used. 
#' This argument can be set as alternative or in addition to tol, 
#' useful notably when the desired rank is considerably 
#' smaller than the dimensions of the matrix.  This argument is 
#' exactly the same as in \code{\link{prcomp}}
#' @param newdata An optional data frame of values for the same 
#' variables of Y to be projected onto 
#' aligned components.  This is only possible with OLS 
#' (transform. = FALSE).
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{ordinate} is a list containing 
#' the following
#' \item{x}{Aligned component scores for all observations}
#' \item{xn}{Optional projection of new data onto components.}
#' \item{d}{The portion of the squared singular values attributed to 
#' the aligned components.}
#' \item{sdev}{Standard deviations of d; i.e., the scale of the components.}
#' \item{rot}{The matrix of variable loadings, i.e. the singular 
#' vectors, \bold{V}.}
#' \item{center}{The OLS or GLS means vector used for centering.}
#' \item{transform}{Whether GLS transformation was used in projection 
#' of residuals 
#' (only possible in conjunction with GLS-centering).}
#' \item{scale}{The scaling used, or FALSE.}
#' \item{alignment}{Whether data were aligned to principal axes or 
#' the name of another matrix.}
#' \item{GLS}{A logical value to indicate if GLS-centering and 
#' projection was used.}
#' @references Collyer, M.L. and D.C. Adams.  2021. 
#' Phylogenetically-aligned Component Analysis. Methods 
#' in Ecology and evolution. In press.
#' @references Revell, L. J.  2009. Size-correction and principal 
#' components for interspecific comparative 
#' studies. Evolution, 63:3258-3268.
#' 
#' @seealso \code{\link{plot.ordinate}}, \code{\link{prcomp}}, 
#' \code{gm.prcomp} within \code{geomorph}
#' @examples
#' 
#' # Examples use residuals from a regression of salamander 
#' # morphological traits against body size (snout to vent length, SVL).
#' # Observations are species means and a phylogenetic covariance matrix
#' # describes the relatedness among observations.
#'
#' data("PlethMorph")
#' Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
#' "Snout.eye", "BodyWidth", 
#' "Forelimb", "Hindlimb")])
#' Y <- as.matrix(Y)
#' R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
#' iter = 0, print.progress = FALSE)$LM$residuals
#' 
#' # PCA (on correlation matrix)
#' 
#' PCA.ols <- ordinate(R, scale. = TRUE)
#' PCA.ols$rot
#' prcomp(R, scale. = TRUE)$rotation # should be the same
#' 
#' # phyPCA (sensu Revell, 2009)
#' # with projection of untransformed residuals (Collyer & Adams 2020)
#' 
#' PCA.gls <- ordinate(R, scale. = TRUE, 
#' transform. = FALSE, 
#' Cov = PlethMorph$PhyCov)
#' 
#' # phyPCA with transformed residuals (orthogonal projection, 
#' # Collyer & Adams 2020)
#' 
#' PCA.t.gls <- ordinate(R, scale. = TRUE, 
#' transform. = TRUE, 
#' Cov = PlethMorph$PhyCov)
#'  
#'  # Align to phylogenetic signal (in each case)
#'  
#'  PaCA.ols <- ordinate(R, A = PlethMorph$PhyCov, scale. = TRUE)
#'  
#'  PaCA.gls <- ordinate(R, A = PlethMorph$PhyCov, 
#'  scale. = TRUE,
#'  transform. = FALSE, 
#'  Cov = PlethMorph$PhyCov)
#'  
#'  PaCA.t.gls <- ordinate(R, A = PlethMorph$PhyCov, 
#'  scale. = TRUE,
#'  transform. = TRUE, 
#'  Cov = PlethMorph$PhyCov)
#'  
#'  # Summaries
#'  
#'  summary(PCA.ols)
#'  summary(PCA.gls)
#'  summary(PCA.t.gls)
#'  summary(PaCA.ols)
#'  summary(PaCA.gls)
#'  summary(PaCA.t.gls)
#'  
#'  # Plots
#'  
#'  par(mfrow = c(2,3))
#'  plot(PCA.ols, main = "PCA OLS")
#'  plot(PCA.gls, main = "PCA GLS")
#'  plot(PCA.t.gls, main = "PCA t-GLS")
#'  plot(PaCA.ols, main = "PaCA OLS")
#'  plot(PaCA.gls, main = "PaCA GLS")
#'  plot(PaCA.t.gls, main = "PaCA t-GLS")
#'  par(mfrow = c(1,1))
#' 
ordinate <- function(Y, A = NULL, Cov = NULL, transform. = TRUE, 
                     scale. = FALSE, 
                     tol = NULL, rank. = NULL, newdata = NULL) {
  
  Ycheck <- try(as.matrix(Y), silent = TRUE)
  if(inherits(Ycheck, "try-error"))
    stop("Y must be a matrix or data frame.\n", call. = FALSE)
  
  id <- if(is.vector(Y)) names(Y) else rownames(Y)
  Y <- as.matrix(Y)
  
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  
  if(is.null(id)) id <- 1:n
  rownames(Y) <- id
  
  alignment <- if(!is.null(A)) 
    try(deparse(substitute(A)), silent = TRUE) else 
      "principal"
  if(length(alignment) != 1) alignment <- "A"

  I <- diag(n)
  
  if(is.null(A)) A <- I
  if(!is.matrix(A))
    stop("A must be a matrix with the same number of rows as Y\n", 
         call. = FALSE)
  if(NROW(A) != n)
    stop("A must be a matrix with the same number of rows as data\n", 
         call. = FALSE)
  if(is.null(rownames(A))) rownames(A) <- rownames(Y) else {
    nnames <- length(intersect(rownames(Y), rownames(A)))
    if(nnames > n)
      stop("The row names of A are not the same as the row names of Y\n",
           call. = FALSE)
    if(isSymmetric(A)) A <- A[id, id] else
      A <- A[id,]
  }
  
  X <- matrix(1, n)
  rownames(X) <- id
  if(!is.null(Cov)) Pcov <- Cov.proj(Cov, rownames(Y))
  cen <- if(is.null(Cov)) colMeans(Y) else 
    lm.fit(as.matrix(Pcov %*% X), as.matrix(Pcov %*% Y))$coefficients
  Z <- as.matrix(scale(Y, center = cen, scale = scale.))
  cen <- attr(Z, "scaled:center")
  sc <- attr(Z, "scaled:scale")
  if (any(sc == 0)) 
    stop("cannot rescale a constant/zero column to unit variance")
  
  k <- if (!is.null(rank.)) {
    stopifnot(length(rank.) == 1, is.finite(rank.), as.integer(rank.) > 
                0)
    min(as.integer(rank.), n, p)
  } else min(n, p)
  
  tf <- if(!is.null(Cov) && transform.) TRUE else FALSE
  if(tf) Z <- as.matrix(Pcov %*% Z)

  rownames(Z) <- id
  
  Saz <-  if(tf || is.null(Cov)) crossprod(A, Z) else crossprod(A, 
                                                        Pcov %*% Z)
  Saz <- as.matrix(Saz)
  
  s <- svd(Saz, nu = 0, nv = k)
  
  if(alignment != "principal") {
    Sa <- crossprod(A)
    Sz <- crossprod(Z)
    RV <- svd(crossprod(A, Z))$d^2 / sqrt(sum(Sa^2) * sum(Sz^2))
  } else RV <- NULL
  
  j <- seq_len(k)
  s$v <- s$v[,j]
  s$d <- s$d[j]
  
  sy <- if(tf || is.null(Cov)) sum(svd(Z)$d^2) else 
    sum(svd(Pcov %*% Z)$d^2)
  s$d <- s$d^2/sum(s$d^2) * sy / max(1, n - 1)
  s$sdev <- sqrt(s$d)

  if (!is.null(tol)) {
    rank <- sum(s$sdev > (s$sdev[1L] * tol))
    if (rank < k) {
      j <- seq_len(k <- rank)
      s$v <- s$v[, j, drop = FALSE]
      s$d <- s$d[j]
      s$sdev <- s$sdev[j]
    }
  }
  
  s$v <- as.matrix(s$v)
  dimnames(s$v) <- list(colnames(Z), paste0("Comp", j))
  
  r <- list(d = s$d, sdev = s$sdev, 
            rot = s$v, 
            center = cen, 
            scale = if(is.null(sc)) FALSE else sc,
            GLS = if(is.null(Cov)) FALSE else TRUE,
            transform = tf,
            alignment = alignment,
            RV = RV)
  
  r$x <- as.matrix(Z %*% s$v)
  rownames(r$x) <- id
  names(r$d) <- names(r$sdev) <- colnames(r$rot)
  if(!is.null(r$RV)) names(r$RV) <- colnames(r$rot)

  if(!is.null(newdata)){
    if(tf) 
    stop("\nNew data cannot be projected onto vectors that were calculated\n
         with respect to the covariances among original observations\n", 
                call. = FALSE)
    xn <- as.matrix(newdata)
    xn <- scale(xn, center = cen, scale = scale.)
    if(NCOL(Z) != NCOL(xn)) 
      stop("Different number of variables in newdata\n", call. = FALSE) 
    r$xn <- xn %*% s$v
    if(!is.null(rownames(newdata)))
      rownames(r$xn) <- rownames(newdata)
  }

  class(r) <- "ordinate"
  r
  
}