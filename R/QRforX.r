#' QR decomposition of linear model design matrices
#' 
#' This function performs a QR decomposition (factorization) on a linear 
#' model design matrix (X) and returns useful results for subsequent analysis.
#' This is intended as an internal function but can be used externally.  Because
#' base::qr and Matrix::qr have different options for QR algorithms, this
#' function assures that results are consistent for other RRPP function use, 
#' whether X is a dense or sparse matrix.
#'
#' @param X A linear model design matrix, but can be any object coercible to matrix.
#' @param returnQ A logical value whether to return the Q matrix.  Generating a
#' Q matrix can be computationally intense for large matrices.  If it is not
#' explicitly needed, this argument can be FALSE.
#' @param ... Further arguments passed to base::qr.
#' @return An object of class \code{QR} is a list containing the 
#' following:
#' \item{Q}{The Q matrix, if requested.}
#' \item{R}{The R matrix.}
#' \item{X}{The X matrix, which could be changed from dense to sparse,
#' or vice versa, and redundant columns removed.}
#' \item{dimnames}{The dimnames associated with Q, R, or X}
#' \item{rank}{The rank of the X matrix.}
#' \item{S4}{Logical value for whether Q, R, and X are S4 class objects.}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples
#' ## Simple Example
#' data(Pupfish)
#' fit <- lm.rrpp(coords ~ Pop, data = Pupfish, print.progress = FALSE)
#' QR <- QRforX(model.matrix(fit))
#' QR$Q
#' QR$R
#' QR$rank
#' QR$S4
#' 
#' ## Not run, but one could get base::qr and Matrix::qr results as
#' 
#' # base::qr(as.matrix(QR$X))
#' # Matrix::qr(QR$X)
#' 
#' ## Complex example
#' 
#' data("PupfishHeads")
#' fit <- suppressWarnings(lm.rrpp(headSize ~ sex + 
#' locality/year, data = PupfishHeads))
#' X <- model.matrix(fit)
#' dim(X) # Already reduced
#' colnames(X)
#' X <- model.matrix(terms(fit), fit$LM$data)
#' dim(X) # Retains redundant parameters
#' colnames(X)
#' QR <- QRforX(X)
#' dim(QR$X)
#' dim(QR$Q)
#' dim(QR$R)
#' QR$rank
#' 
QRforX <- function(X, returnQ = TRUE,
                   ...){
  
  cs <- colSums(as.matrix(X))
  
  if(NCOL(X) > 1 && any(cs == 0))
    X <- X[, cs != 0, drop = FALSE]
  
  if(!isS4(X) && NCOL(X) > 1){
    if(length(which(X != 0)) <= 0.9 * length(X))
      X <- Matrix(X, sparse = TRUE)
  }
  
  if(isS4(X)){
    if(length(which(X != 0)) > 0.9 * length(X))
      X <- as.matrix(X)
  }
  
  dims <- dim(X)
  n <- dims[1]
  k <- dims[2]
  
  use1 <- TRUE
  if(NCOL(X) > 1 && isS4(X)) use1 <- FALSE
  
  QR1 <- function(X){
    QR <- qr(X)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    X <- X[, QR$pivot[1:QR$rank], drop = FALSE]
    out <- list(Q = Q, R = R, X = X)
    out
  }
  
  QR2 <- function(X){
    QR <- suppressWarnings(qr(X))
    R <- suppressWarnings(qrR(QR))
    Q <- suppressWarnings(qr.Q(QR))
    QR <- suppressWarnings(qr(as.matrix(R)))
    if(QR$rank < NCOL(X)) {
      X <- X[,QR$pivot[1:QR$rank], drop = FALSE]
      QR <- suppressWarnings(qr(X))
      R <- suppressWarnings(qrR(QR))
      Q <- suppressWarnings(qr.Q(QR))
    }
    Q <- drop0(Q, tol = 1e-10)
    R <- drop0(R, tol = 1e-10)
    X <- drop0(X, tol = 1e-10) 
    out <- list(Q = Q, R = R, X = X)
    out
  }
  
  QR <- try(if(use1) QR1(X) else QR2(X),
            silent = TRUE)
  if(inherits(QR, "try-error")) {
  
    qrt <- qr(as.matrix(X))
    X <- X[, qrt$pivot[1:qrt$rank], drop = FALSE]
    if(length(which(X != 0)) <= 0.9 * length(X))
      X <- Matrix(X, sparse = TRUE)
    if(NCOL(X) > 1 && isS4(X)) use1 <- FALSE
    QR <- if(use1) QR1(as.matrix(X)) else QR2(X)
  
  }
    
  QR$rank <- NCOL(QR$X)
  if(sum(QR$R) == 0) QR$rank <- 0
  QR$S4 <- isS4(QR$X)
  QR$dimnames <- dimnames(QR$X)
  
  if(!returnQ)
    QR$Q <- NULL
  
  class(QR) <- "QR"
  QR
}
