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
#' @param reduce A logical value for whether redundant parameters in X should be 
#' removed.  This should be TRUE (default) for most cases.
#' @param reQR A logical value for whether to re-perform QR if reduce = TRUE,
#' and X has been reduced.
#' @param ... Further arguments passed to base::qr.
#' @return An object of class \code{QR} is a list containing the 
#' following:
#' \item{Q}{The Q matrix, if requested.}
#' \item{R}{The R matrix.}
#' \item{X}{The X matrix, which could be changes from dense to sparse,
#' or vice versa, and redundant columns removed.}
#' \item{rank}{The rank of the X matrix.}
#' \item{fix}{Logical value for whether redundant columns were removed
#' form X.  TRUE means columns were removed.}
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
#' QR$fixed
#' dim(QR$X) # Reduced again
#' colnames(QR$X)
#' 
QRforX <- function(X, returnQ = TRUE,
                   reduce = TRUE, reQR = TRUE,
                   ...){
  fix <- FALSE
  S4 <- FALSE
  rank <- NULL
  pivot <- NULL
  X <- as.matrix(X)
  p <- NCOL(X)
  if(is.null(colnames(X)) && p >= 1) 
    colnames(X) <- paste("V", 1:p, sep = "")
  
  if(p <= 1){
    QR <- qr(X)
    rank <- QR$rank
    pivot <- QR$pivot
    Q <- if(returnQ) qr.Q(QR) else NULL
    R <- qr.R(QR)
  }
  
  if(p > 1) {
    if(reduce){
      Xs <- Matrix(X, sparse = TRUE)
      Xs@x <- round(Xs@x, 8)
      Xs <- Matrix(Xs, sparse = TRUE)
      if(length(Xs@x) < length(X)) X <- Xs
      rm(Xs)
      S4 <- inherits(X, "Matrix")
      QR <- if(S4) try(suppressWarnings(qr(X, order = 0L)), 
                       silent = TRUE) else
        suppressWarnings(qr(X, ...))
      if(inherits(QR, "try-error")) QR <- qr(as.matrix(X), ...)
      if(inherits(QR, "qr")) S4 <- FALSE
      
      if(S4){
        R <- suppressWarnings(qrR(QR))
        d <- abs(round(diag(R), 8))
        pivot <- which(d > 0)
        rank <- length(pivot)
        if(rank < NCOL(X)) {
          fix <- TRUE
          nms <- dimnames(R)[[2]][pivot]
          X <- X[, nms]
        }
      } else {
        R <- suppressWarnings(qr.R(QR))
        pivot <- with(QR, pivot[1:rank])
        rank <- QR$rank
        if(rank < NCOL(X)) {
          fix <- TRUE
          nms <- dimnames(R)[[2]][pivot]
          X <- X[, nms]
        }
      }
      Q <- if(returnQ) qr.Q(QR) else NULL
    }
    
    Xnms <- colnames(X)
    Xs <- Matrix(X, sparse = TRUE)
    Xs@x <- round(Xs@x, 8)
    Xs <- Matrix(Xs, sparse = TRUE)
    if(length(Xs@x) < length(X)) X <- Xs
    rm(Xs)
    S4 <- inherits(X, "Matrix")
    
    if(!reduce) reQR <- TRUE
    
    if(reQR) {
      QR <- qr(X)
      if(returnQ){
        Q <- Matrix(qr.Q(QR), sparse = TRUE)
        Q@x <- round(Q@x, 8)
        Q <- Matrix(Q, sparse = TRUE)
        if(length(Q@x) == length(Q)) Q <- as.matrix(Q)
      } else Q <- NULL
      
      R <- if(S4) qrR(QR) else qr.R(QR)
      Rs <- Matrix(R, sparse = TRUE)
      Rs@x <- round(Rs@x, 8)
      Rs <- Matrix(Rs, sparse = TRUE)
      if(length(Rs@x) < length(R)) R <- Rs
      rm(Rs)
      
      if(!all.equal(dimnames(R)[[2]], Xnms)){
        nnms <- match(dimnames(R)[[2]], Xnms)
        R <- R[nnms, nnms]
        if(returnQ) Q <- Q[, nnms]
      }
    }
  }
  
  if(is.null(rank)) rank <- NCOL(X)
  if(is.null(pivot)) pivot <- 1:rank
  
  out <- list(Q = Q, R = R, X = X,
              rank = rank, fixed = fix, S4 = S4)
  
  out$dimnames <- list(dimnames(Q)[[1]], dimnames(R)[[2]])
  
  class(out) <- "QR"
  out
  
}
