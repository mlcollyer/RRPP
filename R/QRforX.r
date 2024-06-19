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
#' @param ... Further arguments passed to base::qr.
#' @return An object of class \code{QR} is a list containing the 
#' following:
#' \item{Q}{The Q matrix.}
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
QRforX <- function(X, ...){
  fix <- FALSE
  S4 <- FALSE
  X <- as.matrix(X)
  p <- NCOL(X)
  if(is.null(colnames(X))) 
    colnames(X) <- paste("V", 1:p, sep = "")
  if(p> 1) {
    Xs <- Matrix(X, sparse = TRUE)
    Xs@x <- round(Xs@x, 12)
    Xs <- Matrix(Xs, sparse = TRUE)
    if(length(Xs@x) < length(X)) X <- Xs
    rm(Xs)
    S4 <- inherits(X, "Matrix")
    QR <- if(S4) suppressWarnings(qr(X, order = 0L)) else 
      suppressWarnings(qr(X, ...))
    if(inherits(QR, "qr")) S4 <- FALSE
    if(S4){
      R <- suppressWarnings(qrR(QR))
      d <- abs(round(diag(R), 12))
      pivot <- which(d > 0)
      rank <- length(pivot)
      if(rank < NCOL(X)) {
        fix <- TRUE
        nms <- dimnames(R)[[2]][pivot]
        X <- X[, nms]
        QR <- suppressWarnings(qr(X, order = 0L))
        R <- suppressWarnings(qrR(QR))
      }
    } else {
      R <- suppressWarnings(qr.R(QR))
      pivot <- with(QR, pivot[1:rank])
      rank <- QR$rank
      if(rank < NCOL(X)) {
        fix <- TRUE
        nms <- dimnames(R)[[2]][pivot]
        X <- X[, nms]
        QR <- suppressWarnings(qr(as.matrix(X)))
        R <- suppressWarnings(qr.R(QR))
      }
    }
  } else {
    QR <- qr(X)
    rank <- 1
    pivot <- 1
    R <- qr.R(QR)
  }
  
  Xnms <- colnames(X)
  
  if(S4) {
    QR <- qr(X)
    R <- qrR(QR)
  }
    
  Xnms <- colnames(X)
  Q <- qr.Q(QR)
  Q@x <- round(Q@x, 12)
  Q <- Matrix(Q, sparse = TRUE)
  
  if(!all.equal(dimnames(R)[[2]], Xnms)){
    nnms <- match(dimnames(R)[[2]], Xnms)
    R <- R[nnms, nnms]
    Q <- Q[, nnms]
  }

  out <- list(Q = Q, R = R, X = X,
              rank = rank, fixed = fix, S4 = S4)
  
  out$dimnames <- list(dimnames(Q)[[1]], dimnames(R)[[2]])
  
  class(out) <- "QR"
  out
  
}
