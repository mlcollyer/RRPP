
## logLik.lm.rrpp

#' Calculate the log-likelihood of a lm.rrpp fit
#' 
#' \code{logLik.lm.rrpp} returns the log-likelihood of
#' an \code{lm.rrpp} object.  Ridge regularization will be performed for 
#' ill-conditioned or singular residual covariance matrices, but dimension
#' reduction could be augmented via projection, using the arguments, tol
#' and pc.no.  See \code{\link{ordinate}} for details.
#'
#' @param object Object from \code{\link{lm.rrpp}}
#' @param tol A value indicating the magnitude below which 
#' components should be omitted, following projection.  See \code{\link{ordinate}} 
#' for details.
#' @param pc.no Optionally, a number specifying the maximal number of 
#' principal components, passed onto \code{\link{ordinate}}, as argument, rank.
#' @param Z A logical value for whether to calculate Z scores based on RRPP.
#' @param gls.null A logical value for if a fit has a GLS estimation, should
#' the null model (intercept) also have a GLS estimation, for estimating Z.  
#' Should be FALSE if the log-likelihood is measured to compare different GLS
#' estimations for a covariance matrices
#' @param verbose A logical value for whether to return random log-likelihood values,
#' if Z-scores are calculated.
#' @param ...	further arguments passed to or from other methods
#' @export
#' @author Michael Collyer
#' @keywords utilities

logLik.lm.rrpp <- function(object, tol = NULL, 
                            pc.no = NULL, Z = TRUE,
                            verbose = FALSE,
                            gls.null = FALSE, ...){
  
  if(isTRUE(object$LM$LMM))
    stop("A log-likelihood for lmm.rrpp objects is not yet available.\n",
         call. = FALSE)
  
  ll <- .logLik.lm.rrpp(object, tol = NULL,
                        pc.no = NULL, Z = FALSE,
                        verbose = FALSE,
                        gls.null = FALSE, ...)
  val <- ll$logL
  k <- attr(ll, "k") 
  p <- attr(ll, "p")
  attr(val, "df") <- p * k + 0.5 * p* (p + 1)
  attr(val, "nall") <- attr(ll, "nall") 
  attr(val, "nobs") <- attr(ll, "nobs") 
  class(val) <- "logLik"
  val
}

.logLik.lm.rrpp <- function(object, tol = NULL, 
                           pc.no = NULL, Z = TRUE,
                           verbose = FALSE,
                           gls.null = FALSE, ...){
  
  if(is.null(tol)) tol = 0
  
  PermInfo <- getPermInfo(object, attribute = "all")
  ind <- PermInfo$perm.schedule
  gls <- object$LM$gls
  n <- object$LM$n
  p <- object$LM$p.prime
  X <- as.matrix(object$LM$X)
  Y <- as.matrix(object$LM$Y)
  R <- if(gls) object$LM$gls.residuals else object$LM$residuals
  PCA <- ordinate(R, tol = tol, rank. = min(c(pc.no, p)))
  rnk <- length(PCA$d)
  w <- object$LM$weights
  Cov <- object$LM$Cov  
  if(!is.null(object$LM$Cov) && is.null(object$LM$Pcov)) 
    object$LM$Pcov <- Cov.proj(object$LM$Cov)
  Pcov <- NULL
  if(!is.null((object$LM$Pcov)))
    Pcov <- object$LM$Pcov
  
  PY <- if(gls) {
    if(!is.null(Pcov)) Pcov %*% Y else Y * sqrt(w)
    } else NULL
  PX <- if(gls) {
    if(!is.null(Pcov)) Pcov %*% X else X * sqrt(w)
  } else NULL
    
  U <- if(gls) QRforX(PX, reduce = FALSE)$Q else 
    QRforX(X, reduce = FALSE)$Q
  
  TY <- if(gls.null) PY else Y
  
  ll <- logL(object) 
  
  if(Z) {
    
    perms <- length(ind)
    
    if(!is.null(w)) {
      excl <- w <= 0
      logdetC <- sum(log(w[!excl]))
    } else {
      logdetC <- if(gls) determinant(Cov, logarithm = TRUE)$modulus else 0
    }
    
    ll.list <- sapply(1:perms, function(j){
      y <- TY[ind[[j]], ]
      if(gls && !gls.null) {
        ty <- Pcov %*% y
        r <- ty - fastFit(U, ty, n, rnk)
      } else r <- y - fastFit(U, y, n, rnk)
      r <- as.matrix(r)
      if(NCOL(r) > rnk) r <- ordinate(r, rank. = rnk)$x
      Sig <- as.matrix(crossprod(r) / n)
      if(kappa(Sig) > 1e10) Sig <- RiReg(Sig, r)
      logdetSig <- determinant(Sig, logarithm = TRUE)$modulus
      -0.5 * (n * rnk * log(2 * pi) + rnk * logdetC +
                n * logdetSig + n * rnk) 
    })
     
    z <- effect.size(ll.list)
    ll <- c(ll, Z = z)
    if(verbose) ll <- list(ll = ll, Z = z, 
                           random.logL = ll.list)
  }
  
  attr(ll, "k") <- ncol(X)
  attr(ll, "p") <- p
  attr(ll, "nall") <- n
  attr(ll, "nobs") <- n
  
  return(ll)
    
}
