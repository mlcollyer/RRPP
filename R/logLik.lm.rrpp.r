
## logLik.lm.rrpp

#' Calculate the log-likelihood of a lm.rrpp fit
#' 
#' \code{logLik.lm.rrpp} returns the log-likelihood of
#' an \code{lm.rrpp} object or an \code{lmm.rrpp} object.  
#' Ridge regularization will be performed for 
#' ill-conditioned or singular residual covariance matrices, but dimension
#' reduction could be augmented via projection, using the arguments, tol
#' and pc.no.  See \code{\link{ordinate}} for details.
#' 
#' Note that options are limited for \code{lmm.rrpp} objects.  It is currently
#' not possible to use RRPP to find Z scores.
#' 
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
  
  ll <- if(isTRUE(object$LM$LMM))
    .logLik.lmm.rrpp(object, tol = NULL,
                     pc.no = NULL, Z = FALSE,
                     verbose = FALSE,
                     gls.null = FALSE, ...) else
            .logLik.lm.rrpp(object, tol = NULL,
                        pc.no = NULL, Z = FALSE,
                        verbose = FALSE,
                        gls.null = FALSE, ...)
  
  if(is.list(ll)) val <- ll[[1]] else val <- as.numeric(ll)
  k <- attr(ll, "k") 
  p <- attr(ll, "p")
  attr(val, "df") <- p * k + 0.5 * p* (p + 1)
  attr(val, "nall") <- attr(ll, "nall") 
  attr(val, "nobs") <- attr(ll, "nobs") 
  attr(val, "p") <- attr(ll, "p") 
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


.logLik.lmm.rrpp <- function(object, tol = NULL, 
                             pc.no = NULL, Z = TRUE,
                             verbose = FALSE,
                             gls.null = FALSE, ...){
  
  if(is.null(tol)) tol = 0
  
  PermInfo <- getPermInfo(object, attribute = "all")
  ind <- PermInfo$perm.schedule
  gls <- object$LM$gls
  n <- object$LM$n
  p <- object$LM$p.prime
  X <- QRforX(object$LM$X)$X
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
  R <- if(gls) {
    if(!is.null(Pcov)) Pcov %*% R else R * sqrt(w)
  } else R
  
  RtR_det <- try(det(crossprod(R)), silent = TRUE)
  if(!inherits(RtR_det, "try-error")) 
    object$devcomp$RtR_det <- RtR_det else
      stop("The residuals appear to be singular using the current attributes.\n",
           call. = FALSE)
  # update U
  # get hresid
  XZ <- if(gls) cbind(PX, object$LM$Z %*% object$LM$Lambda) else
    cbind(X, object$LM$Z %*% object$LM$Lambda)
  Hb <- if(gls) getLMM_Hb(PX, object$LM$Z %*% object$LM$Lambda) else
    getLMM_Hb(X, object$LM$Z %*% object$LM$Lambda)
  H <- XZ %*% Hb  - diag(n)
  
  TY <- if(gls.null) PY else Y
    
  k <- NCOL(X)
  
  est <- object$LM$estimation
  devcomp <- object$devcomp
  r <- devcomp$RtR_det 
  if(!is.na(devcomp$UtU_det)) r <- r + devcomp$UtU_det
  lres <- log(2 * pi * r)
  n <- object$LM$n
  p <- object$LM$p.prime
  np <- n * p
  if(est == "REML") npk <- np - 
    length(object$LM$coef.fixed) * object$LM$p.prime else 
      npk <- np
  
  ll <- -0.5 * p * ((devcomp$ldLZ + devcomp$ldLX - devcomp$ldO) +
                      npk * (1 + log((2 * pi * r) / npk)))
  if(Z) {
    
    perms <- length(ind)
    
    ll.list <- sapply(1:perms, function(j){
      y <- TY[ind[[j]], ]
      if(gls && !gls.null) {
        ty <- Pcov %*% y
        r <- H %*% ty
      } else r <- H %*% y
      r <- as.matrix(r)
      if(NCOL(r) > rnk) r <- ordinate(r, rank. = rnk)$x
      r <- det(crossprod(r))
      if(!is.na(devcomp$UtU_det)) r <- r + devcomp$UtU_det
      -0.5 * p * ((devcomp$ldLZ + devcomp$ldLX - devcomp$ldO) +
                    npk * (1 + log((2 * pi * r) / npk)))
    })
    
    z <- effect.size(ll.list)
    ll <- c(ll, Z = z)
    if(verbose) ll <- list(ll = ll, Z = z, 
                           random.logL = ll.list)
  }
  
  
  attr(ll, "k") <- ncol(X) + 
    length(unique(as.vector(object$LM$Lambda))) -1
  attr(ll, "p") <- p
  attr(ll, "nall") <- n
  attr(ll, "nobs") <- n 
  
  ll
}
