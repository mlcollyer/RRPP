#' Test of beta parameters for lm.rrpp model fits
#'
#' @description For any \code{\link{lm.rrpp}} object, a vector of coefficients
#' can be used for a specific test of a vector of betas (specific population parameters).
#' This test follows the form of (b - beta) in the numerator of a t-statistic, where 
#' beta can be a value other than 0 (or 0).  However, for this test, a vector (Beta) of length, p,
#' is used for the p variables in the \code{\link{lm.rrpp}} fit.  If Beta is a vector of 0s, this
#' test is essentially the same as the test performed for \code{\link{coef.lm.rrpp}}.  
#' However, there is one key difference: a Mahalanobis distance rather than an Euclidean distance
#' is calculated as a test statistic.  This assures some consistency across random permutations,
#' by accounting for variation in residual covariance matrices, especially because alternative
#' null models are permitted.
#' 
#' This function is not recommended for high-dimensional data, as it will invoke
#' use of a generalized inverse of the residual covariance matrix.  Additionally, 
#' hypothesizing the beta values for many variables should be challenging.  However,
#' this test will work well with one or few variables.
#' 
#'
#' @param fit Object from \code{\link{lm.rrpp}}
#' @param fit.null Optional object from \code{\link{lm.rrpp}} to use as a null model.
#' @param coef.no. The row of a matrix of coefficients for which to perform the test.  
#' This can be learned by performing coef(fit), prior to the test.  If left NULL, 
#' the last coefficients vector will be used.
#' @param Beta A single value (for univariate data) or a numeric vector with length equal to
#' the number of variables used in the fit object.  If left NULL, 0 is used for each parameter.
#' @return Function returns a list with the following components: 
#'   \item{obs.d}{Length of observed b - Beta vector} 
#'   \item{obs.md}{The observed b - Beta vector length, after accounting for 
#'   residual covariance matrix; the Mahalanobis distance}
#'   \item{Beta}{Hypothesized beta values in the Beta vector.}
#'   \item{obs.B.mat}{The observed matrix of coefficients (before subtracting Beta).}
#'   \item{coef.no}{The row of the observed matrix of coefficients, for which to subtract Beta.}
#'   \item{random.md}{Random Mahalanobis distances produced with RRPP.}
#'   \item{Z}{The effect size of the observed Mahalanbois distance, based on RRPP.}
#'   \item{P}{The P-value of the observed Mahalanobis distance, based on RRPP.}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' \dontrun{
#' data(PlethMorph)
#' fit <- lm.rrpp(TailLength ~ SVL, 
#' data = PlethMorph,
#' verbose = TRUE)
#' 
#' # Allometry test (Beta = 0)
#' T1 <- betaTest(fit, coef.no = 2, Beta = 0)
#' summary(T1)
#' 
#' # compare to
#' coef(fit, test = TRUE)
#' 
#' # Isometry test (Beta = 1)
#' # Failure to reject H0 suggests isometric-like association.
#' 
#' T2 <- betaTest(fit, coef.no = 2, Beta = 1)
#' summary(T2)
#' 
#' # More complex tests
#' 
#' fit2 <- lm.rrpp(HeadLength ~ SVL + TailLength, 
#' data = PlethMorph,
#' SS.type = "II",
#' verbose = TRUE)
#' 
#' ## allometries
#' T3 <- betaTest(fit2, coef.no = 2, Beta = 0)
#' T4 <- betaTest(fit2, coef.no = 3, Beta = 0)
#' summary(T3)
#' summary(T4)
#' 
#' # compare to
#' coef(fit2, test = TRUE)
#' 
#' ## isometries
#' T5 <- betaTest(fit2, coef.no = 2, Beta = 1)
#' T6 <- betaTest(fit2, coef.no = 3, Beta = 1)
#' summary(T5)
#' summary(T6) 
#' 
#' PlethMorph$Y <- cbind(PlethMorph$HeadLength, PlethMorph$TailLength)
#' fit3 <- lm.rrpp(Y ~ SVL, 
#' data = PlethMorph,
#' verbose = TRUE)
#' 
#' T7 <- betaTest(fit3, coef.no = 2, Beta = c(0, 0))
#' T8 <- betaTest(fit3, coef.no = 2, Beta = c(1, 1))
#' 
#' summary(T7)
#' summary(T8)
#' 
#' }

betaTest <- function(fit, fit.null = NULL,
                     coef.no = NCOL(fit$LM$X),
                     Beta = rep(0, fit$LM$p)){
  
  n <- fit$LM$n
  p <- fit$LM$p
  k <- NCOL(fit$LM$X)
  
  QR <- if(!is.null(fit.null))
    getModels(fit.null, "qr")$full else
      getModels(fit, "qr")$reduced
  
  QR <- QR[[length(QR)]]
  
  gls <- fit$LM$gls
  if(gls && !is.null(fit$LM$Cov)) {
    Pcov <- getModelCov(fit, type = "Pcov")
    TY <- Pcov %*% fit$LM$Y
  } else if(gls && is.null(fit$LM$Cov)) {
    TY <- fit$LM$Y * sqrt(fit$LM$weights)
  } else TY <- fit$LM$Y
  
  Fitted <- as.matrix(fastFit(QR$Q, TY, n, p))
  Resid <- as.matrix(TY - Fitted)
  rm(QR)
  
  PI <- getPermInfo(fit, "all")
  ind <- PI$perm.schedule
  
  Bobs <- coef(fit)

  if(coef.no > k)
    stop("The coef.no specificed is greater than the number of available coefficient vectors.\n",
         call. = FALSE)
  if(coef.no <= 0)
    stop("The coef.no must be a positive integer.\n",
         call. = FALSE)
  
  if(!is.vector(Beta))
    stop("Beta must be a numeric vector.\n", call. = FALSE)
  if(length(Beta) != length(Bobs[coef.no, ]))
    stop("Beta must have the same length as number of variables in the model fit.\n", 
         call. = FALSE)
  
  names(Beta) <- colnames(fit$LM$Y)
  
  Qf <- getModels(fit, "qr")$full
  Qf <- Qf[[length(Qf)]]
  Hb <- getHb(Qf)
  Hbs <- drop0(Matrix(Hb, sparse = TRUE), tol = 1e-7)
  if(length(Hbs@x) < length(Hb)) Hb <- Hbs
  rm(Hbs)
  Qf <- Qf$Q
  
  Result <- sapply(ind, function(x){
    Yr <- Fitted + Resid[x,]
    Bi <- as.matrix(Hb %*% Yr)
    b <- Bi[coef.no,] 
    R <- Yr - fastFit(Qf, Yr, n, p)
    S <- fast.solve(crossprod(R) / (n - k))
    as.numeric(sqrt(crossprod(b, S) %*% b))
  })
  
  obs.d <- sqrt(crossprod(Bobs[coef.no, ]))
  b <- Bobs[coef.no, ] - Beta
  S.obs <- fast.solve(crossprod(Resid) / (n - k))
  obs.md <- as.numeric(sqrt(crossprod(b, S.obs) %*% b))
  Result[[1]] <- obs.md

    out <- list(obs.d = obs.d,
                obs.md = obs.md, 
                Beta = Beta,
                obs.B.mat = Bobs,
                coef.no = coef.no,
                random.md = Result, 
                Z = effect.size(Result),
                P = pval(Result))
    attr(out, "class") <- "betaTest"
    if(confidence <= 0 || confidence > 1)
      confidence = 0.95
    out$confidence = confidence
    out
}
