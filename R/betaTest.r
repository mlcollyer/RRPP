#' Test of beta parameters for lm.rrpp model fits
#'
#' @description For any \code{\link{lm.rrpp}} object, a vector of coefficients
#' can be used for a specific test of a vector of betas (specific population parameters).
#' This test follows the form of (b - beta) in the numerator of a t-statistic, where 
#' beta can be a value other than 0 (or 0).  However, for this test, a vector (Beta) of length, p,
#' is used for the p variables in the \code{\link{lm.rrpp}} fit.  If Beta is a vector of 0s, this
#' test is essentially the same as the test performed for \code{\link{coef.lm.rrpp}}.  However,
#' it is possible to test null hypotheses for beta values other than 0, sensu Cicuéndez et al. (2023).
#' 
#' This function can use either the square-root of the inner-product of vectors of coefficients (distance, d)
#' or generalized inner-product based on the inverse of the residual covariance matrix 
#' (Mahalanobis distance, md) as statistics.  In most cases, either will likely yield similar (or same) 
#' P-values.  However, Mahalanobis distance might be preferred for generalized least squares fits, which 
#' do not have consistent residual covariance matrices for null (intercept only) models (the distances are
#' standarized by the residual covariances).  If high-dimensional data are analyzed, a generalized inverse 
#' of the residual covariance matrix will be used because of singular covariance matrices.  Results are less 
#' trustworthy with Mahalanbois distances, in these cases.
#' 
#' The coefficient number should be provided for specific tests.  One can determine this with, e.g., 
#' coef(fit).  If it is not provided (NULL), tests will be performed on all possible vectors of coefficients
#' (rows of coefficients matrix).  These tests will be performed sequentially.  If a null model is not specified,
#' then for each vector of coefficients, the corresponding parameter is dropped from the linear model
#' design matrix to make a null model.  This process is analogous in some ways to a leave-one-out 
#' cross-validation (LOOCV) analysis, testing each coefficient against models containing parameters for all other
#' coefficients.  For example, for a linear model fit, y ~ x1 + x2 + 1, where x1 and x2 are single-parameter 
#' covariates,
#' the analysis would first drop the intercept, then x1, then x2, performing three sequential analyses.  This 
#' option could require large amounts of computation time for large models, high-dimensional data, many RRPP
#' permutations, or any combination of these.  It is important to realize that if x1 and x2 are multi-parameter 
#' factors, betaTest will treat each parameter as unique, unless a specific null model is provided.
#' This is different than coef.lm.rrpp, which uses terms of the model based on SS type, and is, therefore, more
#' complex yet less flexible.
#' 
#'  \subsection{Difference between coef.lm.rrpp test and betaTest}{ 
#'  The test for coef.lm.rrpp uses the square-root of inner-products of vectors (d) as a 
#'  test statistic and only tests the null hypothesis that the length of the vector is 0.  
#'  The significance of the test is based on random values produced by RRPP, based on the
#'  matrices of coefficients that are produced in all permutations.  The null models for generating
#'  RRPP distributions are consistent with those used for ANOVA, as specified in the 
#'  \code{\link{lm.rrpp}} fit by choice of SS type.  Therefore, the random coefficients are
#'  consistent with those produced by RRPP for generating random estimates used in ANOVA.
#'  
#'  The betaTest analysis allows different null hypotheses to be used (vector length is not necessarily 0) 
#'  and unless otherwise specified, uses a null model that lacks one vector of parameters and a full
#'  model that contains all vectors of parameters, for the parameter for which coefficients are estimated.
#'  This is closest to a type III SS method of estimation, but each parameter is dropped from the model,
#'  rather than terms potentially comprising several parameters.  Additionally, betaTest calculates
#'  Mahalanobis distance, in addition to Euclidean distance, for vectors of coefficients.  
#'  This statistic is probably 
#'  better for more types of models (like generalized least squares fits).
#' }
#' 
#' @param fit Object from \code{\link{lm.rrpp}}
#' @param X.null Optional object that is either a linear model design matrix or a model
#' fit from \code{\link{lm.rrpp}}, from which a linear model design matrix can be extracted. Note
#' that if any transformation of a design matrix is required (GLS estimation), 
#' it is assumed that the matrix was transformed prior to analysis.  If X.null is a \code{\link{lm.rrpp}}
#' object, transformation is inherent.
#' @param coef.no. The row or rows of a matrix of coefficients for which to perform the test.  
#' This can be learned by performing coef(fit), prior to the test.  If left NULL, 
#' the analysis will cycle through every possible vector of coefficients (rows of a coefficients matrix).
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
#' @references Tejero-Cicuéndez, H., I. Menéndez, A. Talavera, G. Riaño, B. Burriel-Carranza, 
#' M. Simó-Riudalbas, S. Carranza, and D.C. Adams. 2023. 
#' Evolution along allometric lines of least resistance: 
#' Morphological differentiation in Pristurus geckos. Evolution. 77:2547–2560.
#' @seealso \code{\link{coef.lm.rrpp}}
#' @examples 
#' \dontrun{
#' data(PlethMorph)
#' fit <- lm.rrpp(TailLength ~ SVL, 
#' data = PlethMorph,
#' verbose = TRUE)
#' 
#' ## Allometry test (Beta = 0)
#' T1 <- betaTest(fit, coef.no = 2, Beta = 0)
#' summary(T1)
#' 
#' # compare to
#' coef(fit, test = TRUE)
#' 
#' # Note that is Beta is not provided
#' 
#' T1 <- betaTest(fit, coef.no = 2)
#' summary(T1)
#' 
#' # Note that is coef.no is not provided
#' 
#' T1 <- betaTest(fit)
#' summary(T1)
#' 
#' 
#' ## Isometry test (Beta = 1)
#' # Failure to reject H0 suggests isometric-like association.
#' 
#' T2 <- betaTest(fit, coef.no = 2, Beta = 1)
#' summary(T2)
#' 
#' 
#' ## More complex tests
#' 
#' # Multiple covariates
#' 
#' fit2 <- lm.rrpp(HeadLength ~ SVL + TailLength, 
#' data = PlethMorph,
#' SS.type = "II",
#' verbose = TRUE)
#' 
#' fit.null1 <- lm.rrpp(HeadLength ~ SVL, 
#' data = PlethMorph,
#' verbose = TRUE)
#' 
#' fit.null2 <- lm.rrpp(HeadLength ~ TailLength, 
#' data = PlethMorph,
#' verbose = TRUE)
#' 
#' ## allometries
#' T3 <- betaTest(fit2, fit.null2, coef.no = 2, Beta = 0)
#' T4 <- betaTest(fit2, fit.null1, coef.no = 3, Beta = 0)
#' summary(T3)
#' summary(T4)
#' 
#' # compare to
#' coef(fit2, test = TRUE)
#' 
#' ## isometries
#' T5 <- betaTest(fit2, fit.null2, coef.no = 2, Beta = 1)
#' T6 <- betaTest(fit2, fit.null1, coef.no = 3, Beta = 1)
#' summary(T5)
#' summary(T6) 
#' 
#' # General test of slopes = 0 for all coefficients
#' T7 <- betaTest(fit)
#' summary(T7)
#' 
#' # Compare to
#' coef(fit2, test = TRUE)
#' 
#' # multivariate data
#' 
#' PlethMorph$Y <- cbind(PlethMorph$HeadLength, PlethMorph$TailLength)
#' fit3 <- lm.rrpp(Y ~ SVL, 
#' data = PlethMorph,
#' verbose = TRUE)
#' 
#' T8 <- betaTest(fit3, coef.no = 2, Beta = c(0, 0))
#' T9 <- betaTest(fit3, coef.no = 2, Beta = c(1, 1))
#' 
#' summary(T8)
#' summary(T9)
#' 
#' }

betaTest <- function(fit, X.null = NULL,
                     coef.no = NULL,
                     Beta = NULL,
                     print.progress = FALSE
                     ){
  
  n <- fit$LM$n
  p <- fit$LM$p
  k <- NCOL(fit$LM$X)
  if(is.null(Beta))
    Beta <- rep(0, p)
  
  if(is.null(coef.no)) coef.no <- 1:k
  coef.no <- as.vector(coef.no)
  kk <- length(coef.no)
  
  if(any(coef.no > k))
    stop("coef.no contains values greater than the number of coefficients.\n",
         call. = FALSE)
  if(any(coef.no <= 0))
    stop("coef.no must consiste of positive integers.\n",
         call. = FALSE)

  Bobs <- coef(fit)
  
  if(!is.vector(Beta))
    stop("Beta must be a numeric vector.\n", call. = FALSE)
  if(length(Beta) != p)
    stop("Beta must have the same length as number of variables in the model fit.\n", 
         call. = FALSE)
  
  names(Beta) <- colnames(fit$LM$Y)
  
  Xf <- as.matrix(fit$LM$X)
  Y <- as.matrix(fit$LM$Y)
  gls <- fit$LM$gls
  Pcov <- if(gls && !is.null(fit$LM$Cov)) 
    getModelCov(fit, type = "Pcov") else NULL
  w <- if(gls && is.null(fit$LM$Cov)) 
    sqrt(fit$LM$weights) else NULL
  
  if(gls){
    if(is.null(w)) {
      Xf <- Pcov %*% Xf 
      TY <- Pcov %*% Y
    } else {
      Xf <- Xf * w
      TY <- Y * w 
      }
    } else TY <- Y
  
  QRf <- qr(Xf)
  Qf <- qr.Q((QRf))
  Rf <- qr.R(QRf)
  Hb <- as.matrix(tcrossprod(fast.solve(Rf), Qf))
  Hbs <- drop0(Matrix(Hb, sparse = TRUE), 1e-7)
  if(length(Hbs@x) < length(Hb)) Hb <- Hbs
  rm(Hbs)
  
  coef.nms <- colnames(Xf)[coef.no]
  
  userNULL <- FALSE
  
  if(!is.null(X.null)) {
    if(inherits(X.null, "lm.rrpp")){
      QRr <- getModels(X.null, "qr")
      QRr <- QRr$full[[length(QRr$full)]]
    } else {
      QRr <- QRforX(X.null)
    }
    userNULL <- TRUE
    Qr <- QRr$Q
    Q.list <- lapply(1:kk, function(.) Q = Qr)
  }
  
  if(!userNULL){
    Q.list <- lapply(1:kk, function(j) Q = Qf[, -(coef.no[j])])
  }
  
  getBstats <- function(coef.no, Beta, Hb, Qr, Y, n, p, ind){
    k <- ncol(Qf)
    Fitted <- fastFit(Qr, Y, n, p)
    Resid <- Y - Fitted
    Result <- sapply(ind, function(x) {
      yp <- Fitted + Resid[x,]
      B <- as.matrix(Hb %*% yp)
      R <- yp - fastFit(Qf, yp, n, p)
      S <- fast.solve(crossprod(R) / (n - k))
      b <- B[coef.no,]
      d <- sqrt(sum(b^2))
      md <- sqrt(crossprod(b, S) %*% b)
      c(d = d, md = md)
      })
    B <- as.matrix(Hb %*% Y)
    R <- Y - fastFit(Qf, Y, n, p)
    S <- fast.solve(crossprod(R) / (n - k))
    b <- B[coef.no,] - Beta
    d <- sqrt(sum(b^2))
    md <- sqrt(crossprod(b, S) %*% b)
    Result[, 1] <- c(d = d, md = md)
    colnames(Result) <- names(ind)
    Result
  }
  
  PI <- getPermInfo(fit, "all")
  ind <- PI$perm.schedule
  perms <- length(ind)

  
  Result <- lapply(1:kk, function(j){
    if(print.progress){
      cat("\nCalculating statistics for", coef.nms[j], ",",
          perms, "permutations...")
    }
    cn <- coef.no[j]
    qr <- Q.list[[j]]
    getBstats(cn, Beta, Hb, qr, TY, n, p, ind)
  })
  
  names(Result) <- coef.nms
  
  obs.d <- sapply(Result, function(x) x[, 1][1])
  obs.md <- sapply(Result, function(x) x[, 2][2])
  
  out <- list(
    obs.d = obs.d,
    obs.md = obs.md, 
    Beta = Beta,
    obs.B.mat = Bobs,
    coef.no = coef.no,
    random.stats = Result          
  )
  
    attr(out, "class") <- "betaTest"
    out
}
