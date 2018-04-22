#' Linear Model Evaluation with a Randomized Residual Permutation Procedure
#'
#' Function performs a linear model fit over many random permutations of data, using
#' a randomized residual permutation procedure.
#'
#' The function fits a linear model using ordinary least squares (OLS) or generalized
#' least squares (GLS) estimation of coefficients over any number of random permutations of
#' the data.  A permutation procedure that randomizes vectors of residuals is employed.  This
#' procedure can randomize two types of residuals: residuals from null models or residuals from
#' an intercept model.  The latter is the same as randomizing full values, and is referred to as
#' as a full randomization permutation procedure (FRPP); the former uses the residuals from null
#' models, which are defined by the type of sums of squares and cross-products (SSCP) sought in an
#' analysis of variance (ANOVA), and is referred to as a randomized residual permutation procedure (RRPP).
#' Types I, II, and III SSCPs are supported.
#'
#' Users define the SSCP type, the permutation procedure type, whether a covariance matrix is included
#' (GLS estimation), and a few arguments related to computations.  Analytical results comprise observed linear model
#' results (coefficients, fitted values, residuals, etc.), random sums of squares (SS) across permutation iterations,
#' and other parameters for performing ANOVA and other hypothesis tests, using
#' empirically-derived probability distributions.
#'
#' \code{lm.rrpp} emphasizes estimation of standard deviates of observed statistics as effect sizes
#' from distributions of random outcomes.  When performing ANOVA, using the \code{\link{anova}} function,
#' the effect type (statistic choice) can be varied.  See \code{\link{anova.lm.rrpp}} for more details.  Please
#' recognize that the type of SS must be chosen prior to running \code{lm.rrpp} and not when applying \code{\link{anova}}
#' to the \code{lm.rrpp} fit, as design matrices for the linear model must be created first.  Therefore, SS.type
#' is an argument for \code{lm.rrpp} and effect.type is an argument for \code{\link{anova.lm.rrpp}}.
#' 
#' The \code{\link{coef.lm.rrpp}} function can be used to test the specific coefficients of an lm.rrpp fit.  The test
#' statistics are the distances (d), which are also standardized (Z-scores).  The Z-scores might be easier to compare,
#' as the expected values for random distances can vary among coefficient vectors (Adams and Collyer 2016).
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2).  Can also be a linear model fit
#' from \code{\link{lm}}.
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing
#' @param SS.type A choice between type I (sequential), type II (hierarchical), or type III (marginal)
#' sums of squares and cross-products computations.
#' @param Cov An optional argument for including a covariance matrix to address the non-independence
#' of error in the estimation of coefficients (via GLS).
#' @param data A data frame for the function environment, see \code{\link{rrpp.data.frame}}
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.
#' This is helpful for long-running analyses.
#' @param Parallel A logical value to indicate whether parallel processing should be used.  If TRUE, this argument
#' invokes forking of processor cores, using the \code{parallel} library.  This option is only available to unix systems
#' and should only be used for rather long analyses (that would normally take over 10 seconds on a single core).  Currently,
#' parallel processing is performed on all but one core with no option to change the number of cores.  Systems with Windows
#' platforms will automatically default to a single-core application of this function.
#' @param ... Arguments typically used in \code{\link{lm}}, such as weights or offset, passed on to
#' \code{rrpp.fit} for estimation of coefficients.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{lm.rrpp} is a list containing the following
#' \item{call}{The matched call.}
#' \item{LM}{Linear Model objects, including data (Y), coefficients, design matrix (X), sample size
#' (n), number of dependent variables (p), dimension of data space (p.prime),
#' QR decomposition of the design matrix, fitted values, residuals,
#' weights, offset, model terms, model frame (data), random coefficients (through permutations),
#' random vector distances for coefficients (through permutations), and whether OLS or GLS was performed}
#' \item{ANOVA}{Analysis of variance objects, including the SS type, random SS outcomes, random MS outcomes,
#' random R-squared outcomes, random F outcomes, random Cohen's f-squared outcomes, P-values based on random F
#' outcomes, effect sizes for random outcomes, sample size (n), number of variables (p), and degrees of freedom for
#' model terms (df).  These objects are used to construct ANOVA tables.}
#' \item{PermInfo}{Permutation procedure information, including the number of permutations (perms), The method
#' of residual randomization (perm.method), and each permutation's sampling frame (perm.schedule), which
#' is a list of reordered sequences of 1:n, for how residuals were randomized.}
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance.
#'    Austral Ecology 26: 32-46.
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Computation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C and M.L. Collyer. 2018. Multivariate phylogenetic anova: group-clade aggregation, biological 
#' challenges, and a refined permutation procedure. Evolution. In press.
#' @seealso \code{procD.lm} and \code{procD.pgls} within \code{geomorph}; \code{\link[stats]{lm}} for more on linear model fits.
#' @examples
#' 
#' # Examples use geometric morphometric data
#' # See the package, geomorph, for details about obtaining such data
#'
#' data("PupfishHeads")
#' names(PupfishHeads)
#' 
#' # Head Size Analysis (Univariate)-------------------------------------------------------
#'
#' # Note: lm.rrpp works best if one avoids functions within formulas
#' # Thus,
#' 
#' PupfishHeads$logHeadSize <- log(PupfishHeads$headSize)
#' names(PupfishHeads)
#'
#' fit <- lm.rrpp(logHeadSize ~ sex + locality/year, SS.type = "I", data = PupfishHeads)
#' summary(fit)
#' anova(fit, effect.type = "F") # Maybe not most appropriate
#' anova(fit, effect.type = "Rsq") # Change effect type, but still not most appropriate
#'
#' # Mixed-model approach (most appropriate, as year sampled is a random effect:
#' 
#' anova(fit, effect.type = "F", error = c("Residuals", "locality:year", "Residuals"))
#'
#' # Change to Type III SS
#' 
#' fit <- lm.rrpp(logHeadSize ~ sex + locality/year, SS.type = "III", data = PupfishHeads)
#' summary(fit)
#' anova(fit, effect.type = "F", error = c("Residuals", "locality:year", "Residuals"))
#'
#' # Coefficients Test
#' 
#' coef(fit, test = TRUE)
#'
#' # Predictions (holding alternative effects constant)
#' 
#' sizeDF <- data.frame(sex = c("Female", "Male"))
#' rownames(sizeDF) <- c("Female", "Male")
#' sizePreds <- predict(fit, sizeDF)
#' summary(sizePreds)
#' plot(sizePreds)
#' 
#' # Diagnostics plots of residuals
#' 
#' plot(fit)
#' 
#' # Body Shape Analysis (Multivariate)----------------------------------------------------
#' 
#' data(Pupfish)
#' names(Pupfish)
#' 
#' # Note:
#' 
#' dim(Pupfish$coords) # highly multivariate!
#' 
#' Pupfish$logSize <- log(Pupfish$CS) # better to not have functions in formulas
#' names(Pupfish)
#'
#' # Note: one should increase RRPP iterations but they are not used at all
#' # here for a fast example.  Generally, iter = 999 will take less
#' # than 1s for this example with a modern computer.
#' 
#' fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 0) 
#' summary(fit, formula = FALSE)
#' anova(fit) 
#' coef(fit, test = TRUE)
#'
#' # Predictions (holding alternative effects constant)
#' 
#' shapeDF <- expand.grid(Sex = levels(Pupfish$Sex), Pop = levels(Pupfish$Pop))
#' rownames(shapeDF) <- paste(shapeDF$Sex, shapeDF$Pop, sep = ".")
#' shapeDF
#' 
#' shapePreds <- predict(fit, shapeDF)
#' summary(shapePreds)
#' summary(shapePreds, PC = TRUE)
#' 
#' # Plot prediction
#' 
#' plot(shapePreds, PC = TRUE)
#' plot(shapePreds, PC = TRUE, ellipse = TRUE)
#' 
#' # Diagnostics plots of residuals
#' 
#' plot(fit)
#' 
#' # PC-plot of fitted values
#' 
#' groups <- interaction(Pupfish$Sex, Pupfish$Pop)
#' plot(fit, type = "PC", pch = 19, col = as.numeric(groups))
#' 
#' # Regression-like plot
#' 
#' plot(fit, type = "regression", reg.type = "PredLine", 
#'     predictor = Pupfish$logSize, pch=19,
#'     col = as.numeric(groups))
#'
#' # Body Shape Analysis (Distances)----------------------------------------------------
#' 
#' D <- dist(Pupfish$coords) # inter-observation distances
#' length(D)
#' Pupfish$D <- D
#' 
#' # Note: one should increase RRPP iterations but they are not used at all
#' # here for a fast example.  Generally, iter = 999 will take less
#' # than 1s for this example with a modern computer.
#' 
#' fitD <- lm.rrpp(D ~ logSize + Sex*Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 0) 
#' 
#' # These should be the same:
#' summary(fitD, formula = FALSE)
#' summary(fit, formula = FALSE) 
#'
#' # GLS Example (Univariate) ----------------------------------------------------------
#' 
#' data(PlethMorph)
#' fitOLS <- lm.rrpp(TailLength ~ SVL, data = PlethMorph)
#' fitGLS <- lm.rrpp(TailLength ~ SVL, data = PlethMorph, Cov = PlethMorph$PhyCov)
#' 
#' anova(fitOLS)
#' anova(fitGLS)
#' 
#' sizeDF <- data.frame(SVL = sort(PlethMorph$SVL))
#' plot(predict(fitOLS, sizeDF)) # Correlated error
#' plot(predict(fitGLS, sizeDF)) # Independent error
#' 
#' #' # GLS Example (Multivariate) ----------------------------------------------------------
#' 
#' Y <- as.matrix(cbind(PlethMorph$TailLength,
#' PlethMorph$HeadLength,
#' PlethMorph$TailLength,
#' PlethMorph$Snout.eye,
#' PlethMorph$BodyWidth,
#' PlethMorph$Forelimb,
#' PlethMorph$Hindlimb))
#' PlethMorph <- rrpp.data.frame(PlethMorph, Y=Y)
#' fitOLSm <- lm.rrpp(Y ~ SVL, data = PlethMorph)
#' fitGLSm <- lm.rrpp(Y ~ SVL, data = PlethMorph, Cov = PlethMorph$PhyCov)
#' 
#' anova(fitOLSm)
#' anova(fitGLSm)
#' 
#' plot(predict(fitOLSm, sizeDF), PC= TRUE) # Correlated error
#' plot(predict(fitGLSm, sizeDF), PC= TRUE) # Independent error

lm.rrpp <- function(f1, iter = 999, seed = NULL, int.first = FALSE,
                    RRPP = TRUE, SS.type = c("I", "II", "III"),
                    data = NULL, Cov = NULL,
                    print.progress = TRUE, Parallel = FALSE, ...) {
  if(int.first) ko = TRUE else ko = FALSE
  if(!is.null(data) && inherits(data, "data.frame")) data <- rrpp.data.frame(data)
  if(!is.null(data)) data <- droplevels.rrpp.data.frame(data)
  if(print.progress){
    cat("\nPreliminary Model Fit...\n")
    pb <- txtProgressBar(min = 0, max = 4, initial = 0, style=3)
    step <- 1
  }
  SS.type <- match.arg(SS.type)
  form <- formula(f1)
  dep <- eval(form[[2]], data, parent.frame())
  if(inherits(dep, "dist")) {
    if(any(dep < 0)) stop("Distances in distance matrix cannot be less than 0")
    D <- dep
  } else if((is.matrix(dep) || is.data.frame(dep))
            && isSymmetric(unname(as.matrix(dep)))) {
    D <- as.dist(dep)
  } else D <- NULL
  if(print.progress) {
    step <- 2
    setTxtProgressBar(pb,step)
  }
  fit.o <- rrpp.fit(f1, data = data, keep.order = ko, pca = FALSE, SS.type = SS.type,...)
  if(print.progress) {
    step <- 3
    setTxtProgressBar(pb,step)
  }
  dimsY <- dim(as.matrix(fit.o$Y))
  n <- dimsY[1]; p <- p.prime <- dimsY[2]
  if(p > n){
    fit <- rrpp.fit(f1, data = data, keep.order = ko, pca = TRUE, SS.type = SS.type, ...)
    Y <- as.matrix(fit$Y)
    dimsY <- dim(Y)
    n <- dimsY[1]; p.prime <- dimsY[2]
  } else {
    fit <- fit.o
    Y <- as.matrix(fit$Y)
  }
  if(print.progress) {
    step <- 3
    setTxtProgressBar(pb,step)
  }
  k <- length(fit$term.labels)
  id <- rownames(Y)
  ind <- perm.index(n, iter = iter, seed = NULL)
  perms <- iter + 1
  if(print.progress) {
    step <- 4
    setTxtProgressBar(pb,step)
    close(pb)
  }
  SS.args <- list(fit = fit, ind = ind, P = NULL,
                  RRPP = RRPP, print.progress = print.progress)
  beta.args <- SS.args
  beta.args$fit <- fit.o
  if(!is.null(Cov)){
    Cov.name <- deparse(substitute(Cov))
    Cov.match <- match(Cov.name, names(data))
    if(length(Cov.match) > 1) stop("More than one object matches covariance matrix name")
    if(all(is.na(Cov.match))) Cov <- Cov else Cov <- data[[Cov.match]]
    if(!is.matrix(Cov)) stop("The covariance matrix must be a matrix.")
    dimsC <- dim(Cov)
    if(!all(dimsC == n))
      stop("Either one or both of the dimensions of the covariance matrix do not match the number of observations.")
    Pcov <- Cov.proj(Cov, id)
    SS.args$P <- Pcov
  } else Pcov <- NULL 
  if(k > 0){
    if(Parallel) {
      if(.Platform$OS.type == "windows") betas <- do.call(beta.iter, beta.args)
      else betas <- do.call(beta.iterPP, beta.args)
    } else betas <- do.call(beta.iter, beta.args)
    if(Parallel) {
      if(.Platform$OS.type == "windows") SS <- do.call(SS.iter, SS.args)
      else SS <- do.call(SS.iterPP, SS.args)
      } else SS <- do.call(SS.iter, SS.args)
    ANOVA <- anova.parts(fit, SS)
    LM <- list(coefficients = fit.o$wCoefficients.full[[k]],
               Y = fit.o$Y,  X = fit.o$X, n = n, p = p, p.prime = p.prime,
               QR = fit.o$QRs.full[[k]],
               fitted = fit.o$fitted.full[[k]],
               residuals = fit.o$residuals.full[[k]],
               weights = fit.o$weights, offset = fit.o$offset,
               wQR = fit.o$wQRs.full[[k]],
               wFitted = fit.o$wFitted.full[[k]],
               wResiduals = fit.o$wResiduals.full[[k]],
               Terms = fit.o$Terms, term.labels = fit.o$term.labels,
               data = fit.o$data, 
               random.coef = betas$random.coef,
               random.coef.distances = betas$random.coef.distances,
               ols = TRUE, gls = FALSE)
    PermInfo <- list(perms = perms,
                     perm.method = ifelse(RRPP==TRUE,"RRPP", "FRPP"), 
                     perm.schedule = ind, perm.seed = seed)
    out <- list(call = match.call(), LM = LM, ANOVA = ANOVA, PermInfo = PermInfo)
    if(!is.null(Cov)) {
      out$LM$Cov <- Cov
      out$LM$Pcov <- Pcov
      PY <- crossprod(Pcov, fit.o$wY); PX <- crossprod(Pcov, fit.o$wX)
      w <- fit.o$weights
      fit.cov <- lm.fit(PX, PY)
      out$LM$gls = TRUE; out$LM$ols = FALSE
      out$LM$gls.coefficients = fit.cov$coefficients
      out$LM$gls.fitted = fit.o$X%*%fit.cov$coefficients
      out$LM$gls.residuals = fit.o$Y - fit.o$X%*%fit.cov$coefficients
      out$LM$gls.mean <- colMeans(fit.o$X%*%fit.cov$coefficients) # Fix this
    }
  } else
  {
    if(print.progress) cat("\n No terms for ANOVA; only RSS calculated in each permutation\n")
    if(!is.null(Cov)){
      betas <- beta.iter.null(fit,  P = Pcov, ind = ind,  
                              RRPP = RRPP, print.progress = print.progress)
      SS <- SS.iter.null(fit,  P = Pcov, ind = ind,  
                         RRPP = RRPP, print.progress = print.progress)
    }  else {
      betas <- beta.iter.null(fit, ind = ind, RRPP=RRPP, print.progress = print.progress)
      SS <- SS.iter.null(fit, ind = ind, RRPP=RRPP, print.progress = print.progress)
    }
    SSY <- SS[1]
    df <- n - fit$wQRs.full[[1]]$rank
    LM <- list(coefficients=fit$wCoefficients.full[[1]],
               Y=fit$Y,  X=fit$X, n = n, p = p, p = p.prime,
               QR = fit.o$QRs.full[[1]], fitted = fit$fitted.full[[1]],
               residuals = fit$residuals.full[[1]],
               weights = fit$weights, offset = fit$offset,
               wQR = fit.o$wQRs.full[[1]],
               wFitted = fit.o$wFitted.full[[1]],
               wResiduals = fit.o$wResiduals.full[[1]],
               Terms = fit$Terms,
               term.labels = fit$term.labels, 
               data = fit.o$data, 
               random.coef = betas$random.coef,
               random.coef.distances = betas$random.coef.distances,
               ols = TRUE, gls = FALSE)
    ANOVA <- list(df = df, SS = SS, MS = SS/df)
    PermInfo <- list(perms = perms,
                     perm.method = ifelse(RRPP==TRUE,"RRPP", "FRPP"), 
                     perm.schedule = ind, perm.seed = seed)
    out <- list(call = match.call(), LM = LM, ANOVA = ANOVA, PermInfo = PermInfo)
    if(!is.null(Cov)){
      out$LM$Cov <- Cov
      out$LM$Pcov <- Pcov
      PY <- Pcov%*%fit.o$wY; PX <- Pcov%*%fit.o$wX
      fit.cov <- lm.fit(PX, PY)
      out$LM$gls = TRUE; out$LM$ols = FALSE
      out$LM$gls.coefficients = fit.cov$coefficients
      out$LM$gls.fitted = fit.o$X%*%fit.cov$coefficients
      out$LM$gls.residuals = fit.o$Y - fit.o$X%*%fit.cov$coefficients
      out$LM$gls.mean <- colMeans(fit.o$X%*%fit.cov$coefficients)
    }
  }
  if(!is.null(D)) {
    qrf <- fit$QRs.full[[k]]
    wqrf <- fit$wQRs.full[[k]]
    D.coef <- qr.coef(qrf, D)
    D.wcoef <- qr.coef(wqrf, D)
    out$LM$dist.coefficients <- D.coef
    out$LM$dist.wCoefficients <- D.wcoef
  }
  class(out) = "lm.rrpp"
  out
}

