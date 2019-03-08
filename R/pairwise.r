#' Pairwise comparisons of lm.rrpp fits
#'
#' Function generates distributions of pairwise statistics for a lm.rrpp fit and
#' returns important statistics for hypothesis tests.
#'
#' Based on an lm.rrpp fit, this function will find fitted values over all permutations and based 
#' on a grouping factor, calculate either least squares (LS) means or slopes, and pairwise statistics among
#' them.  Pairwise statistics have two flavors: distances and vector correlations (or angles).  The distance
#' statistics calculate either the length of vectors between LS mean vectors or the absolute difference between 
#' slope vector lengths.  The vector correlations are the inner product of vectors that have been transformed to unit length.
#' The arccosine (acos) of this value is the angle between vectors, which can be expressed in radians or degrees, and is
#' used as a test statistic (with the null hypothesis that vectors are parallel; angle = 0).
#' Over all permutations, these values can be calculated to generate random distributions using the null model.  The 
#' null model is defined via \code{\link{lm.rrpp}}, but one can also use an alternative null model as an optional argument.
#' In this case, residual randomization in the permutation procedure (RRPP) will be performed using the alternative null model 
#' to generate fitted values.  If full randomization of values (FRPP) is preferred,
#' it must be established in the lm.rrpp fit and an alternative model should not be chosen. 
#' 
#' Observed statistics, effect sizes, P-values, and one-tailed confidence limits based on the confidence requested will
#' be summarized with the \code{\link{summary.pairwise}} function.  The \code{\link{summary.pairwise}} function will allow one
#' to select between distance or vector correlation tests, whether angles are measured in radians or degrees, and the level of
#' confidence for the test.  Confidence limits are inherently one-tailed as
#' the statistics are similar to absolute values.  For example, a distance is analogous to an absolute difference.  Therefore,
#' the one-tailed confidence limits are more akin to two-tailed hypothesis tests.  (A comparable example is to use the absolute 
#' value of a t-statistic, in which case the distribution has a lower bound of 0.)  If rather than comparing the LS means or slopes,
#' one wishes to compare the dispersion of residuals among groups, given the model, an option for comparing variances is also
#' available.  Variance degrees of freedom equal n, the group size, rather than n-1, as the purpose is to compare mean dispersion 
#' in the sample.  (Additionally, tests with one subject in a group are possible, or at least not a hindrance to the analysis.)
#' 
#' If data are univariate, test.type = 'cor' should not be chosen because the vector correlation between univariate 
#' vectors is always 1.  Rather, cor.type = 'dist' will return the absolute difference between slopes or between means.  
#' Please note that this function will generate results if test.type = 'cor' for univariate data, but the results will 
#' not make much sense.
#' 
#' @param fit A linear model fit using \code{\link{lm.rrpp}}.
#' @param fit.null An alternative linear model fit to use as a null model for RRPP, if the null model
#' of the fit is not desired.  Note, for FRPP this argument should remain NULL and FRPP
#' must be established in the lm.rrpp fit (RRPP = FALSE).  If the null model is uncertain, 
#' using \code{\link{reveal.model.designs}} will help elucidate the inherent null model used.
#' @param groups A factor or vector that is coercible into a factor, describing the levels of
#' the groups for which to find LS means or slopes.  Normally this factor would be part of the 
#' model fit, but it is not necessary for that to be the case in order to obtain results.
#' @param covariate A numeric vector for which to calculate slopes for comparison  If NULL, 
#' LS means will be calculated instead of slopes.  Normally this variable would be part of the 
#' model fit, but it is not necessary for that to be the case in order to obtain results.
#' @param print.progress If a null model fit is provided, a logical value to indicate whether analytical 
#' results progress should be printed on screen.  Unless large data sets are analyzed, this argument 
#' is probably not helpful.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{pairwise} is a list containing the following
#' \item{LS.means}{LS means for groups, across permutations.}
#' \item{slopes}{Slopes for groups, across permutations.}
#' \item{means.dist}{Pairwise distances between means, across permutations.}
#' \item{means.vec.cor}{Pairwise vector correlations between means, across permutations.}
#' \item{slopes.lengths}{Slope lengths, by group, across permutations.}
#' \item{slopes.dist}{Pairwise distances between slope lengths, across permutations.}
#' \item{slopes.vec.cor}{Pairwise vector correlations between slope vectors, across permutations.}
#' \item{n}{Sample size}
#' \item{p}{Data dimensions; i.e., variable number}
#' \item{PermInfo}{Information for random permutations, passed on from lm.rrpp fit and possibly
#' modified if an alternative null model was used.}
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C and M.L. Collyer. 2018. Multivariate phylogenetic anova: group-clade aggregation, biological 
#' challenges, and a refined permutation procedure. Evolution. In press.
#' @seealso \code{advanced.procD.lm} within \code{geomorph}; \code{\link{lm.rrpp}} for model fits
#' @examples 
#' 
#' # Examples use geometric morphometric data on pupfishes
#' # See the package, geomorph, for details about obtaining such data
#'
#' # Body Shape Analysis (Multivariate)----------------------------------------------------
#' 
#' data("Pupfish")
#' 
#' # Note:
#' 
#' dim(Pupfish$coords) # highly multivariate!
#' 
#' Pupfish$logSize <- log(Pupfish$CS) # better to not have functions in formulas
#'
#' # Note: one should use all dimensions of the data but with this example, there are many
#' # Thus, only three principal components will be used for demonstration purposes.
#' 
#' Pupfish$Y <- prcomp(Pupfish$coords)$x[, 1:3]
#' 
#' ## Pairwise comparisons of LS means
#' 
#' fit1 <- lm.rrpp(Y ~ logSize + Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 999) 
#' summary(fit1, formula = FALSE)
#' anova(fit1) 
#'
#' pup.group <- interaction(Pupfish$Sex, Pupfish$Pop)
#' pup.group
#' PW1 <- pairwise(fit1, groups = pup.group)
#' PW1
#' summary(PW1, confidence = 0.95, test.type = "dist") # distances between means
#' summary(PW1, confidence = 0.95, test.type = "dist", stat.table = FALSE)
#' summary(PW1, confidence = 0.95, test.type = "VC", 
#'    angle.type = "deg") # correlation between mean vectors (angles in degrees)
#'
#' # Can also compare the dispersion around means
#' 
#' summary(PW1, confidence = 0.95, test.type = "var")
#' 
#' ## Pairwise comparisons of slopes
#' 
#' fit2 <- lm.rrpp(Y ~ logSize * Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 999) 
#' summary(fit2, formula = FALSE)
#' anova(fit1, fit2)
#'
#' # Using a null fit that excludes all factor-covariate interactions, not just the last one  
#' 
#' PW2 <- pairwise(fit2, fit.null = fit1, groups = pup.group, covariate = Pupfish$logSize) 
#' PW2
#' summary(PW2, confidence = 0.95, test.type = "dist") # distances between slope vector lengths
#' summary(PW2, confidence = 0.95, test.type = "dist", stat.table = FALSE)
#' summary(PW2, confidence = 0.95, test.type = "VC",
#'    angle.type = "deg") # correlation between slope vectors (and angles)
#'    
#' # Can also compare the dispersion around group slopes
#' 
#' summary(PW2, confidence = 0.95, test.type = "var")
#' 
pairwise <- function(fit, fit.null = NULL, groups, covariate = NULL, 
                      print.progress = FALSE) {
  fitf <- fit
  ind <- fit$PermInfo$perm.schedule
  perms <- length(ind)
  gls <- fit$LM$gls
  if(!inherits(fit, "lm.rrpp")) stop("The model fit must be a lm.rrpp class object")
  
  if(!is.null(fit.null)) {
    fitf$PermInfo$perm.method <- "RRPP"
    if(!inherits(fit.null, "lm.rrpp")) stop("The null model fit must be a lm.rrpp class object")
    
    if(print.progress){
      cat(paste("\nCoefficients estimation from RRPP on null model:", perms, "permutations.\n"))
      pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
    }
    
    gls <- fit.null$LM$gls
    rrpp.args <- list()
    fitted <- fit.null$LM$wFitted
    res <- fit.null$LM$wResiduals
    n <- fit.null$LM$n
    p <- fit.null$LM$p
    k <- length(fit.null$LM$term.labels)
    o <- fit.null$LM$offset
    if(sum(o) != 0) offset = TRUE else offset = FALSE
    rrpp.args <- list(fitted = list(fitted), residuals = list(res),
                      ind.i = NULL, w = NULL, o = NULL)
    Y <- fit.null$LM$Y * sqrt(fit$LM$weights)
    X0 <- fit.null$LM$X * sqrt(fit$LM$weights)
    Q <- qr(X0)
    U <- qr.Q(Q)
    
    if(gls) {
      P <- fit.null$LM$Pcov
      Y <- crossprod(P, Y)
      X0 <- crossprod(P, X0)
      Q <- qr(X0)
      U <- qr.Q(Q)
      fitted <- fastFit(U, Y, n, p)
      res <- Y - fitted
      rrpp.args$fitted <- list(fitted)
      rrpp.args$residuals <- list(res)
    }
    
    Xf <- fitf$LM$X * sqrt(fitf$LM$weights)
    if(fitf$LM$gls) Xf <- crossprod(fitf$LM$Pcov, Xf)
    Qf <- qr(Xf)
    H <- tcrossprod(solve(qr.R(Qf)), qr.Q(Qf))
    getCoef <- function(y) H %*% y
    coef.n <- lapply(1:perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      rrpp.args$ind.i <- ind[[j]]
      y <- do.call(rrpp, rrpp.args)[[1]]
      getCoef(y)
    })
    
    kf <- length(fitf$LM$term.labels)
    fitf$LM$random.coef[[kf]] <- coef.n
    step <- perms + 1
    if(print.progress) {
      setTxtProgressBar(pb,step)
      close(pb)
    }
  }
  
  groups <- as.factor(groups)
  gp.rep <- by(groups, groups, length)
  
  if(!all(gp.rep > 1)) stop("The groups factor does not have replication at each level.")
  if(length(groups) != fitf$LM$n) 
    stop("The length of the groups factor does not match the number of observations in the lm.rrpp fit")
  if(!is.null(covariate)) {
    if(!is.numeric(covariate)) stop("The covariate argument must be numeric")
    if(inherits(covariate, "matrix")) stop("The covariate argument can only be a single covariate, not a matrix")
    if(length(unique(covariate)) < 2) stop("The covariate vector appears to have no variation.")
    if(length(covariate) != fitf$LM$n)
      stop("The length of the covariate does not match the number of observations in the lm.rrpp fit")
  }
  if(is.null(covariate)) {
    means <- getLSmeans(fit = fitf, g = groups)
    slopes <- NULL
    means.dist <- lapply(means, function(x) as.matrix(dist(x)))
    means.vec.cor <- lapply(means, vec.cor.matrix)
    slopes.length <- NULL
    slopes.dist <- NULL
    slopes.vec.cor <- NULL
  } else {
    means <- NULL
    slopes <- getSlopes(fit = fitf, x = covariate, g = groups)
    means.dist <- NULL
    means.vec.cor <- NULL
    slopes.length <- lapply(slopes, function(x) sqrt(rowSums(x^2)))
    slopes.dist <- lapply(slopes.length, function(x) as.matrix(dist(as.matrix(x))))
    slopes.vec.cor <- lapply(slopes, vec.cor.matrix)
  }
  
  if(gls) res <- fitf$LM$gls.residuals else res <- fitf$LM$wResiduals
  disp.args <- list(res = as.matrix(res), ind.i = NULL, x = model.matrix(~groups + 0))
  if(gls) disp.args$x <- crossprod(fitf$LM$Pcov, disp.args$x)
  g.disp <- function(res, ind.i, x) {
    r <- res[ind.i,]
    if(NCOL(r) > 1) d <- apply(r, 1, function(x) sum(x^2)) else
      d <- r^2
    coef(lm.fit(x, d))
  }

  vars <- sapply(1:perms, function(j){
    disp.args$ind.i <- ind[[j]]
    do.call(g.disp, disp.args)
  })
  
  rownames(vars) <- levels(groups)
  colnames(vars) <- c("obs", paste("iter", 1:(perms -1), sep = "."))

  out <- list(LS.means = means, slopes = slopes, means.dist = means.dist, means.vec.cor = means.vec.cor,
              slopes.length = slopes.length, slopes.dist = slopes.dist, slopes.vec.cor = slopes.vec.cor,
              vars = vars, n = fit$LM$n, p = fit$LM$p, PermInfo = fitf$PermInfo)
  
  class(out) <- "pairwise"
  out
}
