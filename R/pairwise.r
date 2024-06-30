#' Pairwise comparisons of lm.rrpp fits
#'
#' Function generates distributions of pairwise statistics for a lm.rrpp fit and
#' returns important statistics for hypothesis tests.
#'
#' Based on an lm.rrpp fit, this function will find fitted values over all 
#' permutations and based 
#' on a grouping factor, calculate either least squares (LS) means or 
#' slopes, and pairwise statistics among
#' them.  Pairwise statistics have multiple flavors, related to vector attributes: 
#' 
#' \describe{
#' \item{\bold{Distance between vectors, "dist"}}{ Vectors for LS means or 
#' slopes originate at the origin and point to some location, having both a 
#' magnitude
#' and direction.  A distance between two vectors is the inner-product of of 
#' the vector difference, i.e., the distance between their endpoints.  For
#' LS means, this distance is the difference between means.  For multivariate 
#' slope vectors, this is the difference in location between estimated change 
#' for the dependent variables, per one-unit change of the covariate considered.
#' For univariate slopes, this is the absolute difference between slopes.}
#' \item{\bold{Standardized distance between vectors, "stdist"}}{Same as the 
#' distance between vectors, but distances are divided by the model standard error 
#' (square-root of the trace of the residual covariance matrix, adjusted by sample size). 
#' Pairwise tests with this statistic should be consistent with ANOVA results.}
#' \item{\bold{Mahalanobis distance between vectors, "mdist"}}{Similar to the
#' standardized distance between vectors but the inverse of the 
#' residual covariance matrix is used in calculation of the distance, rather than
#' dividing the Euclidean distance between means and dividing by the standard error. 
#' If the residual covariance matrix is singular, Mahalanobis distances will not be 
#' calculated.  Pairwise tests with this statistic should be consistent with 
#' MANOVA results.}
#' \item{\bold{Vector correlation, "VC"}}{ If LS mean or slope vectors are 
#' scaled to unit size, the vector correlation is the inner-product of the 
#' scaled vectors.
#' The arccosine (acos) of this value is the angle between vectors, which can 
#' be expressed in radians or degrees.  Vector correlation indicates the 
#' similarity of 
#' vector orientation, independent of vector length.}
#' \item{\bold{Difference in vector lengths, "DL"}}{  If the length of a vector 
#' is an important attribute -- e.g., the amount of multivariate change per 
#' one-unit
#' change in a covariate -- then the absolute value of the difference in 
#' vector lengths is a practical statistic to compare vector lengths.  Let 
#' d1 and
#' d2 be the distances (length) of vectors.  Then |d1 - d2| is a statistic 
#' that compares their lengths.  For slope vectors, this is a comparison of rates.}
#' \item{\bold{Variance, "var}}{  Vectors of residuals from a linear model 
#' indicate can express the distances of observed values from fitted values.  
#' Mean
#' squared distances of values (variance), by group, can be used to measure 
#' the amount of dispersion around estimated values for groups.  Absolute
#' differences between variances are used as test statistics to compare mean 
#' dispersion of values among groups.  Variance degrees of freedom equal n, 
#' the group size, rather than n-1, as the purpose is to compare mean dispersion 
#' in the sample.  (Additionally, tests with one subject in a group are 
#' possible, or at least not a hindrance to the analysis.)}
#' }
#' 
#' The \code{\link{summary.pairwise}} function is used to select a test 
#' statistic for the statistics described above, as
#' "dist", "VC", "DL", and "var", respectively.  If vector correlation is tested, 
#' the \code{angle.type} argument can be used to choose between radians and
#' degrees.
#' 
#' The null model is defined via \code{\link{lm.rrpp}}, but one can 
#' also use an alternative null model as an optional argument.
#' In this case, residual randomization in the permutation procedure 
#' (RRPP) will be performed using the alternative null model 
#' to generate fitted values.  If full randomization of values (FRPP) 
#' is preferred,
#' it must be established in the lm.rrpp fit and an alternative model 
#' should not be chosen. If one is unsure about the inherent
#' null model used if an alternative is not specified as an argument, 
#' the function \code{\link{reveal.model.designs}} can be used.
#' 
#' Observed statistics, effect sizes, P-values, and one-tailed confidence 
#' limits based on the confidence requested will
#' be summarized with the \code{\link{summary.pairwise}} function.  
#' Confidence limits are inherently one-tailed as
#' the statistics are similar to absolute values.  For example, a 
#' distance is analogous to an absolute difference.  Therefore,
#' the one-tailed confidence limits are more akin to two-tailed 
#' hypothesis tests.  (A comparable example is to use the absolute 
#' value of a t-statistic, in which case the distribution has a lower 
#' bound of 0.)  
#' 
#'  \subsection{Notes for RRPP 0.6.2 and subsequent versions}{ 
#'  In previous versions of pairwise, \code{\link{summary.pairwise}} had three 
#'  test types: "dist", "VC", and "var".  When one chose "dist", for LS mean 
#'  vectors, the statistic was the inner-product of the vector difference.  
#'  For slope vectors, "dist" returned the absolute value  of the difference 
#'  between vector lengths, which is "DL" in 0.6.2 and subsequent versions.  This
#'  update uses the same calculation, irrespective of vector types.  Generally,
#'  "DL" is the same as a contrast in rates for slope vectors, but might not have
#'  much meaning for LS means.  Likewise, "dist" is the distance between vector
#'  endpoints, which might make more sense for LS means than slope vectors.  
#'  Nevertheless, the user has more control over these decisions with version 0.6.2
#'  and subsequent versions.
#' }
#' 
#'  \subsection{Notes for RRPP 2.0.4 and subsequent versions}{ 
#'  The test types, standardized distance between vectors, "stdist", and Mahalanobis distances
#'  between vectors were added.  The former
#'  simply divides the distance between vectors by model standard error (square-root 
#'  of the trace of the residual covariance matrix, adjusted by sample size).  
#'  This is a multivariate generalization
#'  of a t-statistic for the difference between means.  It is not the same as Hotelling 
#'  squared-T-statistic, which requires incorporating the inverse of the residual covariance 
#'  matrix in the calculation of the distance, a calculation that also requires a non-singular
#'  covariance matrix.  However, the Mahalanobis distances are similar (and proportional) 
#'  to the Hotelling squared-T-statistic. Pairwise tests using Mahalanobis distances represent
#'  a non-parametric analog to the parametric Hotelling squared-T test.  Both tests should be 
#'  better for GLS model fits compared
#'  to Euclidean distances between means, as the total sums of squares are more likely to vary
#'  across random permutations. In general, if ANOVA is 
#'  performed a pairwise test with "stdist" is a good association; if MANOVA is performed,
#'  a pairwise test with "mdist" is a good association.
#' }
#'  
#' @param fit A linear model fit using \code{\link{lm.rrpp}}.
#' @param fit.null An alternative linear model fit to use as a null 
#' model for RRPP, if the null model
#' of the fit is not desired.  Note, for FRPP this argument should 
#' remain NULL and FRPP
#' must be established in the lm.rrpp fit (RRPP = FALSE).  If the 
#' null model is uncertain, 
#' using \code{\link{reveal.model.designs}} will help elucidate the 
#' inherent null model used.
#' @param groups A factor or vector that is coercible into a factor, 
#' describing the levels of
#' the groups for which to find LS means or slopes.  Normally this 
#' factor would be part of the 
#' model fit, but it is not necessary for that to be the case in order 
#' to obtain results.
#' @param covariate A numeric vector for which to calculate slopes for 
#' comparison  If NULL, 
#' LS means will be calculated instead of slopes.  Normally this variable 
#' would be part of the 
#' model fit, but it is not necessary for that to be the case in order 
#' to obtain results.
#' @param verbose A logical value for whether to include extra results, specifically 
#' the Mahalanobis distances among means, which requires the calculation of residual 
#' covariance matrices for each permutation.  This should be generally
#' FALSE, unless Mahalanobis distances are desired, in which case it must be TRUE.  
#' Verbose computations can require much additional time.
#' @param print.progress If a null model fit is provided, a logical 
#' value to indicate whether analytical 
#' results progress should be printed on screen.  Unless large data 
#' sets are analyzed, this argument 
#' is probably not helpful.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' 
#' @return An object of class \code{pairwise} is a list containing the following
#' \item{LS.means}{LS means for groups, across permutations.}
#' \item{slopes}{Slopes for groups, across permutations.}
#' \item{means.dist}{Pairwise distances between means, across permutations.}
#' \item{std.means.dist}{Pairwise distances between means, across permutations, standardized.}
#' \item{mah.means.dist}{Pairwise Mahalanobis distances between means, across permutations.}
#' \item{means.vec.cor}{Pairwise vector correlations between 
#' means, across permutations.}
#' \item{means.lengths}{LS means vector lengths, by group, across permutations.}
#' \item{means.diff.length}{Pairwise absolute differences between 
#' mean vector lengths, across permutations.}
#' \item{slopes.dist}{Pairwise distances between slopes (end-points), 
#' across permutations.}
#' \item{slopes.vec.cor}{Pairwise vector correlations between slope 
#' vectors, across permutations.}
#' \item{slopes.lengths}{Slope vector lengths, by group, across permutations.}
#' \item{slopes.diff.length}{Pairwise absolute differences between 
#' slope vector lengths, across permutations.}
#' \item{n}{Sample size}
#' \item{p}{Data dimensions; i.e., variable number}
#' \item{PermInfo}{Information for random permutations, passed on 
#' from lm.rrpp fit and possibly
#' modified if an alternative null model was used.}
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A 
#' method for analysis of phenotypic change for phenotypes described
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C and M.L. Collyer. 2018. Multivariate 
#' phylogenetic ANOVA: group-clade aggregation, biological 
#' challenges, and a refined permutation procedure. Evolution. In press.
#' @seealso \code{\link{lm.rrpp}}
#' @examples 
#' \dontrun{
#' # Examples use geometric morphometric data on pupfishes
#' # See the package, geomorph, for details about obtaining such data
#'
#' # Body Shape Analysis (Multivariate) --------------
#' 
#' data("Pupfish")
#' 
#' # Note:
#' 
#' dim(Pupfish$coords) # highly multivariate!
#' 
#' Pupfish$logSize <- log(Pupfish$CS) 
#'
#' # Note: one should use all dimensions of the data but with this 
#' # example, there are many.  Thus, only three principal components 
#' # will be used for demonstration purposes.
#' 
#' Pupfish$Y <- ordinate(Pupfish$coords)$x[, 1:3]
#' 
#' ## Pairwise comparisons of LS means
#' 
#' # Note: one should increase RRPP iterations but a 
#' # smaller number is used here for demonstration 
#' # efficiency.  Generally, iter = 999 will take less
#' # than 1s for these examples with a modern computer.
#' 
#' fit1 <- lm.rrpp(Y ~ logSize + Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 199) 
#' summary(fit1, formula = FALSE)
#' anova(fit1) 
#'
#' pup.group <- interaction(Pupfish$Sex, Pupfish$Pop)
#' pup.group
#' PW1 <- pairwise(fit1, groups = pup.group)
#' PW1
#' 
#' # distances between means
#' summary(PW1, confidence = 0.95, test.type = "dist") 
#' summary(PW1, confidence = 0.95, test.type = "dist", stat.table = FALSE)
#' 
#' # absolute difference between mean vector lengths
#' summary(PW1, confidence = 0.95, test.type = "DL") 
#' 
#' # correlation between mean vectors (angles in degrees)
#' summary(PW1, confidence = 0.95, test.type = "VC", 
#'    angle.type = "deg") 
#'
#' # Can also compare the dispersion around means
#' summary(PW1, confidence = 0.95, test.type = "var")
#' 
#' ## Pairwise comparisons of slopes
#' 
#' fit2 <- lm.rrpp(Y ~ logSize * Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 199) 
#' summary(fit2, formula = FALSE)
#' anova(fit1, fit2)
#'
#' # Using a null fit that excludes all factor-covariate 
#' # interactions, not just the last one  
#' 
#' PW2 <- pairwise(fit2, fit.null = fit1, groups = pup.group, 
#' covariate = Pupfish$logSize, print.progress = FALSE) 
#' PW2
#' 
#' # distances between slope vectors (end-points)
#' summary(PW2, confidence = 0.95, test.type = "dist") 
#' summary(PW2, confidence = 0.95, test.type = "dist", stat.table = FALSE)
#' 
#' # absolute difference between slope vector lengths
#' summary(PW2, confidence = 0.95, test.type = "DL") 
#' 
#' # correlation between slope vectors (and angles)
#' summary(PW2, confidence = 0.95, test.type = "VC",
#'    angle.type = "deg") 
#'    
#' # Can also compare the dispersion around group slopes
#' summary(PW2, confidence = 0.95, test.type = "var")
#' }
#' 
pairwise <- function(fit, fit.null = NULL, groups, covariate = NULL, 
                     verbose = FALSE, print.progress = FALSE) {
  fitf <- fit
  term.labels <- fit$LM$term.labels
  n <- fit$LM$n
  p <- fit$LM$p
  
  res.df <- fit$ANOVA$df[(length(fit$ANOVA$df) - 1)]
  SE <- sqrt(fit$ANOVA$RSS.model[NROW(fit$ANOVA$RSS.model), ] / res.df)
  
  if(verbose){
    if(res.df < fit$LM$p) {
      cat("\nVerbose results were requested but the residual covariance matrix is singular.\n",
          "Therefore, residual covariance matrices will not be estimated for every random permutation.\n")
      verbose = FALSE
    }
  }
  
  PermInfo <- getPermInfo(fit, attribute = "all")
  Models <- getModels(fit, attribute = "all")
  ind <- PermInfo$perm.schedule
  perms <- length(ind)
  gls <- fitf$LM$gls
  if(!inherits(fitf, "lm.rrpp")) 
    stop("The model fit must be a lm.rrpp class object")
  
  if(!is.null(fit.null)) {
    PermInfo$perm.method <- "RRPP"
    if(!inherits(fit.null, "lm.rrpp")) 
      stop("The null model fit must be a lm.rrpp class object")
    
    if(print.progress){
      cat(paste("\nCoefficients estimation from RRPP on null model:", 
                perms, "permutations.\n"))
      pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
    }
  }
  
  fit <- NULL
  
  k <- length(term.labels)
  kk <- length(Models$full)
  if(k != kk) {
    cat("Because the linear model design matrix is not full rank, there might be an issue.\n")
    cat("If there is a subsequent error, this is probably why.\n\n")
  }
  
  Y <- fitf$LM$Y
  if(gls) {
    if(!is.null(fitf$LM$Cov) && is.null(fitf$LM$PCov))
       fitf$LM$Pcov <- Cov.proj(fitf$LM$Cov)
    Y <- if(!is.null(fitf$LM$Pcov)) fitf$LM$Pcov %*% Y else
      Y %*% sqrt(fitf$LM$weights)
  }
  
  X <- Models$reduced[[max(1, kk)]]$X
  if(!is.null(fit.null)) {
    X <- fit.null$LM$X
    if(fit.null$LM$gls) {
      if(!is.null(fit.null$LM$Cov) && is.null(fit.null$LM$PCov))
        fit.null$LM$Pcov <- Cov.proj(fit.null$LM$Cov)
      X <- if(!is.null(fit.null$LM$Pcov)) fit.null$LM$Pcov %*% X else
        X %*% sqrt(fit.null$LM$weights)
    }
  }
  
  if(k > 0) {
    lmf <- LM.fit(X, Y)
    fitted <-lmf$fitted.values
    res <- lmf$residuals
  } else {
    fitted <- matrix(0, n, p)
    res <- Y
  }
  
  rrpp.args <- list(fitted = as.matrix(fitted), residuals = as.matrix(res),
                    ind.i = NULL, o = NULL)
  
  o <- if(!is.null(fitf$LM$offset)) fitf$LM$offset else NULL
  offset <- if(!is.null(o)) TRUE else FALSE
  rrpp.args$o <- o
  
  rrpp <- function(fitted, residuals, ind.i, offset = FALSE, o = NULL) {
    if(offset) fitted + residuals[ind.i,] - o else
      fitted + residuals[ind.i,]
  }
  
  rrpp.args$o <- if(!is.null(fitf$LM$offset)) fitf$LM$offset else NULL
  rrpp.args$offset <- if(!is.null(o)) TRUE else FALSE
  
  if(gls) {
    if(!is.null(fitf$LM$Cov) && is.null(fitf$LM$PCov)){
      Pcov <- Cov.proj(fitf$LM$Cov)
      Qf <- QRforX(Pcov %*% fitf$LM$X, reduce = FALSE)
    } else {
      Qf <- QRforX(fitf$LM$X * sqrt(fit$LM$weights), reduce = FALSE)
    }
  } else Qf <- QRforX(fitf$LM$X, reduce = FALSE)

  H <- tcrossprod(solve(Qf$R), Qf$Q)
  getCoef <- function(y) H %*% y
  
  if(is.null(fitf$LM$random.coef) || !is.null(fit.null)) {
    
    coef.n <- lapply(1:perms, function(j){
      step <- j
      if(print.progress && !is.null(fit.null)) 
        setTxtProgressBar(pb,step)
      rrpp.args$ind.i <- ind[[j]]
      y <- do.call(rrpp, rrpp.args)
      getCoef(y)
    })
    
    fitf$LM$random.coef <- vector("list", length = kk)
    names(fitf$LM$random.coef) <- term.labels
    fitf$LM$random.coef[[max(1, kk)]] <- coef.n
  }
  
  step <- perms + 1
  if(print.progress && !is.null(fit.null)) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  
  groups <- factor(groups)
  n_groups <- tapply(groups, groups, length)
  
  
  if(length(groups) != n) 
    stop("The length of the groups factor does not match the number of observations 
         in the lm.rrpp fit")
  if(!is.null(covariate)) {
    if(!is.numeric(covariate)) 
      stop("The covariate argument must be numeric")
    if(inherits(covariate, "matrix")) 
      stop("The covariate argument can only be a single covariate, not a matrix")
    if(length(unique(covariate)) < 2) 
      stop("The covariate vector appears to have no variation.")
    if(length(covariate) != n)
      stop("The length of the covariate does not match the number of observations 
           in the lm.rrpp fit")
  }
  if(is.null(covariate)) {
    
    means <- getLSmeans(fit = fitf, g = groups)
    slopes <- NULL
    means.dist <- lapply(means, function(x) as.matrix(dist(x)))
    pw_matches <- combn(n_groups, 2, simplify = FALSE)
    pw.sample.adjust <- lapply(pw_matches, function(x) 1/x[1] + 1/x[2])
    ng <- length(n_groups)
    pw_correct <- as.dist(matrix(0, ng, ng))
    pw_correct[1:length(pw.sample.adjust)] <- sqrt(unlist(pw.sample.adjust))
    pw_correct <- as.matrix(pw_correct)
    diag(pw_correct) <- 1
    std.means.dist <- Map(function(m, s) m/(s * pw_correct), means.dist, as.list(SE))
    means.vec.cor <- lapply(means, vec.cor.matrix)
    means.length <- lapply(means, function(x) sqrt(rowSums(x^2)))
    means.diff.length <- lapply(means.length, function(x) 
      as.matrix(dist(as.matrix(x))))
    if(verbose){
      if(print.progress){
        cat(paste("\nEstimation of residual covariance matrices:", 
                  perms, "permutations.\n"))
        pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
      }
      
      ResCovs <- lapply(1:perms, function(j){
        rrpp.args$ind.i <- ind[[j]]
        y <- do.call(rrpp, rrpp.args)
        yh <- fastFit(Qf$Q, y, n, p)
        R <- y - yh
        step <- j
        if(print.progress) 
          setTxtProgressBar(pb,step)
        crossprod(R) / res.df
      })
      
      step <- perms + 1
      if(print.progress) {
        setTxtProgressBar(pb,step)
        close(pb)
      }
        
      if(print.progress){
        cat(paste("\nCalculation of pairwise Mahalanobis distances:", 
                  perms, "permutations.\n"))
        pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
      }
      
      mah.means.dist <- lapply(1:perms, function(j){
        m <- means[[j]]
        s <- ResCovs[[j]]
        step <- j
        md <- mahal_dist(m, s, method = "euclidean")
        if(print.progress) 
          setTxtProgressBar(pb,step)
        as.matrix(md)
      })
      step <- perms + 1
      if(print.progress) {
        setTxtProgressBar(pb,step)
        close(pb)
      }
    } else mah.means.dist <- NULL
    
    slopes.length <- NULL
    slopes.dist <- NULL
    slopes.diff.length <- NULL
    slopes.vec.cor <- NULL
    
  } else {
    
    means <- NULL
    slopes <- getSlopes(fit = fitf, x = covariate, g = groups)
    means.dist <- NULL
    std.means.dist <- NULL
    mah.means.dist <- NULL
    means.vec.cor <- NULL
    means.length <- NULL
    means.diff.length <- NULL
    
    slopes.length <- lapply(slopes, function(x) sqrt(rowSums(x^2)))
    slopes.dist <- lapply(slopes, function(x) as.matrix(dist(as.matrix(x))))
    slopes.diff.length <- lapply(slopes.length, function(x) as.matrix(dist(as.matrix(x))))
    slopes.vec.cor <- lapply(slopes, vec.cor.matrix)
  }
  
  if(gls) res <- fitf$LM$gls.residuals else res <- fitf$LM$residuals
  disp.args <- list(res = as.matrix(res), ind.i = NULL, x = model.matrix(~groups + 0))
  if(gls) disp.args$x <- fitf$LM$Pcov %*% disp.args$x
  g.disp <- function(res, ind.i, x) {
    r <- res[ind.i,]
    if(NCOL(r) > 1) d <- apply(r, 1, function(x) sum(x^2)) else
      d <- r^2
    coef(lm.fit(as.matrix(x), d))
  }
  
  vars <- sapply(1:perms, function(j){
    disp.args$ind.i <- ind[[j]]
    do.call(g.disp, disp.args)
  })
  
  rownames(vars) <- levels(groups)
  colnames(vars) <- c("obs", paste("iter", 1:(perms -1), sep = "."))
  
  out <- list(LS.means = means, slopes = slopes, 
              means.dist = means.dist, std.means.dist = std.means.dist,
              mah.means.dist = mah.means.dist, 
              means.vec.cor = means.vec.cor, 
              means.length = means.length, 
              means.diff.length = means.diff.length,
              slopes.dist = slopes.dist, slopes.vec.cor = slopes.vec.cor,
              slopes.length = slopes.length, 
              slopes.diff.length = slopes.diff.length,
              vars = vars, SE = SE, n = n, p = p, PermInfo = fitf$PermInfo)
  
  class(out) <- "pairwise"
  out
}
