#' Classification probabilities for lm.rrpp model fits
#'
#' Function calculates posterior probabilities for test data, based on factors in a lm.rrpp fit and 
#' values for continuous covariates.
#'
#' This function uses a \code{\link{lm.rrpp}} fit, optional covariate values, optional prior probabilities,
#' and an optional cross-validation procedure to estimate posterior probabilities.  One can choose the 
#' principal component subspace in which to do the classification.  For most cases, 100 percent of the PC space is preferred; 
#' however, when there is a strong signal - groups have different means and observations within groups are close to means - it 
#' might be desirable to ascertain if the signal remains strong in fewer dimensions.  The pc.space
#' argument can be toggeled, and one might wish to run several analyses with different selections, to get a feel
#' for how classification improves with increased data dimensionality.
#' 
#' For high dimensional data, training and testing data are first projected onto the principal components of the 
#' appropriate error covariance matrix from the \code{\link{lm.rrpp}} fit.  (The function will detect whether ordinary 
#' least squares or generalized least squares was used.)  The model covariance is then re-estimated in the 
#' defined amount of the PC space.  A ridge regularization step is also employed if it is determined that the covariance 
#' matrix is ill-conditioned.  In this space, generalized (Mahalanobis) distances are calculated and used to calculate 
#' Bayesian posterior probabilities, based on a multivariate normal density function.  An option for user-supplied 
#' prior probabilities is available, or a user can choose between equal or sample-size corrected prior probabilities.
#' 
#' This function should produce similar results as \code{\link[MASS]{lda}} for models that can be used in that 
#' function, but it will also calculate least squares means if there are covariates in the \code{\link{lm.rrpp}} fit.  
#' Furthermore, the user can choose specific covariate values for estimation.  For example, one could make classifications
#' at small organism size or large size (rather than mean size), if size is a covariate.  This is helpful if there are 
#' factor-covariate interactions to consider, as classification could differ across the span of a covariate.
#' The function does not find canonical (discriminant) axes, like \code{\link[MASS]{lda}}, however these can be calculated from 
#' \code{\link{lm.rrpp}} output.  Whereas as \code{\link[MASS]{lda}} only allows a single factor model, \code{classify} can
#' be used with various model types, and "discriminant" axes would not make as much sense for complex models.  Finally,
#'  \code{\link[MASS]{lda}} scales data and model design matrices to find (transformed) discriminant functions for within-group
#'  covariances that are forced to be spherical; i.e., the discriminant functions are normalized.  This is not a step
#'  taken with \code{classify}, as model design matrices might be more complex than single-factor designs, and obtaining 
#'  discriminant functions is not the goal.  Nevertheless, posterior probabilities will differ slightly between the two
#'  functions as a result.
#' 
#' Regardless of variables input, data are projected onto PCs.  The purpose of this function is to predict 
#' group association, and working in PC space facilitates this objective.
#' 
#' This is a new function and not all limits and scenarios have been tested before its release.  Please report 
#' any issues or limitations or strange results to the maintainer.  (The function will likely be updated frequently
#' until all kinks are worked out.)  
#' 

#' @param fit A linear model fit using \code{\link{lm.rrpp}}.
#' @param test.data An optional matrix (or object coercible to a matrix) for classification.  If NULL,
#' all observed data are used.
#' @param covariate.values An optional data frame (with one value per covariate) to control covariate values.  
#' For example, covariate.values = data.frame(x1 = 2, x2 = 4).  This option is used if estimation at a point on a
#' covariate axis other than the mean value is desired.
#' @param prior An option vector with values between 0 and 1 in the order of groups.  If the groups extracted from 
#' the \code{\link{lm.rrpp}} fit or there order are not clear, it is wise to use the inherent.groups = TRUE option first.
#' This will display the group information.  If NULL, prior probabilities will be based on group size.  A choice of 
#' "equal" will force equal prior probabilities.
#' @param CV A logical argument to indicate whether a leave-one-out cross-validation (Jackknife) procedure should be
#' performed.  This is only a logical choice if test data are the same data used for the \code{\link{lm.rrpp}} fit.
#' In other cases, the argument will be ignored.  The default is CV = TRUE.
#' @param pc.space A value between 0 and 1 to indicate the portion of PC space to use for classification.  The default is 
#' 1.0 (the entire PC space).  Note that this is an upper limit, and the PCs included might (and probably do) comprise
#' less space.  
#' @param inherent.groups A logical argument in case one wishes to have the inherent groups in the model fit revealed.  If 
#' TRUE, no other analysis will be done than to reveal the groups.  This argument should always be FALSE to perform 
#' a classification analysis.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{classify} is a list containing the following
#' \item{means}{LS means (PC scores) for inherent groups.}
#' \item{group}{A factor derived for groups.}
#' \item{group.n}{Sample size for each group.}
#' \item{test.data}{PC scores for the test data used.}
#' \item{Mah.dist.sq}{Generalied Mahalanobis squared distance from test data to group means.}
#' \item{prior}{Prior probabilities.}
#' \item{posterior}{Posterior probabilities.}
#' \item{class}{Classification based on highest posterior probability}
#' \item{class.matrix}{An indicator matrix for classifications.}
#' @examples 
#' 
#' # Using the Pupfish data (see lm.rrpp help for more detail)
#' 
#' data(Pupfish)
#' Pupfish$logSize <- log(Pupfish$CS) 
#' fit <- lm.rrpp(coords ~ logSize + Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 0)
#' 
#' classify(fit, inherent.groups = TRUE) # see groups available
#' class1 <- classify(fit, CV = TRUE)
#' summary(class1)
#' 
#' # Perfect classification
#' 
#' table(interaction(Pupfish$Sex, Pupfish$Pop), class1$class)
#' 
#' class2 <- classify(fit, CV = TRUE, pc.space = 0.6)
#' summary(class2)
#' table(interaction(Pupfish$Sex, Pupfish$Pop), class2$class)
#' 
#' # Allow for heterogenous slopes in the model
#' fit2 <- lm.rrpp(coords ~ logSize * Sex * Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 0)
#' 
#' # Classification for small-sized fish (previous example for mean-sized fish)
#' 
#' class3 <- classify(fit2, CV = TRUE, 
#' covariate.values = data.frame(logSize = 4))
#' summary(class3)
#' 
#' # Classification for large-sized fish 
#' 
#' class4 <- classify(fit2, CV = TRUE, 
#' covariate.values = data.frame(logSize = 4.5))
#' summary(class4)

classify <- function(fit, test.data = NULL, covariate.values = NULL, 
                     prior = NULL, CV = TRUE, pc.space = 1.0, 
                     inherent.groups = FALSE){
  
  dat <- fit$LM$data
  gls <- fit$LM$gls
  res <- fit$LM$wResiduals
  n <- fit$LM$n
  p <- fit$LM$p
  df <- fit$ANOVA$df
  dfe <- df[length(df) - 1]
  X <- fit$LM$X 
  Y <- fit$LM$Y 
  wt <- fit$LM$weights
  
  if(is.null(test.data) && CV) CV = TRUE else CV = FALSE
  if(is.null(test.data)) td <- test.data <- fit$LM$Y else td <- test.data
  if(is.vector(td)) td <- matrix(td, 1, length(td))
  if(is.data.frame(td))  td <- as.matrix(td)
  dims <- dim(td)
  td.n <- dims[1]
  if(dims[2] != p) stop("\nThe number of variables in the test data does not match the number of variables in the model.\n")
  
  Ymean <- colMeans(Y)
  
  if(pc.space > 1) pc.space = 1
  if(gls) pca <- prcomp(fit$LM$gls.residuals, tol= 0.001) else 
    pca <- prcomp(fit$LM$wResiduals, tol= 0.001)
  d <- pca$sdev^2
  cd <- cumsum(d)/sum(d)
  k <- which(cd < pc.space)
  if(length(k) < 1) k = 1
  if(length(cd) > 1) P <- as.matrix(center(Y) %*% pca$rotation[, k]) 
  PTD <- as.matrix((td - matrix(rep(Ymean, each = td.n), td.n, p)) %*% 
                     pca$rotation[, k])
  
  rdf <- rrpp.data.frame(P = P, X = X)
  
  if(gls) {
    pfit <- lm.rrpp(P ~ X + 0, print.progress = FALSE, 
                    Cov = fit$LM$Cov, data = rdf, weights = wt, iter = 0)
    Cov <- (crossprod(pfit$LM$gls.residuals, fast.solve(pfit$LM$Cov)) %*%
              pfit$LM$gls.residuals) / dfe
    if(kappa(Cov) > 1e10) Cov <- RiReg(Cov, pfit$LM$gls.residuals)
    
  } else {
    
    pfit <- lm.rrpp(P ~ X + 0, print.progress = FALSE, 
                    data = rdf, weights = wt, iter = 0)
    Cov <- crossprod(pfit$LM$wResiduals) / dfe
    if(kappa(Cov) > 1e10) Cov <- RiReg(Cov, pfit$LM$wResiduals)
  }
  
  
  
  if(inherent.groups){
    dat.class <- sapply(dat, class)
    fac.list <- which(dat.class == "factor")
    if(length(fac.list) == 0) group <- factor(rep(1, n)) else {
      datf <- dat[fac.list]
      group <- factor(apply(datf, 1, paste, collapse = "."))
    }
    if(nlevels(group) == 1) 
      cat("\n There appears to be no factors in the lm.rrpp fit\n\n") else{
        cat("\n These are the apparent groups (and sizes) in the lm.rrpp fit\n\n")
        print(as.table(by(group, group, length)))
        
      }
  } else {
    
    invCov <- fast.solve (Cov)
    
    dat.class <- sapply(dat, class)
    fac.list <- which(dat.class == "factor")
    if(length(fac.list) == 0) 
      stop("\nThere are no factors in the lm.rrpp fit from which to create groups.") else {
        datf <- dat[fac.list]
        group <- factor(apply(datf, 1, paste, collapse = "."))
      }
    
    cov.list <- which(dat.class == "numeric")
    if(length(cov.list) > 0) for(i in 1:length(cov.list)) {
      k <- cov.list[i]
      dat[[k]] <- mean(dat[[k]])
    }
    
    if(!is.null(covariate.values)) {
      cv <- covariate.values
      if(is.matrix(cv)) cv <- as.data.frame(cv)
      if(!is.data.frame(cv)) stop("\ncovariate.values must be a data.frame object with one row of values\n")
      cv.match <- names(cv) %in% names(dat)
      if(any(!cv.match)) stop("\nSome covariates provided are not found in the model.\n")
      cv.match <- match(names(cv), names(dat))
      for(i in 1:length(cv)) {
        k <- cv.match[i]
        if(dat.class[[k]] == "factor") dat[[k]] <- dat[[k]] else dat[[k]] <- cv[[i]]
      }
    }
    
    Xls <- model.matrix(fit$LM$Terms, data = dat) * sqrt(fit$LM$weights)
    if(NROW(unique(Xls)) != nlevels(group)) cat("\nWarning: Number of groups does not match design matrix levels")
    if(nlevels(group) == 1) xg = matrix(1,n) else xg <-model.matrix(~group + 0)
    
    if(gls) targ <- means <- as.matrix(lm.fit(xg, Xls %*% pfit$LM$gls.coefficients)$coefficients) else 
      targ <- means <- as.matrix(lm.fit(xg, Xls %*% pfit$LM$coefficients)$coefficients)
    
    if(CV){
      gp.check <- by(group, group, length)
      if(any(gp.check < 2))
        stop("\nNot every group has more than 1 observation; cross-validation not possible.\n")
      if(any(gp.check < 5))
        cat("\nWarning: Some groups have small sample sizes, which could lead to poor results\n\n")
      
      betas <- lapply(1:n, function(j){
        leaveOneOut(X * sqrt(fit$LM$weights), P, j)
      })
      xg <- model.matrix(~ group + 0)
      targ <- lapply(1:n, function(j){
        Xls %*% betas[[j]]
        lm.fit(xg, Xls %*% betas[[j]])$coefficients
      })
      nt <- NROW(targ[[1]])
      MD <- t(sapply(1:n, function(j){
        y <- PTD[j,]
        y <- matrix(rep(y, each = nt), nt, NCOL(P))
        targ.j <- targ[[j]]
        td.j <- y - targ.j
        diag(tcrossprod(td.j, invCov) %*% t(td.j))
      }))
    } else 
      
      MD <- sapply(1:NROW(targ), function(j){
        nt <- NROW(td)
        tgt <- matrix(rep(targ[j, ], each = nt), nt, NCOL(targ))
        td.j <- PTD - tgt
        diag(tcrossprod(td.j, invCov) %*% t(td.j))
      })
    
    if(is.vector(MD)) MD <- matrix(MD, 1, length(MD))
    
    gp.n <- nlevels(group)
    
    if(!is.null(prior) && prior == "equal") prior <- rep(1/gp.n, gp.n)
    if(is.null(prior)) prior <- as.vector(by(group, group, length)/length(group))
    if(!is.null(prior) && is.vector(prior)) {
      if(gp.n != length(prior)) 
        stop(cat("\nNumber of inherent groups =", gp.n, "but number of priors =", length(prior)), "\n")
    }
    if(!is.vector(prior)) stop("\nprior must be a vector with same length as inherent groups.\n")
    if(any(prior > 1)) stop("\nprior values must be between 0 and 1.\n")
    if(any(prior < 0)) stop("\nprior values must be between 0 and 1.\n")
    if(sum(prior) != 1) cat("Warning: prior values do not sum to 1 (and, therefore, seem illogical)\n\n")
    prior <- diag(prior, gp.n)
    
    detCov <- det(Cov)
    
    md.est <- function(d){
      sqrt(2 * pi * detCov) * exp(-0.5*d)
    }
    mdEst <- apply(MD, c(1,2), md.est)
    num <- mdEst %*% prior
    den <- apply(num, 1, sum)
    posts <- num/den

    class <- apply(posts, 1, which.max)
    class.matrix <- t(round(apply(posts, 1, function(x) x/max(x))))
    
    rownames(means) <- levels(group)
    colnames(means) <- colnames(PTD)
    colnames(posts) <- colnames(MD) <- levels(group)
    rownames(posts) <- rownames(MD) <- rownames(PTD)
    group.n <- as.vector(by(group, group, length))
    names(group.n) <- levels(group)
    prior <- diag(prior)
    names(prior) <- colnames(MD)
    out <- list(means = means, group = group, group.n = group.n,
                test.data = PTD, Mah.dist.sq = MD, prior = prior,
                posterior = posts, class = class, class.matrix = class.matrix)
    class(out) <- "classify"
    out
    
  }
}

