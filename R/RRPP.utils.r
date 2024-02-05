## rrpp.data.frame

#' Handle missing values in rrpp.data.frame objects
#'
#' @param object object (from \code{\link{rrpp.data.frame}})
#' @param ... further arguments (currently not used)
#' @method na.omit rrpp.data.frame
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples
#' y <- matrix(rnorm(15), 5, 3)
#' x <- rnorm(5)
#' rdf <- rrpp.data.frame(x = x, y = y, d = dist(y))
#' rdf$x[1] <- NA # create missing data
#' rdf
#' 
#' ndf <- na.omit(rdf)
#' ndf

na.omit.rrpp.data.frame <- function(object, ...) {
  nms <- names(object)
  classes <- unlist(lapply(object, function(x) class(x)[1]))
  subDF <- object[!is.na(match(classes, c("numeric", "factor", "integer", 
                                          "character", "matrix", "logical")))]
  sub.classes <- classes[!is.na(match(classes, c("numeric", "factor", "integer", 
                                                 "character", "matrix", "logical")))]
  oDF <- object[is.na(match(classes, c("numeric", "factor", "integer", 
                                       "character", "matrix", "logical")))]
  o.classes <- classes[is.na(match(classes, c("numeric", "factor", "integer", 
                                              "character", "matrix", "logical")))]
  subDF <- as.data.frame(subDF)
  newDF <- na.omit(subDF)
  omits <- attr(newDF, "na.action")
  
  for(i in 1:length(o.classes)){
    
    if(o.classes[[i]] == "array") {
      dims <- dim(oDF[[i]])
      if(length(dims) != 3)
        stop("Data are neither a vector, matrix, nor appopriate array.\n",
             call. = FALSE)
      oDF[[i]] <- oDF[[i]][,,-omits]
    }
    
    if(o.classes[[i]] == "phylo") {
      
      cat("Part of this data frame is a class phylo object\n")
      cat("It is currently not possible to prune the tree according to missing data\n")
      cat("The following actions are recommended:\n")
      cat("1. Make a data frame with all objects or variables except the phylo object\n")
      cat("2. Omit missing data to create a new data frame\n")
      cat("3. Add a pruned tree to the new data frame; e.g., newDF$tree <- myPrunedTree\n")
      stop("na.omit terminated", call. = FALSE)

        stop("Data are netieher a vector, matrix, nor appopriate array.\n",
             call. = FALSE)
      oDF[[i]] <- oDF[[i]][,,-omits]
    }
    
    if(o.classes[[i]] == "dist") {
      
      d <- as.matrix(oDF[[i]])
      d <- d[-omits, -omits]
      oDF[[i]] <- as.dist(d)
    }
    
  }
  
  outDF <- c(as.list(newDF), oDF)
  class(outDF) <- "rrpp.data.frame"
  attr(outDF, "na.action") <- omits
  outDF[nms]
  
}

## lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{lm.rrpp}})
#' @param ... other arguments passed to print/summary
#' @method print lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.lm.rrpp <- function(x, ...){
  cat("\nLinear Model fit with lm.rrpp\n")
  LM <- x$LM
  PI <- getPermInfo(x, attribute = "all")
  AN <- x$ANOVA
  if(is.null(AN$Fs)) {
    AN2 <- getANOVAStats(x, stat = "all")
    AN[c("MS", "Rsq", "F", "cohenf")] <- AN2[c("MS", "Rsq", "F", "cohenf")]
  }
    
  if(!is.null(x$LM$dist.coefficients)) 
    dv <- "(dimensions of data after PCoA of distance matrix)" else
    dv <- " "
  cat(paste("\nNumber of observations:", LM$n))
  cat(paste("\nNumber of dependent variables:", LM$p, dv))
  cat(paste("\nData space dimensions:", LM$p.prime, dv))
  if(!is.null(AN$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", 
                                     AN$SS.type))
  if(!is.null(PI)) cat(paste("\nNumber of permutations:", PI$perms))
  cat("\nCall: ")
  cat(deparse(x$call), fill=TRUE)
  invisible(x)
}

## lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param object print/summary object (from \code{\link{lm.rrpp}})
#' @param formula Logical argument for whether to include formula in summary table
#' @param ... other arguments passed to print/summary
#' @method summary lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.lm.rrpp <- function(object, formula = TRUE, ...){
  if(inherits(object, "manova.lm.rrpp")) out <- summary.manova.lm.rrpp(object) else {
    x <- object
    LM <- x$LM
    PI <- getPermInfo(x, attribute = "all")
    AN <- x$ANOVA
    if(is.null(AN$Fs)) {
      AN2 <- getANOVAStats(object, stat = "all")
      AN[c("SS", "MS", "Rsq", "Fs", "cohenf")] <-
        AN2[c("SS", "MS", "Rsq", "Fs", "cohenf")]
    }
      
    Models <- getModels(x, attribute = "all")
    perms <- PI$perms
    dv <- LM$dist.coefficients
    n <- LM$n
    p <- LM$p
    p.prime <- LM$p.prime
    SS <- AN$SS
    SSM.obs <- SS[,1]
    RSS <- AN$RSS
    TSS <- AN$TSS
    RSS.model <- AN$RSS.model[nrow(AN$RSS.model),]
    if(is.null(RSS)) RSS <- RSS.model
    if(is.null(RSS)) TSS <- RSS.model
    SS.type <- AN$SS.type
    trms <- LM$term.labels
    k <- length(trms)
    kk <- length(Models$full)
    if(k > kk){
      k <- kk
      trms <- names(Models$full)
    }
    
    if(k > 0) {
      df <- AN$df
      dfe <- df[k+1]
      dfM <- sum(df[1:k])
      SSM <- (TSS[1,] - RSS.model)
      SSM.obs <- SSM[1]
      Rsq <- SSM.obs/TSS[1]
      Fs <- (SSM/dfM)/(RSS.model/dfe)
      Fs.obs <- Fs[1]
      P <- pval(as.vector(Fs))
      Z <- effect.size(as.vector(Fs))
      
    } else {
      df <- 0
      dfe <- AN$df
      dfM <- 1
      SSM <- Fs <- Z <- Rsq <- SSM.obs <- Fs.obs <- P <- ""
    }
    
    tab <- data.frame(dfM = dfM, dfe = dfe,
                      SSM = SSM.obs, RSS = RSS[1], Rsq = Rsq,
                      F = Fs.obs, Z = Z, P = P)
    dimnames(tab)[[2]] <- c("Df",
                            "Residual Df",
                            "SS",
                            "Residual SS",
                            "Rsq",
                            "F",
                            "Z (from F)",
                            "Pr(>F)")
    if(formula) dimnames(tab)[[1]] <- deparse(formula(x$call$f1)[[3]]) else
      dimnames(tab)[[1]] <- deparse(substitute(object))
    
    if(LM$ols){
      pca.fitted <- prcomp(LM$fitted)
      pca.residuals <- prcomp(LM$residuals)
      pca.total <- prcomp(LM$Y)
    }
    
    
    if(LM$gls) {
      int <- attr(LM$Terms, "intercept")
      
      if(!is.null(LM$Cov) && is.null(LM$Pcov))
        LM$Pcov <- Cov.proj(LM$Cov)
      
      if(is.null(LM$Pcov)) {
        w <- sqrt(LM$weights)
        PY <- LM$Y * w
        PX <- LM$X * w
        Pint <- as.matrix(rep(int, n) * w)
        RM <- lm.fit(PX, PY)$residuals
        Sr <- crossprod(RM)
      }
      if(!is.null(LM$Pcov)) {
        Pcov <- LM$Pcov
        PY <- as.matrix(Pcov %*% LM$Y)
        PX <- as.matrix(Pcov %*% LM$X)
        Pint <- as.matrix(Pcov %*% rep(int, n))
        RM <- lm.fit(PX, PY)$residuals
        Sr <- crossprod(RM)
      }
      
      RT <- lm.fit(Pint, PY)$residuals
      Sy <- crossprod(RT)
      Sf <- Sy - Sr
      
      Cf <- Sf/(n - 1)
      Cr <- Sr/(n - 1)
      Cy <- Sy/(n - 1)
      
      pca.fitted <- pca.residuals <- pca.total <- list()
      
      svd.Y <- svd(Cy)
      keep <- which(zapsmall(svd.Y$d) > 0)
      pca.total$sdev <- sqrt(svd.Y$d[keep])
      pca.total$rotation <- as.matrix(svd.Y$v[, keep])
      pca.total$x <- RT %*% as.matrix(svd.Y$v[, keep])
      
      svd.f <- svd(Cf)
      keep <- which(zapsmall(svd.f$d) > 0)
      if(length(keep) == 0) {
        pca.fitted$sdev <- pca.fitted$rotation <- pca.fitted$x <- 0
      } else {
        pca.fitted$sdev <- sqrt(svd.f$d[keep])
        pca.fitted$rotation <- as.matrix(svd.f$v[, keep])
        pca.fitted$x <- LM$gls.fitted %*% as.matrix(svd.f$v[, keep])
      }
      
      svd.r <- svd(Cr)
      keep <- which(zapsmall(svd.r$d) > 0)
      pca.residuals$sdev <- sqrt(svd.r$d[keep])
      pca.residuals$rotation <- as.matrix(svd.r$v[, keep])
      pca.residuals$x <- RM %*% as.matrix(svd.r$v[, keep])
      
    }
    
    
    d.f <- pca.fitted$sdev^2
    d.r <- pca.residuals$sdev^2
    d.t <- pca.total$sdev^2
    
    rank.f <- length(which(zapsmall(d.f) > 0))
    rank.r <- length(which(zapsmall(d.r) > 0))
    rank.t <- length(which(zapsmall(d.t) > 0))
    
    Trace <- zapsmall(c(sum(d.f), sum(d.r), sum(d.t)))
    Proportion <- zapsmall(Trace/sum(d.t))
    Rank <- c(rank.f, rank.r, rank.t)
    
    redund <- data.frame(Trace, Proportion, Rank)
    rownames(redund) <- c("Fitted", "Residuals", "Total")
    eigs.f <- eigs.r <- eigs.t <- rep(NA, rank.t)
    if(rank.f > rank.t) rank.f <- rank.t
    if(rank.r > rank.t) rank.r <- rank.t
    eigs.f[1:rank.f] <- d.f[1:rank.f]
    eigs.r[1:rank.r] <- d.r[1:rank.r]
    eigs.t[1:rank.t] <- d.t[1:rank.t]
    eigs <- as.table(zapsmall(rbind(eigs.f, eigs.r, eigs.t)))
    rownames(eigs) <- c("Fitted", "Residuals", "Total")
    colnames(eigs) <- paste("PC", 1:rank.t, sep="")
    
    reduced <- Models$reduced
    full <- Models$full
    
    if(k > 0) {
      
      if(LM$gls) {
        TY <- if(!is.null(LM$Pcov)) LM$Pcov %*% LM$Y else
          sqrt(LM$weights)
      } else TY <- LM$Y
      
      Ur <- lapply(reduced, function(x) qr.Q(x$qr))
      Uf <- lapply(full, function(x) qr.Q(x$qr))
      
      RR <- Map(function(u) TY - fastFit(u, TY, n, p), Ur)
      RF <- Map(function(u) TY - fastFit(u, TY, n, p), Uf)
      
      SSCP <- lapply(1:length(RF), function(j) crossprod(RR[[j]] - RF[[j]]))
      names(SSCP) <- trms
      SSCP <- c(SSCP, list(Residuals = as.matrix(crossprod(RF[[k]]))))
      
    } else {
      RR <- reduced[[1]]$residuals
      RF <- full[[1]]$residuals
      
      if(LM$gls) {
        if(!is.null(LM$Cov)) {
          if(!is.null(LM$Cov) && is.null(LM$Pcov))
            LM$Pcov <- Cov.proj(LM$Cov)
          RR <- LM$Pcov %*% RR
          RF <- LM$Pcov %*% RF
        } else {
          RR <- RR * sqrt(LM$weights)
          RF <- RF * sqrt(LM$weights)
        }
      }
      
      SSCP <- list(Residuals = as.matrix(crossprod(RF)))
      
    }
    
    out <- list(table = tab, SSCP = SSCP, n = n, p = p, p.prime = p.prime, k = k, 
                perms = perms, dv = dv, SS = SS, SS.type = SS.type, 
                redundancy = redund,
                eigenvalues = eigs, gls = x$LM$gls)
    class(out) <- "summary.lm.rrpp"
  }
  out
}

#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{summary.lm.rrpp}})
#' @param ... other arguments passed to print/summary
#' @method print summary.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.summary.lm.rrpp <- function(x, ...) {
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of observations:", x$n))
  if(!is.null(x$dv)) dv <- "(dimensions of data after PCoA of distance matrix)" else
    dv <- " "
  cat(paste("\nNumber of dependent variables:", x$p, dv))
  cat(paste("\nData space dimensions:", x$p.prime, dv))
  if(!is.null(x$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", 
                                    x$SS.type))
  cat(paste("\nNumber of permutations:", x$perms))
  cat("\n\nFull Model Analysis of Variance\n\n")
  print(x$table)
  cat("\n\nRedundancy Analysis (PCA on fitted values and residuals)\n\n")
  if(x$gls) {
    cat("\nGLS mean used rather than center of gravity.  
        Projection is not orthogonal.\n\n")
  }
  print(x$redundancy)
  cat("\nEigenvalues\n\n")
  print(x$eigenvalues)
  cat("\n")
  invisible(x)
}


## terms.lm.rrpp

#' Extract the terms from an lm.rrpp object
#' 
#' \code{terms.lm.rrpp} returns the terms constructed for an \code{lm.rrpp} object.
#'
#' @param x Object from \code{\link{lm.rrpp}}
#' @param ...	further arguments passed to or from other methods
#' @export
#' @author Michael Collyer
#' @keywords utilities

terms.lm.rrpp <- function(x, ...) return(x$LM$Terms)

## model.frame.lm.rrpp

#' Extract model frame from a lm.rrpp object
#' 
#' \code{model.frame.lm.rrpp} returns the model frame constructed for 
#' an \code{lm.rrpp} object.
#'
#' @param formula Object from \code{\link{lm.rrpp}}
#' @param ...	further arguments passed to or from other methods
#' @export
#' @author Michael Collyer
#' @keywords utilities

model.frame.lm.rrpp <- function(formula, ...) {
  res <- if(inherits(formula, "lm.rrpp")) formula$LM$data else
    model.frame.default(formula)
  return(res)
}

## model.matrix.lm.rrpp

#' Extract the model design matrix from an lm.rrpp object
#' 
#' \code{model.matrix.lm.rrpp} returns the design matrix constructed for 
#' an \code{lm.rrpp} object.
#'
#' @param object Object from \code{\link{lm.rrpp}}
#' @param ...	further arguments passed to or from other methods
#' @export
#' @author Michael Collyer
#' @keywords utilities

model.matrix.lm.rrpp <- function(object, ...) return(object$LM$X)


## coef.lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{coef.lm.rrpp}}
#' @param ... Other arguments passed onto coef.lm.rrpp
#' @method print coef.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.coef.lm.rrpp <- function(x, ...){
  
  if(!x$test) {
    cat("\n\nObserved coefficients\n\n")
    print(x$coef.obs)
    if(!is.null(x$coef.se)) {
      cat("\n\nCoefficient standard errors\n\n")
      print(x$coef.se)
    }
  }
  
  if(x$test){
    
    cat("\nLinear Model fit with lm.rrpp\n")
    cat(paste("\nNumber of observations:", x$n))
    cat(paste("\nNumber of dependent variables:", x$p))
    cat(paste("\nData space dimensions:", x$p.prime))
    cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
    cat(paste("\nNumber of permutations:", x$nperms))
    
    if(is.null(x$stat.table)) {
      cat("\n\nTests on coefficients are not possible (only an intercept).")
      cat("\n\nObserved coefficients\n\n")
      print(x$coef.obs)
    } else {
      rrpp.type <- x$RRPP
      cat("\n\nStatistics (distances) of coefficients with ")
      cat(x$confidence*100, "percent confidence intervals,") 
      cat("\neffect sizes, and probabilities of exceeding observed values 
          based on\n") 
      cat(x$nperms, "random permutations using", rrpp.type, "\n\n")
      print(x$stat.tab)
      cat("\n\n")
    }
  }
  invisible(x)
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{coef.lm.rrpp}}
#' @param ... Other arguments passed onto coef.lm.rrpp
#' @method summary coef.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.coef.lm.rrpp <- function(object, ...){
  print.coef.lm.rrpp(object, ...)
}

## predict.lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{predict.lm.rrpp}}
#' @param PC Logical argument for whether to use predicted values 
#' rotated to their PCs
#' @param ... Other arguments passed onto predict.lm.rrpp
#' @method print predict.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.predict.lm.rrpp <- function(x, PC = FALSE, ...){
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of predictions:", NROW(x$mean)))
  cat(paste("\nConfidence level:", x$confidence*100, "%"))
      cat(paste("\nNumber of bootstrap permutations:", 
                length(x$random.predicted)))
      if(PC)  cat(paste("\nPredicted values are rotated to their PCs"))
      cat("\n\nPredicted values (mean of bootstrapped values):\n\n")
      if(PC) print(x$pc.mean) else print(x$mean)
      cat("\n\n", x$confidence*100, "% Lower confidence limits:\n\n") 
      if(PC) print(x$pc.lcl) else print(x$lcl)
      cat("\n\n", x$confidence*100, "% Upper confidence limits:\n\n") 
      if(PC) print(x$pc.ucl) else  print(x$ucl) 
      cat("\n\n")
      invisible(x)
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{predict.lm.rrpp}}
#' @param ... Other arguments passed onto predict.lm.rrpp
#' @method summary predict.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.predict.lm.rrpp <- function(object, ...){
  print.predict.lm.rrpp(object, ...)
}

## anova.lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{lm.rrpp}})
#' @param ... other arguments passed to print/summary
#' @method print anova.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.anova.lm.rrpp <- function(x, ...) {
  tab <- x$table
  pm <- x$perm.method
  perms <- x$perm.number
  est <- x$est.method
  SS.type <- x$SS.type
  effect.type <- x$effect.type
  
  if(pm == "RRPP") pm <- "Randomization of null model residuals" else
    if(pm == "FMRP") pm <- "Randomization of the full model residuals" else 
      pm <- ("Randomization of raw values (residuals of mean)")
  
  if(est == "GLS") est <- "Generalized Least-Squares (via OLS projection)" else
    est <- "Ordinary Least Squares"
  
  if(!is.null(SS.type)) {
    cat("\nAnalysis of Variance, using Residual Randomization\n")
    cat(paste("Permutation procedure:", pm, "\n"))
    cat(paste("Number of permutations:", perms, "\n"))
    cat(paste("Estimation method:", est, "\n"))
    cat(paste("Sums of Squares and Cross-products: Type", SS.type, "\n"))
    cat(paste("Effect sizes (Z) based on", effect.type, "distributions\n\n"))
    print(tab)
    cat("\nCall: ")
    cat(deparse(x$call), fill=TRUE)
  } else {
    cat("\nAnalysis of Variance, using Residual Randomization\n")
    cat(paste("Permutation procedure:", pm, "\n"))
    cat(paste("Number of permutations:", perms, "\n"))
    cat(paste("Estimation method:", est, "\n"))
    if(NROW(tab) == 1) cat("No model effects; simple summary provided\n\n")
    if(NROW(tab) > 1) cat(paste("Effect sizes (Z) based on", 
                                effect.type, "distributions\n\n"))
    print(tab)
    if(NCOL(tab) == 7) {
      cat("\nCall: ")
      cat(deparse(x$call), fill=TRUE)
    }
  }
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{predict.lm.rrpp}}
#' @param ... Other arguments passed onto predict.lm.rrpp
#' @method summary anova.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.anova.lm.rrpp <- function(object, ...){
  print.anova.lm.rrpp(object, ...)
}



## plot functions

#' Plot Function for RRPP
#' 
#' @param x plot object (from \code{\link{lm.rrpp}})
#' @param type Indicates which type of plot, choosing among diagnostics,
#' regression, or principal component plots.  Diagnostic plots are similar to 
#' \code{\link{lm}} diagnostic plots, but for multivariate data.  Regression plots
#' plot multivariate dispersion in some fashion against predictor values. PC plots
#' project data onto the eigenvectors of the covariance matrix for fitted values.
#' @param resid.type If type = "diagnostics", an optional argument for whether Pearson ("p") or normalized ("n")
#' residuals should be used.  These residuals are the same for ordinary least-squares (OLS) estimation
#' but differ for generalized least-squares (GLS) estimation.  For the latter, normalizing residuals
#' requires multiplying them by the transformation matrix obtained for GLS estimation.
#' @param fitted.type As with resid.type, whether fitted values use observed ("o") or transformed ("t")
#' values.  
#' @param predictor An optional vector if "regression" plot type is chosen, 
#' and is a variable likely used in \code{\link{lm.rrpp}}.
#' This vector is a vector of covariate values equal to the number of observations.
#' @param reg.type If "regression" is chosen for plot type, this argument
#' indicates whether prediction line 
#' (PredLine) or regression score (RegScore) plotting is performed.  
#' For explanation of prediction line,
#' see Adams and Nistri (2010).  For explanation of regression score, see 
#' Drake and Klingenberg (2008).
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' @method plot lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Drake, A. G., and C. P. Klingenberg. 2008. 
#' The pace of morphological change: Historical 
#' transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence 
#' and evolution of foot morphology 
#' in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @examples 
#' 
#' # Univariate example
#' data(PlethMorph)
#' fitGLS <- lm.rrpp(TailLength ~ SVL, data = PlethMorph, Cov = PlethMorph$PhyCov, 
#'  print.progress = FALSE, iter = 0)
#'  
#'  par(mfrow = c(2, 2))
#'  plot(fitGLS)
#'  plot(fitGLS, resid.type = "n") # use normalized (transformed) residuals
#'  plot(fitGLS, resid.type = "n", fitted.type = "t") # use also transformed fitted values
#'  
#'  # Multivariate example
#'  
#' Y <- as.matrix(cbind(PlethMorph$TailLength,
#' PlethMorph$HeadLength,
#' PlethMorph$Snout.eye,
#' PlethMorph$BodyWidth,
#' PlethMorph$Forelimb,
#' PlethMorph$Hindlimb))
#' PlethMorph$Y <- Y
#' fitGLSm <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
#' Cov = PlethMorph$PhyCov,
#' print.progress = FALSE, iter = 0)
#' 
#' par(mfrow = c(2, 2))
#'  plot(fitGLSm)
#'  plot(fitGLSm, resid.type = "n") # use normalized (transformed) residuals
#'  plot(fitGLSm, resid.type = "n", fitted.type = "t") # use also transformed fitted values
#'  par(mfrow = c(1, 1))
#'  
plot.lm.rrpp <- function(x, type = c("diagnostics", "regression",
                                      "PC"), 
                         resid.type = c("p", "n"),
                         fitted.type = c("o", "t"),
                         predictor = NULL,
                          reg.type = c("PredLine", "RegScore"), ...){
  plot_args <- list(...)
  type <- match.arg(type)
  resid.type <- match.arg(resid.type)
  fitted.type <- match.arg(fitted.type)
  
  if(!is.null(x$LM$Cov) && is.null(x$LM$Pcov))
    x$LM$Pcov <- Cov.proj(x$LM$Cov)
  
  if(x$LM$gls) {
    
    r <- as.matrix(x$LM$gls.residuals)
    f <- as.matrix(x$LM$gls.fitted)
    
    if(resid.type == "n"){
      if(!is.null(x$LM$Pcov)) r <- x$LM$Pcov %*% r
      if(!is.null(x$LM$weights)) r <- r * sqrt(x$LM$weights)
    }
    if(fitted.type == "t"){
      if(!is.null(x$LM$Pcov)) f <- x$LM$Pcov %*% f
      if(!is.null(x$LM$weights)) f <- f * sqrt(x$LM$weights)
    }
    
  } else {
    r <- as.matrix(x$LM$residuals)
    f <- as.matrix(x$LM$fitted)
  }
  
  if(is.na(match(type, c("diagnostics", "regression", "PC")))) 
    type <- "diagnostics"
  PL <- Reg.proj <- PC.points <- NULL
  
  if(type == "diagnostics") {
    
    if(x$LM$p == 1) {
      plot_args <- NULL
      
      y <- x$LM$Y
      if(resid.type == "n"){
        if(!is.null(x$LM$Pcov)) y <- x$LM$Pcov %*% y
        if(!is.null(x$LM$weights)) y <- y * sqrt(x$LM$weights)
      }
      
      X <- x$LM$X
      if(fitted.type == "t"){
        if(!is.null(x$LM$Pcov)) X <- x$LM$Pcov %*% X
        if(!is.null(x$LM$weights)) X <- X * sqrt(x$LM$weights)
      }

      X <- as.matrix(X)
      y <- as.matrix(y)
      lm.new <- lm(y ~ X + 0)
      lm.new$call <- x$call
      plot(lm.new, ...)
      
    } else {
      
      plot_args <- NULL
      pca.r <- prcomp(r)
      var.r <- round(pca.r$sdev^2/sum(pca.r$sdev^2)*100,2)
      plot(pca.r$x, pch=19, asp =1,
           xlab = paste("PC 1", var.r[1],"%"),
           ylab = paste("PC 2", var.r[2],"%"),
           main = "PCA Residuals")
      pca.f <- prcomp(f)
      var.f <- round(pca.f$sdev^2/sum(pca.f$sdev^2)*100,2)
      dr <- scale(sqrt(diag(tcrossprod(center(r)))))
      plot_QQ(r)
      plot(pca.f$x[,1], dr, pch=19, 
           xlab = paste("PC 1", var.f[1],"%"),
           ylab = "Standardized Euclidean Distance Residuals",
           main = "Residuals vs. PC 1 fitted")
      abline(h = 0, col = "gray", lty =3)
      if(length(unique(round(pca.f$x[,1], 7))) <= 2) {
        lfr <- list()
        lfr$x <- pca.f$x[,1]
        lfr$y <- dr
        fit <- lm(dr ~ pca.f$x[,1])
        lfr$fitted <- fit$fitted.values
        lfr$residuals <- fit$residuals
      } else {
        options(warn = -1)
        lfr <- loess(dr~pca.f$x[,1], span = 1)
        options(warn = 0)
      } 
      lfr <- cbind(lfr$x, lfr$fitted); lfr <- lfr[order(lfr[,1]),]
      points(lfr, type="l", col="red")
      plot_het(r,f)
      p <- ncol(r)
    }
    
  }
    
  if(type == "regression"){
    PL <- prcomp(f)$x[,1]
    reg.type <- match.arg(reg.type)
    if(is.na(match(reg.type, c("PredLine", "RegScore")))) 
      if(is.null(predictor))
        stop("This plot type is not available without a predictor.")
    n <- NROW(r); p <- NCOL(r)
    if(!is.vector(predictor)) stop("Predictor must be a vector")
    if(length(predictor) != n) 
      stop("Observations in predictor must equal observations if procD.lm fit")
    
    plot_args$x <- predictor
    
    plot_args$ylab <- "Regression Score"
    if(is.null(plot_args$xlab)) plot_args$xlab <- deparse(substitute(predictor))

    xc <- predictor
    if(x$LM$gls) {
      if(!is.null(x$LM$weights)) X <- cbind(xc, x$LM$X) * sqrt(x$LM$weights) else
        X <- cbind(xc, x$LM$Pcov %*% x$LM$X)
    } else X <- cbind(xc, x$LM$X)
    
    X <- as.matrix(X)
    b <- as.matrix(lm.fit(X, f)$coefficients)[1, ]
    Reg.proj <- center(x$LM$Y) %*% b %*% sqrt(solve(crossprod(b)))
    plot_args$y <- Reg.proj
    if(reg.type == "RegScore") {
      do.call(plot, plot_args)
    } else {
      plot_args$y <- PL
      plot_args$ylab <- "PC 1 for fitted values"
      do.call(plot, plot_args)
    }
  }
  if(type == "PC"){
    pca <- prcomp(f)
    eigs <- pca$rotation
    P <- center(x$LM$Y)%*%eigs
    v <- pca$sdev^2
    ev <- round(v[1:2]/sum(v)*100, 2)
    plot_args <- if(NCOL(P) > 1) list(x = P[,1], y = P[,2],  
                        xlab = paste("PC 1 for fitted values: ",ev[1],"%", sep = "") ,
                        ylab = paste("PC 2 for fitted values: ",ev[2],"%", sep = "") ,
                        ...) else list(x = 1:length(P), y = P, 
                                       xlab = "Index", 
                                       ylab = "PC 1 for fitted values: 100%",...)
    do.call(plot, plot_args)
    PC.points <- P
    rownames(P) <- rownames(x$LM$data)
  }
  out <- list(PredLine = PL, RegScore = Reg.proj, PC.points = PC.points, 
              plot_args = plot_args)
  invisible(out)
}

plot_het <- function(r,f){
  r <- center(r)
  f <- center(f)
  r <- sqrt(diag(tcrossprod(r)))
  f <- sqrt(diag(tcrossprod(f)))
  if(length(unique(round(f, 7))) <= 2) {
    lfr <- list()
    lfr$x <- f
    lfr$y <- scale(r)
    fit <- lm(scale(r) ~ f)
    lfr$fitted <- fit$fitted.values
    lfr$residuals <- fit$residuals
  } else {
    options(warn = -1)
    lfr <- loess(scale(r)~f, span = 1)
    options(warn = 0)
  }
  lfr <- cbind(lfr$x, lfr$y, lfr$fitted)
  lfr <- lfr[order(lfr[,1]),]
  plot(lfr, pch=19,  
       xlab = "Euclidean Distance Fitted Values",
       ylab = "Standardized Euclidean Distance Residuals", 
       main = "Residuals vs. Fitted")
  abline(h = 0, col = "gray", lty =3)
  points(lfr[,1], lfr[,3], type="l", col="red")
}

plot_QQ <- function(r){
  r <- center(r)
  r <- sqrt(diag(tcrossprod(r)))
  r <- sort(r)
  n <- length(r)
  tq <- (seq(1,n)-0.5)/n
  tq <- qnorm(tq)
  plot(tq, r, pch=19, xlab = "Theoretical Quantiles",
       ylab = "Euclidean Distance Residuals", 
       main = "Q-Q plot")
}

#' Plot Function for RRPP
#' 
#' @param x plot object (from \code{\link{predict.lm.rrpp}})
#' @param PC A logical argument for whether the data space should be 
#' rotated to its 
#' principal components
#' @param ellipse A logical argument to change error bars to ellipses 
#' in multivariate plots.  
#' It has no function for univariate plots or is abscissa is not NULL.
#' @param abscissa An optional vector (numeric of factor) equal in length 
#' to predictions to use for 
#' plotting as the abscissa (x-axis), in which case predictions are the 
#' ordinate (y-axis).  This might be 
#' helpful if predictions are made for a continuous independent variable.  
#' The abscissa would be the
#' same variable used to make predictions (and can be the data.frame used 
#' for 
#' newdata in \code{\link{predict.lm.rrpp}}).
#' @param label A logical argument for whether points should be labeled 
#' (in multivariate plots).
#' @param ... other arguments passed to plot, arrows, points, or text (helpful 
#' to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}},  \code{\link{arrows}}, \code{\link{points}},
#' \code{\link{par}}, and \code{\link{text}}
#' @method plot predict.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @examples 
#' # See \code{\link{lm.rrpp}} for examples.
plot.predict.lm.rrpp <- function(x, PC = FALSE, ellipse = FALSE,
                                 abscissa = NULL,
                                 label = TRUE, ...){
  oldw <- getOption("warn")
  options(warn = -1)
  
  m <- if(is.matrix(x$mean)) rbind(x$mean, x$lcl, x$ucl) else 
    c(x$mean, x$lcl, x$ucl)
  mpc <- if(is.matrix(x$pc.mean)) rbind(x$pc.mean, x$pc.lcl, x$pc.ucl) else 
    c(x$pc.mean, x$pc.lcl, x$pc.ucl)
  conf <- x$confidence
  k <- NROW(x$mean)
  plot_args <- dots <- list(...)
  arrow_args <- text_args <- eP <- NULL
  
  plot_names <- names(plot_args)
  arrows.names <- c("angle", "length", "code", "col",
                    "lty", "lwd")
  text.names <- c("adj","offset", "pos", "vfont", "labels", "cex")
  
  plot_args <- plot_args[!(plot_names %in% arrows.names)]
  arrow_args <- dots[names(dots) %in% arrows.names]
  arrow_args <- dots[names(dots) %in% text.names]
  
  absx <- FALSE
  
  if(!is.null(abscissa)) {
    absx <- TRUE
    if(is.list(abscissa)) {
      xlabel <- names(abscissa)[[1]]
      abscissa <- unlist(abscissa)
    } else xlabel <- deparse(substitute(abscissa))
    if(length(abscissa) != k)
      stop("\n The length of the abscissa does not match the number of 
           predictions\n",
           call. = FALSE)
  } else xlabel <-  "Predicted values"
  
  plot_type <- if(NCOL(mpc) == 1) "uni" else "multi"
  if(absx) plot_type <- "uni"
  response.type <- if(NCOL(mpc) == 1) "uni" else "multi"
  
  if(plot_type == "uni") {

    if(response.type == "uni") {
      
      xx <- seq(1:k)
      xf <- as.factor(rownames(x$mean))
      if(absx && is.numeric(abscissa)) xx <- abscissa
      
      resp <- if(PC) x$pc.mean[,1] else x$mean[,1]
      mt <- if(PC) paste("Predicted PC1 Values, +/-", 
                         conf*100, "percent confidence levels",
                         sep = " ") else
                           paste("Predicted Values, +/-", 
                                 conf*100, "percent confidence levels",
                                 sep = " ")
      lcl <- if(PC) x$pc.lcl[,1] else x$lcl[,1]
      ucl <- if(PC) x$pc.ucl[,1] else x$ucl[,1]
      
      plot_args$x <- xx
      plot_args$y <- resp
      if(is.null(plot_args$xlab)) plot_args$xlab <- xlabel
      if(is.null(plot_args$ylab)) plot_args$ylab <- if(PC) "PC1: 100%" else  
        x$data.name
      if(is.null(plot_args$ylim)) plot_args$ylim <- c(min(lcl), max(ucl))
      if(is.null(plot_args$main)) plot_args$main <- mt
      if(is.null(plot_args$pch)) plot_args$pch <- 19
      if(is.null(plot_args$cex.main)) plot_args$cex.main <- 0.7
      plot_args$xaxt <- "n"
      
      do.call(plot, plot_args)

      if(absx && is.numeric(abscissa)) axis(1, xx) else
        axis(1, at = xx, labels = as.character(xf)) 
      
      arrow_args$x0 <- xx
      arrow_args$y0 <- resp
      arrow_args$x1 <- xx
      arrow_args$y1 <- lcl
      if(is.null(arrow_args$angle)) arrow_args$angle <- 90
      if(is.null(arrow_args$lemgth)) arrow_args$length <- 0.1
   
      do.call(arrows, arrow_args)
      arrow_args$y1 <- ucl
      do.call(arrows, arrow_args)
      do.call(points, plot_args)
      
    }
    
    if(response.type == "multi") {
      
      if(PC){
        pca <- x$pca
        mr <- mpc
        if(is.matrix(mr)) mr <- mr[,1]
        v <- pca$sdev^2
        ve <- round(v/sum(v)*100, 2)[1:2]
        ylm <- c(min(mr), max(mr))
        ylb <- paste("PC1: ", ve, "%", sep = "")
        type <- "uni"
      } else {
        mr <- m[,1]
        cn <- colnames(x$mean)[[1]]
        if(is.null(cn)) cn <- paste("V", 1, sep=".")
        ylb <- cn
        ylm <- c(min(mr), max(mr))
      }
      
      mt <- if(PC) 
        paste("PC1 predicted Values, +/-", conf*100, 
              "percent confidence levels",
                         sep = " ") else 
        paste("Predicted Values (variable 1), +/-", conf*100, 
              "percent confidence levels",
                        sep = " ")
      xx <- seq(1:k)
      xf <- as.factor(rownames(x$mean))
      if(absx && is.numeric(abscissa)) xx <- abscissa
      
      plot_args$x <- xx
      plot_args$y <- mr[1:k]
      if(is.null(plot_args$xlab)) plot_args$xlab <- xlabel
      if(is.null(plot_args$ylab)) plot_args$ylab <- ylb
      if(is.null(plot_args$ylim)) plot_args$ylim <- c(min(mr), max(mr))
      if(is.null(plot_args$main)) plot_args$main <- mt
      if(is.null(plot_args$pch)) plot_args$pch <- 19
      if(is.null(plot_args$cex.main)) plot_args$cex.main <- 0.7
      plot_args$xaxt <- "n"
      
      do.call(plot, plot_args)
      
      if(absx && is.numeric(abscissa)) axis(1, xx) else
        axis(1, at = xx, labels = as.character(xf)) 
      
      la <- (k + 1):(2 * k)
      ua <- (2 * k +1):(3 * k)
      
      arrow_args$x0 <- xx
      arrow_args$y0 <- mr[1:k]
      arrow_args$x1 <- xx
      arrow_args$y1 <- mr[la]
      if(is.null(arrow_args$angle)) arrow_args$angle <- 90
      if(is.null(arrow_args$lemgth)) arrow_args$length <- 0.1
      
      do.call(arrows, arrow_args)
      arrow_args$y1 <- mr[ua]
      do.call(arrows, arrow_args)
      do.call(points, plot_args)
      
    }
    
  }
  
  if(plot_type == "multi") {
    
    if(PC) {
      
      pca <- x$pca
      mr <- mpc[,1:2]
      v <- pca$sdev^2
      ve <- round(v/sum(v)*100, 2)[1:2]
      xlb <- paste("PC1: ", ve[1], "%", sep = "")
      ylb <- paste("PC2: ", ve[2], "%", sep = "")
      
    }
    
    if(!PC) {
      
      mr <- m[,1:2]
      cn <- colnames(x$mean)
      if(is.null(cn)) cn <- paste("V", 1:2, sep=".")
      xlb <- cn[1]
      ylb <- cn[2]
    }
    
    mt <- if(PC) "Among-prediction PC rotation" else 
      "Plot of first two variables"
    mt <- paste(mt, "; ", conf*100, "% confidence limits", sep = "")
  
    eP <- if(PC) ellipse.points(m = x$pc.mean[,1:2],
                                pr = x$random.predicted.pc, conf) else
                                       ellipse.points(m = x$mean[,1:2],
                                                pr = x$random.predicted, conf)
    span <- apply(eP$ellP, 2, "c")
    xlim <- c(min(span[,1]), max(span[,1]))
    ylim <- c(min(span[,2]), max(span[,2]))
    
    plot_args$x <- mr[1:k, 1]
    plot_args$y <- mr[1:k, 2]
    if(is.null(plot_args$xlab)) plot_args$xlab <- xlb
    if(is.null(plot_args$ylab)) plot_args$ylab <- ylb
    if(is.null(plot_args$xlim)) plot_args$xlim <- xlim
    if(is.null(plot_args$ylim)) plot_args$ylim <- ylim
    if(is.null(plot_args$main)) plot_args$main <- mt
    if(is.null(plot_args$cex.main)) plot_args$cex.main <- 0.7
    if(is.null(plot_args$asp)) plot_args$asp <- 1
    
    do.call(plot, plot_args)
    
    if(ellipse) {
      for(i in 1:(dim(eP$ellP)[[3]])){
        points(eP$ellP[,,i], type = "l", ...)
      }
      if(length(plot_args) == 0) points(eP$means, pch=19, cex = 0.7) else
        points(eP$means, ...)
      
    } else {
      
      la <- (k + 1):(2 * k)
      ua <- (2 * k +1):(3 * k)
      
      if(is.null(arrow_args$angle)) arrow_args$angle <- 90
      if(is.null(arrow_args$lemgth)) arrow_args$length <- 0.05

      arrow_args$x0 <- mr[1:k,1]
      arrow_args$y0 <- mr[1:k,2]
      arrow_args$x1 <- mr[1:k,1]
      arrow_args$y1 <- mr[la, 2]
      
      do.call(arrows, arrow_args)
      
      arrow_args$y1 <- mr[ua, 2]
      do.call(arrows, arrow_args)
      
      arrow_args$x1 <- mr[la, 1]
      arrow_args$y1 = mr[1:k,2]
      do.call(arrows, arrow_args)
      
      arrow_args$x1 <- mr[ua, 1]
      do.call(arrows, arrow_args)
      
      do.call(points, plot_args)
      
      arrow_args$y1 <- mr[la, 2] # return to original
    }
    
    
    if(label) {
      
      text_args <- list(x = NULL, y = NULL, cex = 1, col = 1,
                        adj = NULL, offset = 0.5, pos = 1, vfont = NULL,
                        labels = rownames(x$mean))
      
      dot.match <- intersect(c("x", "y", "cex", "col", "adj", 
                               "offset", "pos", "vfont", "labels"), 
                             names(dots))
      if(length(dot.match) > 0)
        text_args[dot.match] <- dots[dot.match]
      
      text_args$x <- plot_args$x
      text_args$y <- plot_args$y
      
      do.call(text, text_args)
      
    }
    
  }

  options(warn = oldw)
  out <- list(plot_args = plot_args, arrow_args = arrow_args, 
           text_args = text_args, ellipse.points = eP)
  class(out) <- "plot.predict.lm.rrpp"
  invisible(out)
}

## resid and fitted

# residuals.lm.rrpp
# S3 generic for lm.rrpp

#' Extract residuals
#' 
#' @param object plot object (from \code{\link{lm.rrpp}})
#' @param ... Arguments passed to other functions 
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' # See examples for lm.rrpp
residuals.lm.rrpp <- function(object, ...) {
  if(object$LM$gls) out <- object$LM$gls.residuals else 
    out <- object$LM$residuals
  out
}

# fitted.lm.rrpp
# S3 generic for lm.rrpp

#' Extract fitted values
#' 
#' @param object plot object (from \code{\link{lm.rrpp}})
#' @param ... Arguments passed to other functions 
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' # See examples for lm.rrpp
fitted.lm.rrpp <- function(object,  ...) {
  if(object$LM$gls) out <- object$LM$gls.fitted else 
    out <- object$LM$fitted
  out
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{pairwise}}
#' @param ... Other arguments passed onto pairwise
#' @method print pairwise
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.pairwise <- function(x, ...){
  if(is.null(x$LS.means)) type = "slopes"
  if(is.null(x$slopes)) type = "means"
  if(type == "means") {
    means <- x$LS.means
    perms <- length(means)
    groups <- rownames(means[[1]])
  }
  if(type == "slopes") {
    slopes <- x$slopes
    perms <- length(x$slopes)
    groups <- rownames(slopes[[1]])
  }
  perm.type <- x$PermInfo$perm.method
  cat("\nPairwise comparisons\n")
  cat("\nGroups:", groups, "\n\n")
  cat(paste(perm.type,":", sep=""), perms, "permutations\n")
}


#' Print/Summary Function for RRPP
#'
#' See \code{\link{pairwise}} for further description.  
#' 
#' The following summarize the test that can be performed: 
#' 
#' #' \itemize{
#' \item{\bold{Distance between vectors, "dist"}}{ Vectors for LS means or 
#' slopes originate at the origin and point to some location, having both a 
#' magnitude
#' and direction.  A distance between two vectors is the inner-product of of 
#' the vector difference, i.e., the distance between their endpoints.  For
#' LS means, this distance is the difference between means.  For multivariate 
#' slope vectors, this is the difference in location between estimated change 
#' for the dependent variables, per one-unit change of the covariate considered.  
#' For univariate slopes, this is the absolute difference between slopes.}
#' \item{\bold{Vector correlation, "VC"}}{ If LS mean or slope vectors are 
#' scaled to unit size, the vector correlation is the inner-product of the 
#' scaled vectors.
#' The arccosine (acos) of this value is the angle between vectors, which 
#' can be expressed in radians or degrees.  Vector correlation indicates 
#' the similarity of 
#' vector orientation, independent of vector length.}
#' \item{\bold{Difference in vector lengths, "DL"}}{  If the length of a 
#' vector is an important attribute -- e.g., the amount of multivariate 
#' change per one-unit
#' change in a covariate -- then the absolute value of the difference in 
#' vector lengths is a practical statistic to compare vector lengths.  
#' Let d1 and
#' d2 be the distances (length) of vectors.  Then |d1 - d2| is a statistic 
#' that compares their lengths.}
#' \item{\bold{Variance, "var"}}{  Vectors of residuals from a linear 
#' model indicate can express the distances of observed values from 
#' fitted values.  Mean
#' squared distances of values (variance), by group, can be used to 
#' measure the amount of dispersion around estimated values for groups.  
#' Absolute
#' differences between variances are used as test statistics to compare 
#' mean dispersion of values among groups.  Variance degrees of freedom 
#' equal n, 
#' the group size, rather than n-1, as the purpose is to compare mean 
#' dispersion 
#' in the sample.  (Additionally, tests with one subject in a group 
#' are possible, or at least not a hindrance to the analysis.)}
#' }
#' 
#' The argument, \code{test.type} is used to select one of the tests 
#' above.  See \code{\link{pairwise}} for examples.
#' 
#'  \subsection{Notes for RRPP 0.6.2 and subsequent versions}{ 
#'  In previous versions of pairwise, code{\link{summary.pairwise}} had three 
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
#' @param object Object from \code{\link{pairwise}}
#' @param stat.table Logical argument for whether results should be 
#' returned in one table 
#' (if TRUE) or separate pairwise tables (if FALSE)
#' @param test.type the type of statistic to test.  See below
#' should be used in the test.
#' @param angle.type If test.type = "VC", whether angle results are 
#' expressed in radians or degrees.
#' @param confidence Confidence level to use for upper confidence 
#' limit; default = 0.95 (alpha = 0.05)
#' @param show.vectors Logical value to indicate whether vectors 
#' should be printed.
#' @param ... Other arguments passed onto pairwise
#' @method summary pairwise
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.pairwise <- function(object, stat.table = TRUE, 
                             test.type = c("dist", "VC", "DL", "var"),
                             angle.type = c("rad", "deg"),
                             confidence = 0.95, show.vectors = FALSE, ...){
  test.type <- match.arg(test.type)
  angle.type <- match.arg(angle.type)
  x <- object
  if(test.type != "var") {
    if(is.null(x$LS.means)) type = "slopes"
    if(is.null(x$slopes)) type = "means"
  } else {
    type <- "var"
    tab <- NULL
  }
 
  vars <- object$vars
  if(type == "var") {
    var.diff <- lapply(1:NCOL(vars), function(j){
      v <- as.matrix(vars[,j])
      as.matrix(dist(v))
    })
    L <- d.summary.from.list(var.diff, confidence = confidence)
    tab <- makePWDTable(L)
  }
  
  if(type == "means") {
    if(test.type == "dist") {
      L <- d.summary.from.list(x$means.dist, confidence = confidence)
      if(stat.table) tab <- makePWDTable(L)
    }
    
    if(test.type == "VC") {
      L <- r.summary.from.list(x$means.vec.cor, confidence = confidence)
      if(angle.type == "deg") {
        options(warn = -1)
        L$angle <- L$angle * 180 / pi
        L$aCL <- L$aCL * 180 / pi
        options(warn = 0)
      }
      if(stat.table) tab <- makePWCorTable(L)
    }
    
    if(test.type == "DL") {
      L <- d.summary.from.list(x$means.diff.length, confidence = confidence)
      if(stat.table) tab <- makePWDTable(L)
    }
    
  }
  
  if(type == "slopes") {

    if(test.type == "dist") {
      L <- d.summary.from.list(x$slopes.dist)
      if(stat.table) tab <- makePWDTable(L)
    }
    
    if(test.type == "VC") {
      L <- r.summary.from.list(x$slopes.vec.cor) 
      if(angle.type == "deg") {
        options(warn = -1)
        L$angle <- L$angle * 180 / pi
        L$aCL <- L$aCL * 180 / pi
        options(warn = 0)
      }
      if(stat.table) tab <- makePWCorTable(L)
    }
    
    if(test.type == "DL") {
      L <- d.summary.from.list(x$slopes.diff.length, confidence = confidence)
      if(stat.table) tab <- makePWDTable(L)
    }
    
  }
  
  out <- list()
  out$pairwise.tables <- L
  if(stat.table){
    out$stat.table <- TRUE
    out$summary.table <- tab
  } else {
    out$stat.table <- FALSE
    out$summary.table <- NULL
  }

  out$type <- type
  out$test.type <- test.type
  out$angle.type <- angle.type 
  out$confidence <- confidence
  if(!is.null(vars)) out$vars <- vars else out$vars <- NULL
  out$show.vectors <- show.vectors
  out$x <- x
  
  class(out) <- "summary.pairwise"
  out
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{summary.pairwise}}
#' @param ... Other arguments passed onto summary.pairwise
#' @method print summary.pairwise
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.summary.pairwise <- function(x, ...) {
  type <- x$type
  test.type <- x$test.type
  stat.table <- x$stat.table
  print.pairwise(x$x)
  cat("\n")
  
  L <- x$pairwise.tables
  if(stat.table) tab <- x$summary.table

  if(type == "var") {
    
    cat("\nObserved variances by group\n\n")
    print(x$vars[,1])
    
    if(stat.table) {
      cat("\nPairwise distances between variances, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise distances between variances\n")
      print(L$D)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
          "upper confidence limits between variances\n")
      print(L$CL)
      cat("\nPairwise effect sizes (Z) between variances\n")
      print(L$Z)
      cat("\nPairwise P-values between variances\n")
      print(L$P)
    }
  }
  
  if(type == "means") {
    
    cat("LS means:\n")
    if(x$show.vectors) print(x$x$LS.means[[1]]) else 
      cat("Vectors hidden (use show.vectors = TRUE to view)\n")
    
    if(test.type == "dist") {
      if(stat.table) {
        cat("\nPairwise distances between means, plus statistics\n")
        print(tab)
      } else {
        cat("\nPairwise distances between means\n")
        print(L$D)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
            "Upper confidence limits between means\n")
        print(L$CL)
        cat("\nPairwise effect sizes (Z) between means\n")
        print(L$Z)
        cat("\nPairwise P-values between means\n")
        print(L$P)
      }
    }
    
    if(test.type == "VC") {
      if(stat.table) {
        cat("\nPairwise statistics based on mean vector correlations\n")
        print(tab)
      } else {
        cat("\nPairwise vector correlations between mean vectors\n")
        print(L$r)
        cat("\nPairwise angles between mean vectors\n")
        print(L$angle)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
            "Upper confidence limits for angles between mean vectors\n")
        print(L$aCL)
        cat("\nPairwise effect sizes (Z) for angles between mean vectors\n")
        print(L$Z)
        cat("\nPairwise P-values for angles between mean vectors\n")
        print(L$P)
      }
    }
    
    if(test.type == "DL") {
      if(stat.table) {
        cat("\nPairwise absolute difference (d) between vector lengths, 
            plus statistics\n")
        print(tab)
      } else {
        cat("\nPairwise absolute differences (d) between mean vector lengths\n")
        print(L$D)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
            "Upper confidence limits between mean vector lengths\n")
        print(L$CL)
        cat("\nPairwise effect sizes (Z) between mean vector lengths\n")
        print(L$Z)
        cat("\nPairwise P-values between mean vector lengths\n")
        print(L$P)
      }
    }
    
  }
  
  if(type == "slopes") {
    cat("Slopes (vectors of variate change per one unit of covariate 
        change, by group):\n")
    if(x$show.vectors) print(x$x$slopes[[1]]) else 
      cat("Vectors hidden (use show.vectors = TRUE to view)\n")
    
    if(test.type == "dist") {
      if(stat.table) {
        cat("\nPairwise distances between slope vector 
            (end-points), plus statistics\n")
        print(tab)
      } else {
        cat("\nPairwise distances between slope vector (end-points\n")
        print(L$D)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
            "Upper confidence limits between slopes\n")
        print(L$CL)
        cat("\nPairwise effect sizes (Z) between slopes\n")
        print(L$Z)
        cat("\nPairwise P-values between slopes\n")
        print(L$P)
      }
    }
    
    if(test.type == "VC") {
      cat("\nPairwise statistics based on slopes vector correlations (r) 
          and angles, acos(r)")
      cat("\nThe null hypothesis is that r = 1 (parallel vectors).")
      cat("\nThis null hypothesis is better treated as the angle 
          between vectors = 0\n")
      if(stat.table) print(tab)
        else {
        cat("\nPairwise vector correlations between slope vectors\n")
        print(L$r)
        cat("\nPairwise angles between slope vectors\n")
        print(L$angle)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
            "upper confidence limits for angles between mean vectors\n")
        print(L$aCL)
        cat("\nPairwise effect sizes (Z) for angles between slope vectors\n")
        print(L$Z)
        cat("\nPairwise P-values for angles between slope vectors\n")
        print(L$P)
      }
    }
    
    if(test.type == "DL") {
      cat("\nSlope vector lengths\n")
      print(x$x$slopes.length[[1]])
      if(stat.table) {
        cat("\nPairwise absolute difference (d) between vector 
            lengths, plus statistics\n")
        print(tab)
      } else {
        cat("\nPairwise absolute differences (d) between slope vector lengths\n")
        print(L$D)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
            "Upper confidence limits between slope vector lengths\n")
        print(L$CL)
        cat("\nPairwise effect sizes (Z) between slope vector lengths\n")
        print(L$Z)
        cat("\nPairwise P-values between slope vector lengths\n")
        print(L$P)
      }
    }
    
  }
  invisible(x)
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{model.comparison}}
#' @param ... Other arguments passed onto model.comparison
#' @method print model.comparison
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.model.comparison <- function(x ,...){
  print(x$table)
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{model.comparison}}
#' @param ... Other arguments passed onto model.comparison
#' @method summary model.comparison
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.model.comparison <- function(object ,...){
  x <- object
  type <- colnames(x$table)[[1]]
  if(type == "cov.trace") type.name <- "traces of model covariance matrices"
  if(type == "logLik") type.name <- "model log-likelihoods"
  if(type == "Z") type.name <- "model Z-scores"
  n <- NROW(x$table)
  cat("\n\n Summary statistics for", type.name, "\n\n")
  cat(n, "Models compared.\n\n")
  tab <- x$table
  rownames(tab) <- x$names
  print(tab)
}

#' Plot Function for RRPP
#' 
#' @param x plot object (from \code{\link{model.comparison}})
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' @method plot model.comparison
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.model.comparison <- function(x, ...){
  object <- x
  type <- names(object$table)[[1]]
  tab <- object$table
  nms <- object$names
  spln <- if(type == "logLik") ncol(tab) == 5 else 
    if(type == "cov.trace") ncol(tab) == 3 else ncol(tab) == 2
  
  if(type == "cov.trace") type.name <- "Trace"
  if(type == "logLik") {
    type.name <- if(spln) "log-likelihood" else "-2 * log-likelihood"
  }
  if(type == "Z") type.name <- "Z"
  
  if(type == "logLik") {
    x <- if(spln) tab[,5] else tab[,3]
    xLab <- if(spln) names(tab)[5] else names(tab)[3] 
  } else if(type == "cov.trace"){
    x <- if(spln) tab[,3] else tab[,2]
    xLab <- if(spln) names(tab)[3] else names(tab)[2]
  } else {
    x <- if(spln) tab[,2] else NULL
    xLab <- names(tab)[if(spln) 2 else NULL]
  }
    
  if(type == "Z" && is.null(x)){
    cat("No predictor.  Z-scores will be plotted in the order provided.\n")
    x <- 1:nrow(tab)
  }
  y <- tab[,1]
  if(all(is.na(x))) x <- 1:length(y)
  if(type == "logLik" && !spln) y <- -2 * y
  xrange <- range(x)
  dx <- xrange[2] - xrange[1]
  xlim <- xrange + c(-0.1*dx, 0.1*dx)
  
  yrange <- range(y)
  dy <- yrange[2] - yrange[1]
  ylim <- yrange + c(-0.1*dy, 0.1*dy)

  
  if(var(x) == 0)
    stop("There is no plot to make, as the parameter penalty is invariant.",
         call. = FALSE) 
  plot(x, y, xlab = xLab, 
       ylab = type.name, xlim = xlim, ylim = ylim, ...)
  
  if(spln){
    
    if(!any(is.na(x))) {
      sp <- spline(x, y)
      points(sp$x, sp$y, col = 2, type = "l")
      xmax <- sp$x[which.max(sp$y)]
      xmin <-sp$x[which.min(sp$y)]
      text(x, y, nms, pos = 1, cex = 0.4)
      cat("max", type, "found at", xLab, "=", round(xmax, 4))
      cat("\nmin", type, "found at", xLab, "=", round(xmin, 4), "\n")
    }
      
  } else {
    if(is.null(xLab))
      stop("There is no plot to make, as there is no predictor.",
           call. = FALSE) 
    f <- lm(y ~ x)
    abline(f, lty = 3, lwd = 0.8, col = "red")
    text(x, y, nms, pos = 1, cex = 0.4)
               
           
  }

}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{lm.rrpp}}, updated with 
#' \code{\link{manova.update}}
#' @param test Type of multivariate test statistic to use.
#' @param ... Other arguments passed onto manova.lm.rrpp
#' @method summary manova.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.manova.lm.rrpp <- function(object, test = c("Roy", "Pillai", "Hotelling-Lawley", "Wilks"), ...){
  if (!inherits(object, "manova.lm.rrpp")) 
    stop(gettextf("object must be of class %s", dQuote("manova.lm.rrpp"), domain = NA))
  test <- match.arg(test)
  if(test == "Hotelling-Lawley") test <- "Hotelling.Lawley"
  p <- object$LM$p
  p.prime <- object$LM$p.prime
  n <- object$LM$n
  Models <- getModels(object, attribute = "all")
  PermInfo <- getPermInfo(object, attribute = "all")
  perm.method <- object$PermInfo$perm.method
  if(perm.method == "RRPP") RRPP = TRUE else RRPP = FALSE
  ind <- PermInfo$perm.schedule
  perms <- length(ind)
  trms <- object$LM$term.labels
  k <- length(trms)
  kk <- length(Models$full)
  if(k > kk){
    k <- kk
    trms <- names(Models$full)
  }
  df <- object$ANOVA$df
  df.model <- sum(df[1:k])
  df <- c(df[1:k], df.model, df[k+1])
  names(df) <- c(trms, "Full.Model", "Residuals")
  
  MANOVA <- object$MANOVA
  eigs <- MANOVA$eigs
  error <- MANOVA$error
  stats <- as.data.frame(matrix(NA, nrow = k + 2, ncol = 5, byrow = FALSE,
                                dimnames <- list(names(df), 
                                    c("Df", "Rand", test, "Z", "Pr"))))
  stats$Df <- df
  if(!is.null(error)) stats$Rand[1:(k+1)] <- c(error, 
                            "Residuals") else 
                              stats$Rand[1:(k+1)] <- rep("Residuals", k+1)
  
  if(test == "Pillai") {
    
    rand.stats <- sapply(1:perms, function(j){
      y <- eigs[[j]]
      apply(y, 1, pillai)
    })
    
    test.stats <- rand.stats[, 1]
    Z <- apply(rand.stats, 1, effect.size)
    P <- apply(rand.stats, 1, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Pillai[1:(k+1)] <- test.stats
  }
  else if(test == "Hotelling.Lawley") {
    
    rand.stats <- sapply(1:perms, function(j){
      y <- eigs[[j]]
      apply(y, 1, hot.law)
    })
    
    test.stats <- rand.stats[, 1]
    Z <- apply(rand.stats, 1, effect.size)
    P <- apply(rand.stats, 1, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Hotelling.Lawley[1:(k+1)] <- test.stats
  }
  
  else if(test == "Wilks"){
    
    rand.stats <- sapply(1:perms, function(j){
      y <- eigs[[j]]
      apply(y, 1, wilks)
    })
    
    test.stats <- rand.stats[, 1]
    Z <- apply(-1 * rand.stats, 1, effect.size)
    P <- apply(-1 * rand.stats, 1, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Wilks[1:(k+1)] <- test.stats
  }
  else {
    rand.stats <- sapply(1:perms, function(j){
      y <- eigs[[j]]
      apply(y, 1, max)
    })
    
    test.stats <- rand.stats[,1]
    Z <- apply(rand.stats, 1, effect.size)
    P <- apply(rand.stats, 1, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Roy[1:(k+1)] <- test.stats
  }
  
  if(test == "Wilks") 
    names(stats)[[length(stats)]] <- paste("Pr(<", test, ")", sep = "") else
    names(stats)[[length(stats)]] <- paste("Pr(>", test, ")", sep = "")
  
  out <- list(stats.table = stats, rand.stats = rand.stats, stat.type = test,
              n = n, p = p, p.prime = p.prime, e.rank = MANOVA$e.rank, 
              manova.pc.dims = MANOVA$manova.pc.dims, PCA = MANOVA$PCA,
              SS.type = object$ANOVA$SS.type, SS.tot = MANOVA$SS.tot, 
              perms = object$PermInfo$perms)
  class(out) <- "summary.manova.lm.rrpp"
  out
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{summary.manova.lm.rrpp}}
#' @param ... Other arguments passed onto summary.manova.lm.rrpp
#' @method print summary.manova.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.summary.manova.lm.rrpp <- function(x, ...){
  
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of observations:", x$n))
  cat("\nNumber of dependent variables:", x$p)
  cat("\nData space dimensions:", x$p.prime)
  cat("\nResidual covariance matrix rank:", x$e.rank)
  pc.max <- x$manova.pc.dims
  
  if(pc.max < x$e.rank) {
    cat("\n   Data reduced to", pc.max, "PCs, as required or prescribed.")
    PCA <- x$PCA
    d2 <- PCA$sdev^2
    d.p <- sum(d2[1:pc.max])/x$SS.tot
    cat("\n  ", round(d.p*100, 1), 
        "% of overall variation explained by these PCs.")
    cat("\n   See $MANOVA$PCA from manova.lm.rrpp object for more information.")
  }
  
  if(!is.null(x$SS.type)) 
    cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
  cat(paste("\nNumber of permutations:", x$perms), "\n\n")
  
  tab <- as.matrix(x$stats.table)
  class(tab) <- c("anova", class(tab))
  print.table(tab, na.print = "")
  invisible(x)
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{lr_test}}
#' @param ... Other arguments passed onto print.lr_test
#' @method print lr_test
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.lr_test <- function(x, ...){
  x <- if(!is.null(x$tab)) x$tab else as.data.frame(x)
  colnames(x)[4] <- "Pr(>LR)"
  class(x) <- c("anova", class(x))
  
  cat("\n\n")
  
  print(x)
  cat("\n\n")
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{lr_test}}
#' @param ... Other arguments passed onto print.lr_test
#' @method summary lr_test
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.lr_test <- function(object, ...){
  print.lr_test(object)
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{trajectory.analysis}}
#' @param ... Other arguments passed onto 
#' @method print trajectory.analysis
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.trajectory.analysis<- function(x, ...){
  if(is.null(x$SD)) type = "vectors" else type = "trajectories"
  perms <- length(x$MD)
  pca = x$pca
  groups <- rownames(x$LS.means[[1]])
  pca <- x$pca
  cat("\nTrajectory analysis\n\n")
  cat(perms, "permutations.\n\n")
  if(!is.null(pca)) cat("Points projected onto trajectory PCs\n")
  
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{trajectory.analysis}}
#' @param stat.table Logical argument for whether results should be 
#' returned in one table 
#' (if TRUE) or separate pairwise tables (if FALSE)
#' @param attribute Whether magnitude differences (MD, absolute difference 
#' in trajectory path lengths), 
#' trajectory correlations (TC), or trajectory shape differences (SD) are 
#' summarized.
#' @param angle.type If attribute = "TC", whether angle results are 
#' expressed in radians or degrees.
#' @param confidence Confidence level to use for upper confidence limit; 
#' default = 0.95 (alpha = 0.05)
#' @param show.trajectories Logical value to indicate whether trajectories 
#' should be printed.
#' @param ... Other arguments passed onto trajectory.analysis
#' @method summary trajectory.analysis
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.trajectory.analysis <- function(object, stat.table = TRUE, 
                             attribute = c("MD", "TC", "SD"),
                             angle.type = c("rad", "deg"),
                             confidence = 0.95, show.trajectories = FALSE, ...) {
  
  attribute <- match.arg(attribute)
  angle.type <- match.arg(angle.type)
  x <- object
  if(is.null(x$SD)) type = "vectors" else type = "trajectories"
  
  MD <- object$MD
  TC <- object$TC
  SD <- object$SD
  
  if(attribute == "MD"){
    L <- d.summary.from.list(MD, confidence = confidence)
    tab <- makePWDTable(L)
  }
  
  if(attribute == "TC"){
    if(is.null(TC)) stop("Trajectory correlations not available\n", 
                         call. = FALSE)
    L <- r.summary.from.list(TC, confidence = confidence)
    
    if(angle.type == "deg") {
      options(warn = -1)
      L$angle <- L$angle * 180 / pi
      L$aCL <- L$aCL * 180 / pi
      options(warn = 0)
    }
    
    tab <- makePWCorTable(L)
  }
  
  if(attribute == "SD"){
    if(type == "trajectories") {
      L <- d.summary.from.list(SD, confidence = confidence)
      tab <- makePWDTable(L)
    } else {
      L <- NULL
      tab <- NULL
    }
  }
  
 
  out <- list()
  out$pairwise.tables <- L
  if(stat.table){
    out$stat.table <- TRUE
    out$summary.table <- tab
  } else {
    out$stat.table <- FALSE
    out$summary.table <- NULL
  }
  
  out$type <- type
  out$attribute <- attribute
  out$angle.type <- angle.type 
  out$confidence <- confidence
  out$show.trajectories <- show.trajectories
  out$x <- x
  
  class(out) <- "summary.trajectory.analysis"
  out
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{summary.trajectory.analysis}}
#' @param ... Other arguments passed onto summary.trajectory.analysis
#' @method print summary.trajectory.analysis
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.summary.trajectory.analysis <- function(x, ...) {
  attribute <- x$attribute
  type <- x$type
  if(attribute == "SD" && type == "vectors") 
    cat("\n\nCannot summarize trajectory shape differences for vectors.\n\n")
  stat.table <- x$stat.table
  if(attribute != "SD") print.trajectory.analysis(x$x) else 
    if(type == "trajectories") print.trajectory.analysis(x$x)

  cat("\n")
  
  L <- x$pairwise.tables
  if(stat.table) tab <- x$summary.table
  
  if(attribute == "MD") {
    
    cat("Trajectories:\n")
    if(x$show.trajectories) print(x$x$trajectories[[1]]) else 
      cat("Trajectories hidden (use show.trajectories = TRUE to view)\n")
    
    cat("\nObserved path distances by group\n\n")
    print(x$x$PD[[1]])
    
    if(stat.table) {
      cat("\nPairwise absolute differences in path distances, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise absolute differences in path distancess\n")
      print(L$D)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
          "upper confidence limits, absolute differences in path distancess\n")
      print(L$CL)
      cat("\nPairwise effect sizes (Z) for absolute differences in path distancess\n")
      print(L$Z)
      cat("\nPairwise P-values for absolute differences in path distances\n")
      print(L$P)
    }
  }
  
  if(attribute == "TC") {
    
    cat("Trajectories:\n")
    if(x$show.trajectories) print(x$x$trajectories[[1]]) else 
      cat("Trajectories hidden (use show.trajectories = TRUE to view)\n")
    
    if(stat.table) {
      cat("\nPairwise correlations between trajectories, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise statistics based on trajectory vector 
          correlations (r) and angles, acos(r)")
      cat("\nThe null hypothesis is that r = 1 (parallel vectors).")
      cat("\nThis null hypothesis is better treated as the angle 
          between vectors = 0\n")
      
      cat("\nPairwise vector correlations between trajectories\n")
      print(L$r)
      cat("\nPairwise angles between trajectories\n")
      print(L$angle)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
          "upper confidence limits for angles\n")
      print(L$aCL)
      cat("\nPairwise effect sizes (Z) for angles\n")
      print(L$Z)
      cat("\nPairwise P-values for angles\n")
      print(L$P)
    }
  }
  
  if(type == "trajectories"  && attribute == "SD") {
    
    cat("Trajectories:\n")
    if(x$show.trajectories) print(x$x$trajectories[[1]]) else 
      cat("Trajectories hidden (use show.trajectories = TRUE to view)\n")

    if(stat.table) {
      cat("\nPairwise trajectory shape differences, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise trajectory shape differences\n")
      print(L$D)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), 
          "upper confidence limits, shape differnces\n")
      print(L$CL)
      cat("\nPairwise effect sizes (Z) for trajectory shape differencess\n")
      print(L$Z)
      cat("\nPairwise P-values for trajectory shape differences\n")
      print(L$P)
    }
  }
  
  cat("\n\n")
    
  invisible(x)
}


#' Plot Function for RRPP
#' 
#'  Function generates a principal component plot for trajectories
#'
#'  The function calculates and plots principal components of fitted values from 
#'  \code{\link{lm.rrpp}} that are passed onto \code{\link{trajectory.analysis}}, 
#'  and projects
#'  data onto them.  This function is a set.up, and \code{\link{add.trajectories}} 
#'  is needed to 
#'  add trajectories to the plot.  By having two stages of control, the plotting 
#'  functions are more 
#'  flexible.  This function also returns plotting information that can be 
#'  valuable for making
#'  individualized plots, if \code{\link{add.trajectories}} is not preferred.
#' @param x plot object (from \code{\link{trajectory.analysis}})
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' 
#' @return If an object is assigned, it will return:
#'  \item{pca}{Principal component analysis performed using \code{\link{prcomp}}.}
#'  \item{pc.points}{Principal component scores for all data.}
#'  \item{trajectory.analysis}{Trajectory analysis passed on.}
#'  \item{trajectories}{pca Observed trajectories projected onto principal 
#'  components.}
#'  
#' @seealso 
#' \code{\link{plot.default}} and \code{\link{par}}
#' @method plot trajectory.analysis
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Adams, D. C., and M. M. Cerney. 2007. 
#' Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. 
#' The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. 
#' A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. 
#' Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @references Collyer, M. L., and D. C. Adams. 2013. 
#' Phenotypic trajectory analysis: comparison of shape change patterns 
#' in evolution and ecology. Hystrix 24: 75-83.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. 
#' A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' 
#' @examples 
#' # See \code{\link{trajectory.analysis}} for examples
plot.trajectory.analysis <- function(x, ...) {
  
  if(!is.null(x$pca)) {
    pca <- x$pca
    rot <- pca$rotation
    Y <- x$fit$LM$Y
    props <- pca$sdev^2 / sum(pca$sdev^2)
    pc.points <- center(Y) %*% rot
    trajectories <- x$trajectories[[1]]
  }
  
  if(is.null(x$pca) && x$type == "factorial") {
    f <- if(x$fit$LM$gls) x$fit$LM$gls.fitted else x$fit$LM$fitted
    pca <- prcomp(f)
    rot <- pca$rotation
    Y <- x$fit$LM$Y
    Y.cent <- colMeans(Y)
    props <- pca$sdev^2 / sum(pca$sdev^2)
    pc.points <- center(Y) %*% rot
    trajectories <- x$trajectories[[1]]
    if(is.matrix(trajectories)) trajectories <- list(trajectories)
    traj.c <- matrix(Y.cent, NROW(trajectories[[1]]), 
                     NCOL(trajectories[[1]]), byrow = TRUE)
    trajectories <- lapply(trajectories, function(x) (x - traj.c) %*% rot)
  }
  
  if(x$type == "single.factor") {
    f <- if(x$fit$LM$gls) x$fit$LM$gls.fitted else x$fit$LM$fitted
    tp <- x$n.points
    p <- NCOL(f)/tp
    n <- NROW(f)
    ft <- array(f, c(n, p, tp))
    ft2 <- ft[,,1]
    for(i in 2:tp) ft2 <- rbind(ft2, ft[,,i])
    pca <- prcomp(ft2)
    rot <- pca$rotation
    Y <- x$fit$LM$Y
    Y2 <- array(Y, c(n, p, tp))
    Y <- Y2[,,1]
    for(i in 2:tp) Y <- rbind(Y, Y2[,,i])
    Y.cent <- colMeans(Y)
    props <- pca$sdev^2 / sum(pca$sdev^2)
    pc.points <- center(Y) %*% rot
    trajectories <- x$trajectories[[1]]
    if(is.matrix(trajectories)) trajectories <- list(trajectories)
    traj.c <- matrix(Y.cent, NROW(trajectories[[1]]), 
                     NCOL(trajectories[[1]]), byrow = TRUE)
    trajectories <- lapply(trajectories, function(x) (x - traj.c) %*% rot)
  }
  
  
  dots <- list(...)
  if(is.null(dots$xlab))
    xlabel <- paste("PC 1 for fitted values: ", 
                    round(props[1] *100, 2), "%", sep = "")
  if(is.null(dots$ylab))
    ylabel <- paste("PC 2 for fitted values: ", 
                    round(props[2] *100, 2), "%", sep = "")
  
  if(!is.null(dots$xlab) && !is.null(dots$ylab)) 
    plot(pc.points[,1], pc.points[,2], asp = 1, ...)
  if(!is.null(dots$xlab) && is.null(dots$ylab)) 
    plot(pc.points[,1], pc.points[,2], asp = 1, ylab = ylabel, ...)
  if(is.null(dots$xlab) && !is.null(dots$ylab)) 
    plot(pc.points[,1], pc.points[,2], asp = 1, xlab = xlabel, ...)
  if(is.null(dots$xlab) && is.null(dots$ylab))
    plot(pc.points[,1], pc.points[,2], asp = 1, xlab = xlabel, ylab = ylabel, ...)
  out <- list(pca = pca, pc.points = pc.points, 
              trajectoy.analysis = x, trajectories = trajectories)
  invisible(out)
}

#' Plot Function for RRPP
#' 
#'  Function adds trajectories to a principal component plot
#'
#'  The function adds trajectories to a plot made by 
#'  \code{\link{plot.trajectory.analysis}}.
#'  This function has a restricted set of plot parameters 
#'  based on the number of trajectories
#'  to be added to the plot.
#'  
#' @param TP plot object (from \code{\link{plot.trajectory.analysis}})
#' @param traj.pch Plotting "character" for trajectory points.  
#' Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for pch.
#' @param traj.col The color of trajectory lines.  
#' Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for col.
#' @param traj.lty Trajectory line type.  Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for lty.
#' @param traj.lwd Trajectory line width.  Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for lwd.
#' @param traj.cex Trajectory point character expansion.  Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for cex.
#' @param traj.bg Trajectory point background.  Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for bg.
#' @param start.bg Trajectory point background, just the start points.  
#' Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for bg.  Green start points are the default.
#' @param end.bg Trajectory point background, just the end points.  
#' Can be a single value or vector 
#' of length equal to the number of trajectories.  
#' See \code{\link{par}} and its description 
#' for bg.  Red end points are the default.
#' 
#' @seealso 
#' \code{\link{plot.default}} and \code{\link{par}}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Adams, D. C., and M. M. Cerney. 2007. 
#' Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. 
#' The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. 
#' A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. 
#' Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @references Collyer, M. L., and D. C. Adams. 2013. 
#' Phenotypic trajectory analysis: comparison of shape change patterns 
#' in evolution and ecology. Hystrix 24: 75-83.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. 
#' A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
add.trajectories <- function(TP, 
                             traj.pch = 21,
                             traj.col = 1,
                             traj.lty = 1,
                             traj.lwd = 1,
                             traj.cex = 1.5,
                             traj.bg = 1,
                             start.bg = 3,
                             end.bg = 2) {
  
  traj <- TP$trajectories
  nt <- length(traj)
  np <- NROW(traj[[1]])
  
  if(length(traj.pch) != 1 && length(traj.pch) != nt)
    stop("For add.trajectories, traj.pch must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.pch) == 1) 
           traj.pch <- rep(traj.pch, nt)
  
  if(length(traj.col) != 1 && length(traj.col) != nt)
    stop("For add.trajectories, traj.col must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.col) == 1) 
           traj.col <- rep(traj.col, nt)
  
  if(length(traj.lty) != 1 && length(traj.lty) != nt)
    stop("For add.trajectories, traj.lty must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.lty) == 1) 
           traj.lty <- rep(traj.lty, nt)
  
  if(length(traj.lwd) != 1 && length(traj.lwd) != nt)
    stop("For add.trajectories, traj.lwd must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.lwd) == 1) 
           traj.lwd <- rep(traj.lwd, nt)
  
  if(length(traj.cex) != 1 && length(traj.cex) != nt)
    stop("For add.trajectories, traj.cex must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.cex) == 1) 
           traj.cex <- rep(traj.cex, nt)
  
  if(length(traj.bg) != 1 && length(traj.bg) != nt)
    stop("For add.trajectories, traj.bg must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.bg) == 1) 
           traj.bg <- rep(traj.bg, nt)
  
  if(length(start.bg) != 1 && length(start.bg) != nt)
    stop("For add.trajectories, start.bg must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(start.bg) == 1) 
           start.bg<- rep(start.bg, nt)
  
  if(length(end.bg) != 1 && length(end.bg) != nt)
    stop("For add.trajectories, end.bg must be equal in length 
         to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(end.bg) == 1) 
           end.bg<- rep(end.bg, nt)     
  
  for(i in 1:nt){
    x <- traj[[i]][,1]
    y <- traj[[i]][,2]
    lines(x, y, col = traj.col[i], lwd = traj.lwd[i], lty = traj.lty[i])
    points(x, y, col = 1, pch = traj.pch[i], 
           lwd = traj.lwd[i], cex = traj.cex[i], bg = traj.bg[i])
    points(x[1], y[1], col = 1, pch = traj.pch[i], 
           lwd = traj.lwd[i], cex = traj.cex[i], bg = start.bg[i])
    points(x[np], y[np], col = 1, pch = traj.pch[i], 
           lwd = traj.lwd[i], cex = traj.cex[i], bg = end.bg[i])
  }
  
}

# ordinate

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{ordinate}}
#' @param ... Other arguments passed onto print.ordinate
#' @method print ordinate
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.ordinate <- function(x, ...){
  ord.type <- if(x$alignment == "principal") 
    "Principal Component Analysis" else
    "Alignment to an alternative matrix"
  cat("\nOrdination type:", ord.type, "\n")
  if(x$alignment != "principal")
    cat("Alignment matrix:", x$alignment, "\n")
  cen <- if(x$GLS) "GLS" else "OLS"
  cat("Centering by", cen, "mean\n")
  
  if(x$alignment == "principal") {
    if(x$GLS && x$transform)
      cat("GLS residuals transformed for orthogonal projection\n") else
        if(x$GLS) cat("Oblique projection of GLS-centered residuals\n") else
          cat("Orthogonal projection of OLS residuals\n")
  } else {
    if(x$GLS && x$transform)
      cat("GLS residuals transformed\n") else
        if(x$GLS) cat("GLS-centered residuals, not transformed\n") else
          cat("OLS residuals\n")
    cat("Alignment to ", x$alignment, 
        "means residual projection is not orthogonal.\n")
  }
  
  cat("Number of observations:", NROW(x$x), "\n")
  cat("Number of vectors", NCOL(x$x), "\n\n")
} 

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{ordinate}}
#' @param ... Other arguments passed onto print.ordinate
#' @method summary ordinate
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
summary.ordinate <- function(object, ...){
  x <- object
  print.ordinate(x, ...)
  d <- x$d
  p <- d/sum(d)
  cp <- cumsum(d)/sum(d)
  r <- as.data.frame(rbind(d, p, cp))
  r <- r[, 1:min(length(d), NCOL(x$x), NCOL(r))]
  colnames(r) <- colnames(x$x)[1:NCOL(r)]
  rownames(r) <- c("Singular Value", 
                   "Proportion of Covariance", "Cumulative Proportion")
  
  if(x$alignment == "principal") rownames(r)[1:2] <- c("Eigenvalues", 
                                                       "Proportion of Variance")
  
  if(x$alignment != "principal") {
    rv <- x$RV
    r <- rbind(r, rv, cumsum(rv))
    rownames(r)[4:5] <- c("RV by Component", "Cumulative RV")
  }

  cat("Importance of Components:\n")
  print(r)
  out <- r
  invisible(out)
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{summary.ordinate}}
#' @param ... Other arguments passed onto print.ordinate
#' @method print summary.ordinate
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.summary.ordinate <- function(x, ...){
  cat("\nImportance of Components:\n")
  print(x$table)
} 

#' Plot Function for RRPP
#' 
#' @param x An object of class \code{\link{ordinate}}
#' @param axis1 A value indicating which component should be 
#' displayed as the X-axis (default = C1)
#' @param axis2 A value indicating which component should be 
#' displayed as the Y-axis (default = C2)
#' @param flip An argument that if not NULL can be used to flip 
#' components in the plot.  
#' The values need to match axis1 or axis2.  For example, if axis1 = 3 
#' and axis2 = 4, flip = 1 will not
#' change either axis; flip = 3 will flip only the horizontal axis; 
#' flip = c(3, 4) will flip both axes.
#' @param include.axes A logical argument for whether axes should be shown at x = 0 and y = 0.  
#' This is different than the axes argument in the generic \code{\link{plot.default}} function, which
#' controls the edges of the plot (providing a box effect or not).  Using include.axes = TRUE does not 
#' allow aesthetic control of the axes.  If desired, it is better to use include.axes = FALSE and augment
#' the plot object with \code{\link{abline}} (choosing h = 0 and v = 0 in separate applications).
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' @return An object of class "plot.ordinate" is a list with components
#'  that can be used in other plot functions, such as the type of plot, points, 
#'  a group factor, and other information depending on the plot parameters used.
#' @method plot ordinate
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.ordinate <- function(x, axis1 = 1, axis2 = 2, flip = NULL, 
                          include.axes = TRUE, ...) {
  wrn <- options()$warn
  options(warn = -1)
  if(NCOL(x$x) == 1) 
    stop("Only one component.  No plotting capability with this function.\n", 
                          call. = FALSE)
  v <- x$d/sum(x$d)
  if(!is.null(x$RV)) rv <- x$RV
  
  pcdata <- x$x[, c(axis1, axis2)]
  if(!is.null(flip)) {
    if(length(flip) > 2) flip <- flip[1:2]
    flip <- flip[(flip %in% 1:ncol(pcdata))]
    if(length(flip > 0)) pcdata[, flip] <- pcdata[, flip] * -1
  }
  
  plot_args <- list(...)
  plot_args$x = pcdata[,1]
  plot_args$y = pcdata[,2]
  
  if(x$alignment == "principal") {
    xlabel <- paste("PC ", axis1, ": ", round(v[axis1] * 100, 2), "%", sep = "")
    ylabel <- paste("PC ", axis2, ": ", round(v[axis2] * 100, 2), "%", sep = "")
  } else {
    xlabel <- paste("C ", axis1,  sep = "")
    ylabel <- paste("C ", axis2,  sep = "")
  }

  if(is.null(plot_args$xlab)) plot_args$xlab <- xlabel
  if(is.null(plot_args$ylab)) plot_args$ylab <- ylabel
  if(is.null(plot_args$xlim)) plot_args$xlim <- 1.05*range(plot_args$x)
  if(is.null(plot_args$ylim)) plot_args$ylim <- 1.05*range(plot_args$y)
  if(is.null(plot_args$asp)) plot_args$asp <- 1
  
  do.call(plot.default, plot_args)
  
  if(include.axes){
    abline(h = 0, lty=2, ...)
    abline(v = 0, lty=2, ...)
  }
  
  out <- list(points = pcdata,   
              call = match.call())
  
  out$plot_args <- plot_args
  class(out) <- "plot.ordinate"
  options(warn = wrn)
  invisible(out)
  
}


# looCV

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{looCV}}
#' @param ... Other arguments passed onto print.looCV
#' @method print looCV
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.looCV <- function(x, ...){
  cat("Cross-validated scores for", length(x$d$cv), "components,\n")
  cat("based on", NROW(x$scores$cv), "observations.\n\n")
  cat("Cross-validated scores should not be used as data in subsequent analyses.\n")
} 

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{looCV}}
#' @param ... Other arguments passed onto print.looCV
#' @method summary looCV
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
summary.looCV <- function(object, ...){
  x <- object
  print.looCV(x, ...)
  dobs <- x$d$obs
  dcv <- x$d$cv
  pobs <- dobs/sum(dobs)
  cpobs <- cumsum(dobs)/sum(dobs)
  pcv <- dcv/sum(dcv)
  cpcv <- cumsum(dcv)/sum(dcv)
  
  robs <- as.data.frame(rbind(dobs, pobs, cpobs))
  robs <- robs[, 1:min(length(dobs), NCOL(x$scores$obs), NCOL(robs))]
  
  rcv <- as.data.frame(rbind(dcv, pcv, cpcv))
  rcv <- rcv[, 1:min(length(dcv), NCOL(x$scores$cv), NCOL(rcv))]
  
  colnames(robs) <- colnames(rcv) <- colnames(x$x)[1:NCOL(robs)]
  rownames(robs) <- rownames(rcv) <- c("Eigenvalue", 
                   "Proportion of Variance", "Cumulative Proportion")
  
  cat("\nObserved eigenvalues")
  print(robs)
  cat("\nCross-validated eigenvalues")
  print(rcv)
  
}


#' Plot Function for RRPP
#' 
#' @param x An object of class \code{\link{looCV}}
#' @param axis1 A value indicating which component should be 
#' displayed as the X-axis (default = C1)
#' @param axis2 A value indicating which component should be 
#' displayed as the Y-axis (default = C2)
#' @param flip An argument that if not NULL can be used to flip 
#' components in the plot.  
#' The values need to match axis1 or axis2.  For example, if axis1 = 3 
#' and axis2 = 4, flip = 1 will not
#' change either axis; flip = 3 will flip only the horizontal axis; 
#' flip = c(3, 4) will flip both axes.  Axis will only be flipped in first
#' plot.
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' @method plot looCV
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.looCV<- function(x, axis1 = 1, axis2 = 2, 
                      flip = NULL, ...) {
  options(warn = -1)
  if(NCOL(x$scores$obs) == 1) 
    stop("Only one component.  No plotting capability with this function.\n", 
         call. = FALSE)
  opars <- par()
  plot_args <- list(...)
  
  if(axis1 > length(x$d$obs) || axis2 > length(x$d$obs))
    stop("Choice of at least one axis exceeds total axes possible.\n",
         call. = FALSE)
  
  par(mfrow = c(1, 3))
  
  pcdata <- x$scores$obs[, c(axis1, axis2)]
  
  if(!is.null(flip)) {
    if(length(flip) > 2) flip <- flip[1:2]
    flip <- flip[(flip %in% 1:ncol(pcdata))]
    if(length(flip > 0)) pcdata[, flip] <- pcdata[, flip] * -1
  }
  
  plot_args$main <- NULL
  plot_args$x <- pcdata[, 1]
  plot_args$y <- pcdata[, 2]
  plot_args$xlab <- paste("PC", axis1, "for fitted values:", 
                          round(x$d$obs[axis1]/sum(x$d$obs) * 100, 2),
                          "%")
  plot_args$ylab <- paste("PC", axis2, "for fitted values:", 
                          round(x$d$obs[axis2]/sum(x$d$obs) * 100, 2),
                          "%")
  do.call(plot, plot_args)
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  title("Observed PC values")
  
  plot_args$x <- as.matrix(x$scores$cv)[, axis1]
  plot_args$y <- as.matrix(x$scores$cv)[, axis2]
  plot_args$xlab <- paste("PC", axis1, "for fitted values:", 
                          round(x$d$cv[axis1]/sum(x$d$cv) * 100, 2),
                          "%")
  plot_args$ylab <- paste("PC", axis2, "for fitted values:", 
                          round(x$d$cv[axis2]/sum(x$d$cv) * 100, 2),
                          "%")
  do.call(plot, plot_args)
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  title("Cross-validated PC values")
  
  k <- seq(1, min(c(length(x$d$obs), length(x$d$cv))))
  plot_args$x <-  x$d$obs[k]
  plot_args$y <- x$d$cv[k]
  plot_args$xlab <- "Observed eigenvalues"
  plot_args$ylab <- "Cross-validated eigenvalues"
  plot_args$pch <- 19
  plot_args$cex = 1
  plot_args$col = 1
  emax <- max(c(x$d$obs, x$d$cv))
  plot_args$xlim <- c(0, emax)
  plot_args$ylim <- c(0, emax)
  plot_args$asp <- 1
  do.call(plot, plot_args)
  title("Values close to 1:1 imply robust observed scores", 
        cex.main = 0.6)
  abline(0, 1, lty = 3)
 
  par(opars)
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{measurement.error}}
#' @param ... Other arguments passed onto measurement.error
#' @method print measurement.error
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples
#' 
#' 
#' # Measurement error analysis on simulated data of fish shapes
#' 
#' data(fishy)
#' 
#' # Analysis unconcerned with groups 
#' 
#' ME1 <- measurement.error(
#'   Y = "coords",
#'   subjects = "subj",
#'   replicates = "reps",
#'   data = fishy)
#' 
#' anova(ME1)
#' ICCstats(ME1, subjects = "Subjects", with_in = "Systematic ME")
#' plot(ME1)
#' 
#' # Analysis concerned with groups 
#' 
#' ME2 <- measurement.error(
#'   Y = "coords",
#'   subjects = "subj",
#'   replicates = "reps",
#'   groups = "groups",
#'   data = fishy)
#' 
#' anova(ME2)
#' ICCstats(ME2, subjects = "Subjects", 
#'   with_in = "Systematic ME", groups = "groups")
#' P <- plot(ME2)
#' focusMEonSubjects(P, subjects = 18:20, shadow = TRUE)
#' 
print.measurement.error <- function(x, ...){
  
  print(anova(x))
  
  ev <- Re(eigen(x$SSCP.ME.product)$values)
  ev <- ev[zapsmall(ev) > 0]
  
  cat("\n\nNumber of SNR components (eigenvectors):", length(ev), "\n\n")
  cat("Use plot.measurement.error to visualize results\n\n")
  cat("Use function, ICCstats, to get ICC statistics.\n\n")
  
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{measurement.error}}
#' @param ... Other arguments passed onto measurement.error
#' @method summary measurement.error
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.measurement.error <- function(object, ...){
  
  print.measurement.error(object)
  
}

#' Plot Function for RRPP
#' 
#' This function produces multivariate signal-to-noise ratio plots for  
#' \code{\link{measurement.error}} objects.  See the function, 
#' \code{\link{plot.interSubVar}} for plotting the inter-subject variability
#' from a \code{\link{measurement.error}} object, after applying the function, 
#' \code{\link{interSubVar}}.
#'
#' @param x Object from \code{\link{measurement.error}}
#' @param separate.by.groups A logical value for whether to make separate plots
#' for each group, if different groups are available.  If FALSE, groups
#' are still represented by different symbols in the plot, unless overridden
#' by plot arguments.
#' @param add.connectors A logical value for whether to add connectors, like
#' vectors, between replicate observations of the same subjects.  
#' @param titles An optional vector or list for augmenting the titles of plots produced.
#' The length of the vector or list should match the number of plots produced by 
#' other arguments.
#' @param add.labels A logical value for whether to label subjects.  
#' Labels are either subject name (if available) or number of occurrence
#' in the data set.
#' @param use.std.vectors A logical value for whether to use vectors obtained from
#' a standardized matrix, which are orthogonal.  This is not strictly necessary.
#' @param add.legend A logical value for whether to add a legend to plots.  If
#' separate.by.groups is TRUE, adding a legend to plots will be slightly redundant. If
#' certain parameters are augmented by user (point characters, colors), add.legend
#' will be made to be FALSE to prevent misinterpretation of intended plotting scheme.
#' @param ... Other arguments passed onto plot
#' @method plot measurement.error
#' @export
#' @author Michael Collyer
#' @keywords utilities
plot.measurement.error <- function(x, 
                                   separate.by.groups = TRUE,
                                   add.connectors = TRUE,
                                   add.labels = FALSE,
                                   use.std.vectors = FALSE,
                                   titles = NULL,
                                   add.legend = TRUE, ...){
  
  subj <- x$LM$data$subjects
  reps <- x$LM$data$replicates
  groups <- if(!is.null(x$LM$data$groups))
    x$LM$data$groups else NULL
  
  plot_args <- list(...)
  if(!is.null(titles)) titles <- as.list(titles)
  
  if(is.null(plot_args$pch)) {
    plot_args$pch <- 20 
    if(length(x$LM$data) > 3) {
      plot_args$pch <- plot_args$pch + as.numeric(groups)
      ldf.pch <- unique(data.frame(a = groups, 
                                  b = plot_args$pch))
    } else {
      plot_args$pch <- plot_args$pch + 1
      ldf.pch <- data.frame(a = 1, b = plot_args$pch)
    }
  } else {
    if(length(x$LM$data) > 3) {
      ldf.pch <- unique(data.frame(a = groups, 
                                  b = plot_args$pch))
    } else {
      ldf.pch <-unique(data.frame(a = 1, b = plot_args$pch))
    }
  }
  
  if(is.null(plot_args$bg)) plot_args$bg <- as.numeric(reps)
  ldf.bg <- unique(data.frame(a = reps, b = plot_args$bg))
  
  if(is.null(plot_args$col)) plot_args$col <- as.numeric(reps)
  ldf.col <- unique(data.frame(a = reps, b = plot_args$col))
  
  
  if(!is.null(x$SSCP.ME.product)) {
    S <- if(use.std.vectors) svd(x$SSCP.ME.product.std) else
      svd(x$SSCP.ME.product)
  } else S <- list(d = 1, u = 1, v = 1)
  
  Y <- center(x$LM$Y)
  if(length(x$LM$data) > 3){
    gp <- x$LM$data$groups
    Y <- resid(lm(Y ~ gp))
  }
    
  d <- S$d
  dx <- d[1]/sum(d)
  dy <- if(length(d) > 1) d[2]/sum(d) else 0
  
  pts <- as.matrix(Y) %*% S[[3]] %*% 
    diag(sqrt(S[[1]]))[, 1:(min(length(d), 2))]
  
  if(NCOL(pts) == 1) pts <- cbind(pts, 0)
  
  plot_args$asp <- 1
  if(is.null(plot_args$xlab))
    plot_args$xlab <- paste("EV 1:", round(dx *100, 2), "%")
  if(is.null(plot_args$ylab))
    plot_args$ylab <- paste("EV 2:", round(dy *100, 2), "%")
  if(is.null(plot_args$main))
    plot_args$main <- "Systematic ME / Random ME"
  if(!is.null(titles))
    plot_args$main <- titles[[1]]
  if(is.null(plot_args$cex.main))
    plot_args$cex.main <- 1
  
  add.me.connectors <- function(x, y, subj, d, maxy) {
    
    dindx <- which(zapsmall(d) > 0)
    if(length(dindx) == 1) {
      y <- y + 0.1 * maxy
      subjjig <- seq(0, 0.9 * maxy, (0.9 * maxy) / (nlevels(subj) - 1))
      for(i in 1:length(subjjig)) {
        y[which(subj == levels(subj)[i])] <- 
          y[which(subj == levels(subj)[i])] + subjjig[i]
      }
    }
    subjlv <- levels(subj)
    for(i in 1:length(subjlv)){
      xx <- x[subj == subjlv[i]]
      yy <- y[subj == subjlv[i]]
      xmean <- mean(xx)
      ymean <- mean(yy)
      for(j in 1:length(xx)) {
        arrows(xmean, ymean, xx[j], yy[j], length = 0, 
               angle = 0, lwd = 0.5, ...)
        
      }
    }
  }
  
  add.me.labels <- function(x, y, subj, d, maxy) {
    dindx <- which(zapsmall(d) > 0)
    if(length(dindx) == 1) {
      y <- y + 0.1 * maxy
      subjjig <- seq(0, 0.9 * maxy, (0.9 * maxy) / (nlevels(subj) - 1))
      for(i in 1:length(subjjig)) {
        y[which(subj == levels(subj)[i])] <- 
          y[which(subj == levels(subj)[i])] + subjjig[i]
      }
    }
    subjlv <- levels(subj)
    for(i in 1:length(subjlv)){
      xx <- x[subj == subjlv[i]]
      yy <- y[subj == subjlv[i]]
      xmean <- mean(xx)
      ymean <- mean(yy)
      for(j in 1:length(xx)) {
        text(xmean, ymean, subjlv[i], pos = 3, cex = 0.4, offset = 0.1)
      }
    }
  }
  
  plot_args$x <- as.matrix(pts)[,1]
  plot_args$y <- if(NCOL(as.matrix(pts)) == 1) 
    rep(0, length(pts)) else
      as.matrix(pts)[,2]
  if(is.null(plot_args$ylim)) {
    ymax <- max(c(plot_args$y, max(plot_args$x)))
    ymin <- min(plot_args$y)
    plot_args$ylim <- c(ymin, ymax)
  }
  
  if(separate.by.groups && is.null(x$LM$data$groups))
    separate.by.groups <- FALSE
  
  if(separate.by.groups) {
    nplots <- nlevels(x$LM$data$groups)
    point_args <- plot_args
    plot_args$type <- "n"
    plot_args$main <- NULL
    gps <- levels(x$LM$data$groups)
    if(!is.null(titles)) {
      if(length(titles) != nplots) {
        if(length(titles) > nplots)
          titles <- titles[1:nplots] else
            titles <- rep(titles, nplots)[1:nplots]
      }
    }
    
    for(i in 1:nplots) {
      do.call(plot, plot_args)
      if(is.null(titles)) {
        title(paste("Systematic ME / Random ME for", 
                    "Group", gps[[i]], sep = " "),
              cex.main = plot_args$cex.main)
      } else {
        title(titles[[i]], cex.main = plot_args$cex.main)
      }
      
      point_args.b <- point_args
      point_args.b$pch[which(x$LM$data$groups != gps[[i]])] <- NA
      
      do.call(points, point_args.b)
      
      if(add.connectors) add.me.connectors(
        point_args.b$x[which(!is.na(point_args.b$pch))],
        point_args.b$y[which(!is.na(point_args.b$pch))],
        factor(subj[which(!is.na(point_args.b$pch))]), 
        d, point_args.b$ylim[2])
      if(add.labels) add.me.labels(
        point_args.b$x[which(!is.na(point_args.b$pch))],
        point_args.b$y[which(!is.na(point_args.b$pch))],
        factor(subj[which(!is.na(point_args.b$pch))]), 
        d, point_args.b$ylim[2])
      
      if(add.legend) {
        
        legend("topleft", 
               legend = ldf.col$a,
               pch = unique(na.omit(point_args.b$pch)),
               col = ldf.col$b,
               pt.bg = ldf.bg$b,
               bg = "white")
      }
    }
    
  } else {
    do.call(plot, plot_args)
    
    if(add.connectors)
      add.me.connectors(plot_args$x, plot_args$y, subj, d, plot_args$ylim[2])
    if(add.labels)
      add.me.labels(plot_args$x, plot_args$y, subj, d, plot_args$ylim[2])
    
    if(add.legend) {
      legend("topleft", 
             legend = ldf.col$a,
             pch = min(ldf.pch$b),
             col = ldf.col$b,
             pt.bg = ldf.bg$b,
             bg = "white")

    }
  }
  
  plot_args$data <- x$LM$data
  class(plot_args) <- "plot.measurement.error"
  invisible(plot_args)
  
}

#' Plot Function for RRPP
#' 
#' This function produces a heat map for inter-subject variability, based on 
#' results from a \code{\link{measurement.error}} object.  The function, 
#' \code{\link{interSubVar}}, must first be used on the \code{\link{measurement.error}} 
#' object to obtain variability statistics.  This function use the \code{\link{image}}
#' function to produce plots.  It does little to manipulate such plots, but any
#' argument for \code{\link{image}} can be manipulated here, as well as the graphical
#' parameters that can be adjusted within \code{\link{image}}.
#'
#' @param x Object from \code{\link{interSubVar}}
#' @param ... Arguments passed onto \code{\link{image}} and \code{\link{plot}}.
#' @method plot interSubVar
#' @export
#' @author Michael Collyer
#' @keywords utilities
plot.interSubVar <- function(x, ...){
  subj <- x$subject.order
  n <- length(subj)
  names(subj) <- 1:n
  cat("If not apparent in the plot, the order of subjects on axes from 1 to",
      n, "is:\n")
  print(subj)
  cat("\n\n")
  image.args <- list(...)
  d <- as.matrix(x$var.map)
  image.args$x <- image.args$y <- 1:n
  image.args$z <- d
  image.args$xlab <- image.args$ylab <- "Subjects"
  do.call(image, image.args)
  out <- image.args
  out$subj <- subj
  class(out) <- "plot.interSubVar"
  invisible(out)
}

#' Plot Function for RRPP
#' 
#' Reduces a plot.measurement.error to a single research subject.  This can be 
#' for cases when many overlapping subjects in a plot obscure interpretation 
#' for specific subjects.
#'
#' @param x Plot from \code{\link{plot.measurement.error}}
#' @param subjects The specific subject to plot
#' @param shadow A logical value for whether to show other subject values as
#' shadows of their locations.
#' @param ... Other arguments passed onto plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
focusMEonSubjects <- function(x, subjects = NULL, 
                              shadow = TRUE, ...){
  if(is.null(subjects))
    stop("Please specify at least 
         one subject to plot.\n", call. = FALSE)
  if(length(unique(subjects)) > length(x$x)) 
    stop("Please specify subjects within the vector origianlly used.\n",
         call. = FALSE)
  
  keep <- which(x$data$subj %in% unique(subjects))
  if(length(keep) == 0)
  stop("Subjects names not found within subject factor.\n", call. = FALSE)

  no.trm <- which(names(x) == "data")
  plot_args <- x[-no.trm]
  plot_args$type <- "n"
  do.call(plot, plot_args)
  if(shadow){
    plot_args$type <- "p"
    plot_args$cex <- if(!is.null(plot_args$cex)) 
      min(c(plot_args$cex / 2, 0.5)) else 0.5
    if(!is.null(plot_args$col)) plot_args$col <- alpha(plot_args$col, 0.1)
    if(!is.null(plot_args$bg)) plot_args$bg <- alpha(plot_args$bg, 0.1)
    do.call(points, plot_args)
  }

  point_args <- list(x = x$x[keep], y = x$y[keep])
  if(!is.null(x$pch)) point_args$pch <- if(length(x$pch) > 1)
    x$pch[keep] else x$pch
  if(!is.null(x$col)) point_args$col <- if(length(x$col) > 1)
    x$col[keep] else x$col
  if(!is.null(x$bg)) point_args$bg <- if(length(x$bg) > 1)
    x$bg[keep] else x$bg
  if(!is.null(x$lty)) point_args$lty <- if(length(x$lty) > 1)
    x$lty[keep] else x$lty
  if(!is.null(x$lwd)) point_args$lwd <- if(length(x$lwd) > 1)
    x$lwd[keep] else x$lwd
  if(!is.null(x$cex)) point_args$cex <- x$cex
  do.call(points, point_args)

  for(i in 1:length(subjects)) {
    keep <- which(x$data$subj %in% subjects[i])
    xx <- x$x[keep]
    yy <- x$y[keep]
    xmean <- mean(xx)
    ymean <- mean(yy)
    for(j in 1:length(xx)) {
      arrows(xmean, ymean, xx[j], yy[j], length = 0, 
             angle = 0, lwd = 0.5, ...)
      
    }
  }
  
}

#' Utility Function for RRPP
#' 
#' A function mostly for internal processing but can be used to extract
#' ANOVA statistics for other uses, such as plotting histograms
#'
#' @param fit Object from \code{\link{lm.rrpp}}.
#' @param stat The ANOVA statistic to extract.  Returns every RRPP permutation of 
#' the statistic.  If "all", a list of each statistic is returned.  
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples
#' 
#' data(Pupfish)
#' fit <- lm.rrpp(coords ~ log(CS) + Sex*Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE, iter = 999) 
#' anova(fit)
#' Fstats <- getANOVAStats(fit, stat = "F")
#' par(mfrow = c(2, 2))
#' hist(Fstats$Fs[1,], breaks = 50, main = "log(CS)", xlab = "F")
#' abline(v = Fstats$Fs[1, 1])
#' hist(Fstats$Fs[2,], breaks = 50, main = "Sex", xlab = "F")
#' abline(v = Fstats$Fs[2, 1])
#' hist(Fstats$Fs[3,], breaks = 50, main = "Pop", xlab = "F")
#' abline(v = Fstats$Fs[3, 1])
#' hist(Fstats$Fs[4,], breaks = 50, main = "Sex:Pop", xlab = "F")
#' abline(v = Fstats$Fs[4, 1])
#' 
getANOVAStats <- function(fit, stat = c("SS", "MS", "Rsq", "F", "cohenf", "all")){
  ANOVA <- fit$ANOVA
  verbose <- fit$verbose
  SS <- ANOVA$SS
  MS <- ANOVA$MS
  Fs <- ANOVA$Fs
  cohenf <- ANOVA$cohenf
  Rsq <- ANOVA$Rsq
  Df <- ANOVA$df
  RSS <- ANOVA$RSS
  TSS <- ANOVA$TSS
  SS.type <- ANOVA$SS.type
  k <- length(fit$LM$term.labels)
  
  if(is.null(SS)){
    k <- 0
    SS <- ANOVA$RSS.model
    MS <- SS/Df
    Rsq <- rep(1, length(SS))
    names(Rsq) <- names(SS)
    Fs <- cohenf <- rep(0, length(SS))
    names(Fs) <- names(cohenf) <- names(SS)
    verbose <- TRUE
  }
  
  if(!verbose) {
    MS <- SS / Df[1:k]
    MSE <- RSS / Df[k + 1]
    Fs <- MS / MSE
    Rsq <- SS/TSS
    cohenf <- Rsq / (1 - Rsq)
    if(SS.type != "III"){
      etas <- Rsq
      if(k == 1) unexp <- 1 - etas else unexp <- 1 - apply(etas, 2, cumsum)
      cohenf <- etas/unexp
    }
  }
  out <- list(SS = SS, MS = MS, Rsq = Rsq, Fs = Fs, cohenf = cohenf)
  stat <- match.arg(stat)
  keep <- which(c("SS", "MS", "Rsq", "F", "cohenf", "all") %in% stat)
  if(keep == 6) keep <- 1:5
  out[keep]
}

#' Utility Function for RRPP
#' 
#' A function mostly for internal processing but can be used to extract
#' RRPP permutation information for other reasons.
#'
#' @param fit Object from \code{\link{lm.rrpp}}
#' @param attribute The various attributes that are used to generate RRPP
#'  permutation schedules.  If there are n observations, each iteration has
#'  some randomization of 1:n, restricted by the arguments that match attributes
#'  provided by this function.
#' @export
#' @author Michael Collyer
#' @keywords utilities

getPermInfo <- function(fit, attribute = c("perms", "perm.method",
                                           "block", "perm.seed",
                                           "perm.schedule", "all")){
  
  PermInfo <- fit$PermInfo
  perms <- PermInfo$perms
  perm.method <- PermInfo$perm.method
  block <- PermInfo$block
  perm.seed <- PermInfo$perm.seed
  perm.schedule <- PermInfo$perm.schedule
  n <- fit$LM$n
  
  if(is.null(perm.schedule)) {
    if(perms == 1) perms <- 2
    perm.schedule <- perm.index(n, perms - 1, block = block, 
                                seed = perm.seed)
  }
  
  out <- list(perms = perms, perm.method = perm.method, block = block,
              perm.seed = perm.seed, perm.schedule = perm.schedule)
  
  attribute <- match.arg(attribute)
  keep <- which(c("perms", "perm.method",
                  "block", "perm.seed",
                  "perm.schedule", "all") %in% attribute)
  if(keep == 6) keep <- 1:5
  out[keep]
}

#' Utility Function for RRPP
#' 
#' A function mostly for internal processing but can be used to extract
#' The terms for each reduced and full model used in an \code{\link{lm.rrpp}} fit.
#'
#' @param fit Object from \code{\link{lm.rrpp}}
#' @export
#' @author Michael Collyer
#' @keywords utilities

getTerms <- function(fit){
  Model.Terms <- .getTerms(fit)
  Model.Terms <- lapply(Model.Terms, function(x){
    lapply(x, function(y){
      attr(y, "dataClasses") <- NULL
      y
    })
  })
  
  names(Model.Terms) <- c("reduced", "full")
  Model.Terms
}

#' Utility Function for RRPP
#' 
#' A function mostly for internal processing but can be used to obtain terms,
#' design matrices, or QR decompositions used for each reduced or full
#' model that is fitted in an \code{\link{lm.rrpp}} fit.
#'
#' @param fit Object from \code{\link{lm.rrpp}}
#' @param attribute The various attributes that are used to extract RRPP
#'  permutation schedules.  If there are n observations, each iteration has
#'  some randomization of 1:n, restricted by the arguments that match attributes
#'  provided by this function.
#' @export
#' @author Michael Collyer
#' @keywords utilities

getModels <- function(fit, attribute = c("terms", "X", "qr", "all")) {
  
  attribute <- match.arg(attribute)
  create <- is.null(fit$Models)
  if(!create) Models <- fit$Models else Models <- NULL
  
  Pcov <- if(!is.null(fit$LM$Pcov)) fit$LM$Pcov else
    if(!is.null(fit$LM$Cov)) Cov.proj(fit$LM$Cov) else NULL
  w <- if(!is.null(fit$LM$weights)) sqrt(fit$LM$weights) else NULL

  Model.Terms <- getTerms(fit)
  STerm <- if(!is.null(fit$subjects.var))
    which(fit$LM$term.labels == fit$subjects.var) else NULL
  
  if(create){
    Xs <- suppressWarnings(getXs(Terms = fit$LM$Terms, 
                                  Y = fit$LM$Y, 
                                  SS.type = fit$ANOVA$SS.type,
                                  model = fit$LM$data,
                                 subjects.term = STerm)) 
  } else {
    Xs <- lapply(Models, function(x){
      lapply(x, function(y) y$X )
    })
  }
  names(Xs) <- c("reduced", "full")
  k <- length(Xs[[1]])
  
  if(attribute == "all" || attribute == "qr") {
    
    if(create) {
      Models <-lapply(1:2, function(j){
        res <- lapply(1:max(1, k), function(jj){
          X <- Xs[[j]][[jj]]
          TX <- if(!is.null(Pcov)) Pcov %*% X else if(!is.null(w)) X*w else X
          qr <- qr(TX)
          out <- list(X = X, qr = qr)
        })
        
        res
      })
      
      for(i in 1:2){
        for(j in 1:max(1, k)){
          Models[[i]][[j]]$terms <- Model.Terms[[i]][[j]]
        }
      }
      names(Models) <- c("reduced", "full")
      names(Models$reduced) <- names(Xs[[1]])
      names(Models$full) <- names(Xs[[2]])
    } else Models <- fit$Models
    
  } else Models <- NULL
  
  if(!is.null(Models)){
    Models <- lapply(Models, function(x){
      lapply(x, function(y){
        attr(y$terms, "dataClasses") <- NULL
        y
      })
    })
  }
  
  if(attribute == "all")
    out <- Models
  if(attribute == "terms")
    out <- Model.Terms
  if(attribute == "X")
    out <- Xs
  if(attribute == "qr"){
    QR <- lapply(Models, function(x){
      lapply(x, function(y) y$qr)
    })
    out <- QR
  }
  out
}

#' Utility Function for RRPP
#' 
#' A function mostly for internal processing but can be used to extract either
#' the model covariance matrix (Cov) or the projection matrix for transformations made 
#' from the covariance matrix (Pcov), which is basically the square-root of 
#' the covariance matrix.  This matrix is the model covariance used for estimation,
#' not the residual covariance matrix (see \code{\link{getResCov}}). 
#' There are also options for S3 or S4 format versions, or
#' a forcing of symmetry for Pcov.
#'
#' @param fit Object from \code{\link{lm.rrpp}}
#' @param type Whether the Cov or Pcov matrix is returned
#' @param format Whether an S3 or S4 format is returned
#' @param forceSym Logical value for whether a symmetric matrix should 
#' be returned for Pcov, even if Pcov was triangular as a solution.
#' @export
#' @author Michael Collyer
#' @keywords utilities

getModelCov <- function(fit, type = c("Cov", "Pcov"), 
                        format = c("S3", "S4"),
                        forceSym = TRUE) {
  if(is.null(fit$LM$Cov))
    stop("There was no model covariance matrix used in this lm.rrpp fit.\n",
         call. = FALSE)
  if(is.null(fit$LM$Pcov)) fit$LM$Pcov <- Cov.proj(fit$LM$Cov)
  
  type <- match.arg(type)
  format <- match.arg(format)
  
  res <- if(type == "Pcov") fit$LM$Pcov else fit$LM$Cov
  res <- if(format == "S4") Matrix(res, sparse = TRUE) else
    as.matrix(res)
  isSym <- (res[2, 1] == res[1, 2])
  if(!isSym && forceSym && type == "Pcov"){
    res <- Cov.proj(fit$LM$Cov, symmetric = TRUE)
    res <- if(format == "S4") Matrix(res, sparse = TRUE) else
      as.matrix(res)
  }
    
  
  res
}

#' Utility Function for RRPP
#' 
#' A function mostly for internal processing but can be used to extract the residual
#' covariance matrix.  This matrix is the residual covariance matrix,
#' not the model covariance matrix used for estimation (see \code{\link{getModelCov}}). 
#' Options for averaging over degrees of freedom or number of
#' observations, plus standardization, are also available.
#'
#' @param fit Object from \code{\link{lm.rrpp}}
#' @param useDf Logical value for whether the degrees of freedom from the linear model
#' fit should be used (if TRUE), as opposed to the number of observations (if FALSE).
#' @param standardize Logical value for whether residuals should be standarized.  If TRUE,
#' a correlation matrix is produced.
#' @export
#' @author Michael Collyer
#' @keywords utilities

getResCov <- function(fit, useDf = TRUE, standardize = FALSE) {
  gls <- fit$LM$gls
  if(gls && is.null(fit$LM$Pcov) && !is.null(fit$LM$Cov)) 
    fit$LM$Pcov <- Cov.proj(fit$LM$Cov)
  if(gls && !is.null(fit$LM$Pcov)) R <- fit$LM$Pcov %*% fit$LM$gls.residuals
  if(gls && is.null(fit$LM$Pcov)) R <- fit$LM$gls.residuals * sqrt(fit$LM$weights)
  if(!gls) R <-  fit$LM$residuals
  S <- as.matrix(crossprod(R))
  if(useDf){
    df <- fit$ANOVA$df
    if(length(df) > 1) df <- df[(length(df) - 1)] 
  } else df <- fit$LM$n
  
  S <- S / df
  
  if(standardize && length(S) > 1) {
    w <- diag(1 / sqrt(diag(S))) 
    S <- w %*% S %*% w
  } else if(standardize && length(S) == 1) S <- 1
    
    
  S
  
}

#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{ICCstats}})
#' @param ... other arguments passed to print/summary
#' @method print ICCstats
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.ICCstats <-function(x, ...){
  
  cat("\nICC stats based on dispersion of values\n\n")
  print(x$ICC_disp)
  
  if(!is.null(x$ICC_mult)) {
    cat("\nICC stats by Eignevector of multivariate ICC matrix\n\n")
    print(x$ICC_mult)
  }
    
}

#' Print/Summary Function for RRPP
#'
#' @param object print/summary object (from \code{\link{ICCstats}})
#' @param ... other arguments passed to print/summary
#' @method summary ICCstats
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.ICCstats <-function(object, ...){
  print.ICCstats(object, ...)
}