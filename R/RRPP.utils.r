## lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{lm.rrpp}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.lm.rrpp <- function(x, ...){
  cat("\nLinear Model fit with lm.rrpp\n")
  LM <- x$LM
  PI <- x$PermInfo
  AN <- x$ANOVA
  if(!is.null(x$LM$dist.coefficients)) dv <- "(dimensions of data after PCoA of distance matrix)" else
    dv <- " "
  cat(paste("\nNumber of observations:", LM$n))
  cat(paste("\nNumber of dependent variables:", LM$p, dv))
  if(!is.null(AN$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", AN$SS.type))
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
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.lm.rrpp <- function(object, formula = TRUE, ...){
  x <- object
  LM <- x$LM
  PI <- x$PermInfo
  AN <- x$ANOVA
  perms <- PI$perms
  dv <- LM$dist.coefficients
  n <- LM$n
  p <- LM$p
  SS <- AN$SS
  SS.type <- AN$SS.type
  k <- length(LM$term.labels)
  if(LM$gls) P <- LM$Pcov
  
  if(k > 0) {
    SSE <- unlist(SS[k+1,])
    SST <- unlist(SS[k+2,])
    SSM <- SST - SSE
    df <- AN$df
    dfe <- df[k+1]
    dfM <- sum(df[1:k])
    Fs <- (SSM/dfM)/(SSE/dfe)
    P <- pval(log(Fs))
    Z <- effect.size(Fs)
    Rsq <- SSM[1]/SST[1]
    SSM.obs <- SSM[1]
    Fs.obs <- Fs[1]
  } else {
      df <- 0
      dfe <- AN$df
      dfM <- 1
      SST <- SSE <- unlist(SS)
      SSM <- Fs <- Z <- Rsq <- SSM.obs <- Fs.obs <- P <- ""
    }

  tab <- data.frame(dfM = dfM, dfe = dfe,
                    SSM = SSM.obs, RSS = SSE[1], Rsq = Rsq,
                    F = Fs.obs, Z = Z, P = P)
  dimnames(tab)[[2]] <- c("Df",
                          "Residual Df",
                          "SS",
                          "Residual SS",
                          "Rsq",
                          "F",
                          "Z (from F)",
                          "Pr(>F)")
  if(formula) dimnames(tab)[[1]] <- deparse(x$call[[2]]) else
    dimnames(tab)[[1]] <- deparse(substitute(object))

  rfit <-refit(x)
  RR <- rfit$wResiduals.reduced
  RF <- rfit$wResiduals.full
  SSCP <- lapply(1:length(RF), function(j) crossprod(RR[[j]] - RF[[j]]))
  names(SSCP) <- LM$term.labels
  SSCP <- c(SSCP, list(Residuals = as.matrix(crossprod(RF[[length(RF)]]))))
  out <- list(table = tab, SSCP = SSCP, n = n, p = p, k = k, 
              perms = perms, dv = dv, SS = SS, SS.type = SS.type)
  class(out) <- "summary.lm.rrpp"
  out
}

#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{summary.lm.rrpp}})
#' @param ... other arguments passed to print/summary
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
  if(!is.null(x$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
  cat(paste("\nNumber of permutations:", x$perms))
  cat("\n\nFull Model Analysis of Variance\n\n")
  print(x$table)
  invisible(x)
}
## coef.lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{coef.lm.rrpp}}
#' @param ... Other arguments passed onto coef.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.coef.lm.rrpp <- function(x, ...){
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of observations:", x$n))
  cat(paste("\nNumber of dependent variables:", x$p))
  cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
  cat(paste("\nNumber of permutations:", x$nperms))
  if(!x$test) {
    cat("\n\nObserved coefficients\n\n")
    print(x$coef.obs)
  }
  if(x$test){
    rrpp.type <- x$RRPP
    cat("\n\nStatistics (distances) of coefficients with ")
    cat(x$confidence*100, "percent confidence intervals,") 
    cat("\neffect sizes, and probabilities of exceeding observed values based on\n") 
    cat(x$nperms, "random permutations using", rrpp.type, "\n\n")
    print(x$stat.tab)
    cat("\n\n")
  }
  invisible(x)
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{coef.lm.rrpp}}
#' @param ... Other arguments passed onto coef.lm.rrpp
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
#' @param PC Logical argumnet for whether to use predicted values 
#' rotated to their PCs
#' @param ... Other arguments passed onto predict.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.predict.lm.rrpp <- function(x, PC = FALSE, ...){
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of predictions:", NROW(x$mean)))
  cat(paste("\nConfidence level:", x$confidence*100, "%"))
      cat(paste("\nNumber of bootstrap permutations:", length(x$random.predicted)))
      if(PC)  cat(paste("\nPredicted values are rotated to their PCs"))
      cat("\n\nPredicted values:\n\n")
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
    if(NROW(tab) > 1) cat(paste("Effect sizes (Z) based on", effect.type, "distributions\n\n"))
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
#' @param predictor An optional vector if "regression" plot type is chosen, 
#' and is a variable likely used in \code{\link{lm.rrpp}}.
#' This vector is a vector of covariate values equal to the number of observations.
#' @param reg.type If "regression" is chosen for plot type, this argument
#' indicates whether a common regression component (CRC) plot, prediction line 
#' (Predline) plot, or regression score (RegScore) plotting is performed.  
#' For explanation of CRC, see Mitteroecker et al. (2004).  For explanation of prediction line,
#' see Adams and Nistri (2010).  For explanation of regression score, see 
#' Drake and Klingenberg (2008).
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#' Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#' transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#' in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
plot.lm.rrpp <- function(x, type = c("diagnostics", "regression",
                                      "PC"), predictor = NULL,
                          reg.type = c("CRC", "PredLine", "RegScore"), ...){
  r <- as.matrix(x$LM$wResiduals)
  f <- as.matrix(x$LM$wFitted)
  if(!is.null(x$Pcov)) {
    r <- as.matrix(x$LM$gls.residuals)
    f <- as.matrix(x$LM$gls.fitted)
  }
  type <- match.arg(type)
  if(is.na(match(type, c("diagnostics", "regression", "PC")))) 
    type <- "diagnostics"
  CRC <- PL <- Reg.proj <- NULL
  if(type == "diagnostics") {
    pca.r <- prcomp(r)
    var.r <- round(pca.r$sdev^2/sum(pca.r$sdev^2)*100,2)
    plot(pca.r$x, pch=19, asp =1,
         xlab = paste("PC 1", var.r[1],"%"),
         ylab = paste("PC 2", var.r[2],"%"),
         main = "PCA Residuals")
    pca.f <- prcomp(f)
    var.f <- round(pca.f$sdev^2/sum(pca.f$sdev^2)*100,2)
    dr <- scale(sqrt(diag(tcrossprod(center(r)))))
    plot.QQ(r)
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
    plot.het(r,f)
    p <- ncol(r)
  }
  if(type == "regression"){
    reg.type <- match.arg(reg.type)
    if(is.na(match(reg.type, c("CRC", "PredLine", "RegScore")))) 
      if(is.null(predictor))
        stop("This plot type is not available without a predictor.")
    n <- NROW(r); p <- NCOL(r)
    if(!is.vector(predictor)) stop("Predictor must be a vector")
    if(length(predictor) != n) 
      stop("Observations in predictor must equal observations if procD.lm fit")
    X <- x$LM$X * sqrt(x$LM$weights)
    if(!is.null(x$LM$Pcov)) B <- x$LM$gls.coefficients else B <- x$LM$coefficients
    xc <- predictor
    pred.match <- sapply(1:NCOL(X), function(j){
      any(is.na(match(xc, X[,j])))
    })
    if(all(pred.match)) {
      b <- lm(f ~ xc)$coefficients
      if(is.matrix(b)) b <- b[2,] else b <- b[2]
    } else {
      Xcrc <- as.matrix(X)
      Xcrc[,!pred.match] <- 0
      f <- Xcrc %*% B
      r <- x$LM$Y - f
      b <- lm(f ~ xc)$coefficients
      if(is.matrix(b)) b <- b[2,] else b <- b[2]
    }
    a <- crossprod(r, xc)/sum(xc^2)
    a <- a/sqrt(sum(a^2))
    CRC <- r%*%a  
    resid <- r%*%(diag(p) - matrix(crossprod(a),p,p))
    RSC <- prcomp(resid)$x
    Reg.proj <- x$LM$Y%*%b%*%sqrt(solve(crossprod(b)))
    PL <- prcomp(f)$x[,1]
    if(reg.type == "CRC"){
      par(mfcol = c(1,2))
      par(mar = c(4,4,1,1))
      plot(predictor, CRC,  ...)
      plot(CRC, RSC[,1], asp=1, xlab = "CRC", ylab = "RSC 1", ...)
      par(mar = c(5,4,4,2) + 0.1)
      par(mfcol=c(1,1))
    } else if(reg.type == "RegScore") {
      plot(predictor, Reg.proj, 
           ylab = "Regression Score", ...)
    } else {
      plot(predictor, PL, 
           ylab = "PC 1 for fitted values", ...)
    }
  }
  if(type == "PC"){
    pca <- prcomp(f)
    eigs <- pca$rotation
    P <- x$LM$Y%*%eigs
    v <- pca$sdev^2
    ev <- round(v[1:2]/sum(v)*100, 2)
    plot(P, asp=1,
         xlab = paste("PC 1 for fitted values: ",ev[1],"%", sep = ""),
         ylab = paste("PC 2 for fitted values: ",ev[2],"%", sep = ""),
                      ...)
  }
  out <- list(CRC = CRC, PredLine = PL, RegScore = Reg.proj)
  invisible(out)
}

plot.het <- function(r,f){
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

plot.QQ <- function(r){
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
#' @param x plot object (from \code{\link{lm.rrpp}})
#' @param PC A logical argument for whether the data space should be rotated to its 
#' principal components
#' @param ellipse A logical argument to change error bars to ellipses in multivariate plots.  
#' It has no function for univariate plots.
#' @param label A logical argument for whether points should be labeled 
#' (in multivariate plots).
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.predict.lm.rrpp <- function(x, PC = FALSE, ellipse = FALSE,
                                 label = TRUE, ...){
  if(is.matrix(x$mean)) m <- rbind(x$mean, x$lcl, x$ucl) else 
    m <- c(x$mean, x$lcl, x$ucl)
  if(is.matrix(x$pc.mean)) mpc <- rbind(x$pc.mean, x$pc.lcl, x$pc.ucl) else 
    mpc <- c(x$pc.mean, x$pc.lcl, x$pc.ucl)
  conf <- x$confidence
  dots <- list(...)
  if(NCOL(m) == 1) {
    k <- length(x$mean)
    xx <- seq(1:k)
    xf <- as.factor(rownames(x$mean))
    plot(xf, x$mean, lty = "blank",
         xlab = "Predicted values", 
         ylab = x$data.name,
         ylim = c(min(m), max(m)),
         main = paste("Predicted Values, +/-", conf*100, "percent confidence levels",
         sep = " "), pch=19,
         cex.main = 0.7, ...
         )
    arrows(xx, x$mean, xx, x$lcl, lty = 1, angle = 90, length = 0.10, ...)
    arrows(xx, x$mean, xx, x$ucl, lty = 1, angle = 90, length = 0.10, ...)
    if(length(dots) == 0) points(xx, x$mean, pch = 19, cex = 0.7)
    else points(xx, x$mean, ...)
  }
  if(NCOL(m) > 1) {
    type <- "multi"
    k <- NROW(x$mean)
    if(PC){
      pca <- x$pca
      d <-length(which(zapsmall(pca$sdev) > 0))
      if(d == 1) {
        mr <- mpc
        ve <- 100
        ylm <- c(min(mr), max(mr))
        ylb <- paste("PC1: ", ve, "%", sep = "")
        type <- "uni"
      } else {
        mr <- mpc[,1:2]
        v <- pca$sdev^2
        ve <- round(v/sum(v)*100, 2)[1:2]
        xlm <- c(min(mr[,1]), max(mr[,1]))
        ylm <- c(min(mr[,2]), max(mr[,2]))
        xlb <- paste("PC1: ", ve[1], "%", sep = "")
        ylb <- paste("PC2: ", ve[2], "%", sep = "")
        }
    } else {
      mr <- m[,1:2]
      cn <- colnames(x$mean)
      if(is.null(cn)) cn <- paste("V", 1:2, sep=".")
      xlb <- cn[1]
      ylb <- cn[2]
      xlm <- c(min(mr[,1]), max(mr[,1]))
      ylm <- c(min(mr[,2]), max(mr[,2]))
    }
    mt <- if(PC) "Among-prediction PC rotation" else 
      "Plot of first two variables"
    mt <- paste(mt, "; ", conf*100, "% confidence limits", sep = "")
    
    if(PC && type == "uni"){
      xx <- seq(1:k)
      xf <- as.factor(rownames(x$mean))
      plot(xf, mr[1:k], lty = "blank",
           xlab = "Predicted values", 
           ylab = ylb,
           ylim = c(min(mr), max(mr)),
           main = mt, pch=19,
           cex.main = 0.7, ...
      )
      la <- (k + 1):(2 * k)
      ua <- (2 * k +1):(3 * k)
      arrows(xx, mr[1:k], xx, mr[la], lty = 1, angle = 90, length = 0.10, ...)
      arrows(xx, mr[1:k], xx, mr[ua], lty = 1, angle = 90, length = 0.10, ...)
      if(length(dots) == 0) points(xx, mr[1:k], pch = 19, cex = 0.7)
      else points(xx, mr[1:k], ...)
    }
    if(PC && type == "multi"){
      if(ellipse){
        eP <- ellipse.points(m = x$pc.mean[,1:2],
                             pr = x$random.predicted.pc,
                             conf)
        plot(eP$pc1lim, eP$pc2lim, asp = 1, type = "n",
             main = mt, xlab = xlb, ylab = ylb, ...)
        for(i in 1:(dim(eP$ellP)[[3]])){
          points(eP$ellP[,,i], type = "l", ...)
        }
        if(length(dots) == 0) points(eP$means, pch=19, cex = 0.7) else
          points(eP$means, ...)
        if(label) text(eP$means, rownames(x$mean), 
                       pos=1)
      } else {
        plot(mr[1:k, 1], mr[1:k, 2], type = "n",
             xlim = xlm, ylim = ylm, main = mt,
             asp = 1, xlab = xlb, ylab = ylb, ...)
        la <- (k + 1):(2 * k)
        ua <- (2 * k +1):(3 * k)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[la, 1], mr[1:k, 2], lty = 1, length = 0.10, angle = 90, ...)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[ua, 1], mr[1:k, 2], lty = 1, length = 0.10, angle = 90, ...)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[1:k, 1], mr[la, 2], lty = 1, length = 0.10, angle = 90, ...)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[1:k, 1], mr[ua, 2], lty = 1, length = 0.10, angle = 90, ...)
        dots <- list(...)
        if(length(dots) == 0) points(mr[1:k, 1], mr[1:k, 2], pch = 19, cex = 0.7)
        else points(mr[1:k, 1], mr[1:k, 2], ...)
        if(label) text(mr[1:k, 1], mr[1:k, 2], rownames(x$mean), 
                       pos=1)
      }
    }
    
    if(!PC) {
      if(ellipse) {
        eP <- ellipse.points(m = x$mean[,1:2],
                             pr = x$random.predicted,
                             conf)
        plot(eP$pc1lim, eP$pc2lim, asp = 1, type = "n",
             main = mt, xlab = xlb, ylab = ylb, ...)
        for(i in 1:(dim(eP$ellP)[[3]])){
          points(eP$ellP[,,i], type = "l", ...)
        }
        if(length(dots) == 0) points(eP$means, pch=19, cex = 0.7) else
          points(eP$means, ...)
        if(label) text(eP$means, rownames(x$mean), 
                       pos=1)
      } else {
        plot(mr[1:k, 1], mr[1:k, 2], type = "n",
             xlim = xlm, ylim = ylm, main = mt,
             asp = 1, xlab = xlb, ylab = ylb, ...)
        la <- (k + 1):(2 * k)
        ua <- (2 * k +1):(3 * k)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[la, 1], mr[1:k, 2], lty = 1, length = 0.10, angle = 90, ...)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[ua, 1], mr[1:k, 2], lty = 1, length = 0.10, angle = 90, ...)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[1:k, 1], mr[la, 2], lty = 1, length = 0.10, angle = 90, ...)
        arrows(mr[1:k, 1], mr[1:k, 2],
               mr[1:k, 1], mr[ua, 2], lty = 1, length = 0.10, angle = 90, ...)
        if(length(dots) == 0) points(mr[1:k, 1], mr[1:k, 2], pch = 19, cex = 0.7)
        else points(mr[1:k, 1], mr[1:k, 2], ...)
        if(label) text(mr[1:k, 1], mr[1:k, 2], rownames(x$mean), 
                       pos=1)
      }
    }
  }
}

## resid and fitted

# residuals.lm.rrpp
# S3 generic for lm.rrpp

#' Extract residuals
#' 
#' @param object plot object (from \code{\link{lm.rrpp}})
#' @param weighted A logical argument to return weighted or unweighted residuals.
#' @param ... Arguments passed to other functions 
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' # See examples for lm.rrpp
residuals.lm.rrpp <- function(object, weighted = TRUE, ...) {
  if(!weighted) out <- object$LM$residuals else 
    out <- object$LM$wResiduals
  out
}

# fitted.lm.rrpp
# S3 generic for lm.rrpp

#' Extract fitted values
#' 
#' @param object plot object (from \code{\link{lm.rrpp}})
#' @param weighted A logical argument to return weighted or unweighted residuals.
#' @param ... Arguments passed to other functions 
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' # See examples for lm.rrpp
fitted.lm.rrpp <- function(object, weighted = TRUE, ...) {
  if(!weighted) out <- object$LM$fitted else 
    out <- object$LM$wFitted
  out
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{pairwise}}
#' @param ... Other arguments passed onto predict.lm.rrpp
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
#' @param object Object from \code{\link{pairwise}}
#' @param stat.table Logical argument for whether results should be returned in one table 
#' (if TRUE) or separate pairwise tables (if FALSE)
#' @param test.type Whether distances or vector correlations between vectors should be used.
#' @param angle.type If test.type = "VC", whether angle results are expressed in radians or degrees.
#' @param confidence Confidence level to use for upper confidence limit; default = 0.95 (alpha = 0.05)
#' @param ... Other arguments passed onto predict.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.pairwise <- function(object, stat.table = TRUE, 
                             test.type = c("dist", "VC"),
                             angle.type = c("rad", "deg"),
                             confidence = 0.95, ...){
  test.type <- match.arg(test.type)
  angle.type <- match.arg(angle.type)
  x <- object
  if(is.null(x$LS.means)) type = "slopes"
  if(is.null(x$slopes)) type = "means"
 
  print.pairwise(x)
  cat("\n")
  
  if(type == "means") {
    cat("LS means\n")
    print(x$LS.means[[1]])
    
    if(test.type == "dist") {
      L <- d.summary.from.list(x$means.dist)
      if(stat.table) {
        tab <- makePWDTable(L)
        cat("\nPairwise distances between means, plus statistics\n")
        print(tab)
      } else {
        cat("\nPairwise distances between means\n")
        print(L$D)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits between means\n")
        print(L$CL)
        cat("\nPairwise effect sizes (Z) between means\n")
        print(L$Z)
        cat("\nPairwise P-values between means\n")
        print(L$P)
      }
    }
    
    if(test.type == "VC") {
      L <- r.summary.from.list(x$means.vec.cor)
      if(stat.table) {
        tab <- makePWCorTable(L)
        cat("\nPairwise statistics based on mean vector correlations\n")
        if(angle.type == "deg") {
          tab$angle <- tab$angle*180/pi
          tab[,3] <- tab[,3]*180/pi
        }
        print(tab)
      } else {
        cat("\nPairwise vector correlations between mean vectors\n")
        print(L$r)
        cat("\nPairwise angles between mean vectors\n")
        if(angle.type == "deg") print(L$angle*180/pi) else print(L$angle)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits for angles between mean vectors\n")
        if(angle.type == "deg") print(L$aCL*180/pi) else print(L$aCL)
        cat("\nPairwise effect sizes (Z) for angles between mean vectors\n")
        print(L$Z)
        cat("\nPairwise P-values for angles between mean vectors\n")
        print(L$P)
      }
    }
  }
  
  if(type == "slopes") {
    cat("Slopes (vectors of variate change per one unit of covariate change, by group)\n")
    print(x$slopes[[1]])
    
    if(test.type == "dist") {
      cat("\nSlope vector lengths\n")
      print(x$slopes.length[[1]])
      L <- d.summary.from.list(x$slopes.dist)
      if(stat.table) {
        tab <- makePWDTable(L)
        cat("\nPairwise absolute difference (d) between vector lengths, plus statistics\n")
        print(tab)
      } else {
        cat("\nPairwise absolute differences (d) between slope lengths\n")
        print(L$D)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits between slope lengths\n")
        print(L$CL)
        cat("\nPairwise effect sizes (Z) between slope lengths\n")
        print(L$Z)
        cat("\nPairwise P-values between slope lengths\n")
        print(L$P)
      }
    }
    
    if(test.type == "VC") {
      L <- r.summary.from.list(x$slopes.vec.cor)
      cat("\nPairwise statistics based on slopes vector correlations (r) and angles, acos(r)")
      cat("\nThe null hypothesis is that r = 1 (parallel vectors).")
      cat("\nThis null hypothesis is better treated as 1 - r = 0, or the angle between vectors = 0\n\n")
      if(stat.table) {
        tab <- makePWCorTable(L)
        if(angle.type == "deg") {
          tab$angle <- tab$angle*180/pi
          tab[,3] <- tab[,3]*180/pi
        }
        print(tab)
      } else {
        cat("\nPairwise vector correlations between slope vectors\n")
        print(L$r)
        cat("\nPairwise angles between slope vectors\n")
        if(angle.type == "deg") print(L$angle*180/pi) else print(L$angle)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits for angles between mean vectors\n")
        if(angle.type == "deg") print(L$aCL*180/pi) else print(L$aCL)
        cat("\nPairwise effect sizes (Z) for angles between slope vectors\n")
        print(L$Z)
        cat("\nPairwise P-values for angles between slope vectors\n")
        print(L$P)
      }
    }
  }
}

