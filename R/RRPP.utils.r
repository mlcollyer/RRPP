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
  cat(paste("\nData space dimensions:", LM$p.prime, dv))
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
  if(inherits(object, "manova.lm.rrpp")) out <- summary.manova.lm.rrpp(object) else {
    x <- object
    LM <- x$LM
    PI <- x$PermInfo
    AN <- x$ANOVA
    perms <- PI$perms
    dv <- LM$dist.coefficients
    n <- LM$n
    p <- LM$p
    p.prime <- LM$p.prime
    SS <- AN$SS
    SSM.obs <- SS[,1]
    RSS <- AN$RSS
    TSS <- AN$TSS
    RSS.model <- AN$RSS.model
    if(is.null(RSS)) RSS <- RSS.model
    if(is.null(RSS)) TSS <- RSS.model
    SS.type <- AN$SS.type
    k <- length(LM$term.labels)
    if(LM$gls) P <- LM$Pcov
    
    if(k > 0) {
      df <- AN$df
      dfe <- df[k+1]
      dfM <- sum(df[1:k])
      SSM <- (TSS[1,] - RSS.model)
      SSM.obs <- SSM[1]
      Rsq <- SSM.obs/TSS[1]
      Fs <- (SSM/dfM)/(RSS.model/dfe)
      Fs.obs <- Fs[1]
      P <- pval(Fs)
      Z <- effect.size(log(Fs))
      
      Fs.obs <- Fs[1]
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
    
    pca.fitted <- prcomp(x$LM$wFitted)
    pca.residuals <- prcomp(x$LM$wResiduals)
    pca.total <- prcomp(x$LM$Y)
    
    if(x$LM$gls) {
      Pcov <- x$LM$Pcov
      PY <- crossprod(Pcov, x$LM$Y * sqrt(x$LM$weights))
      PX <- as.matrix(crossprod(Pcov, x$LM$X * sqrt(x$LM$weights)))
      Uf <- qr.Q(qr(PX))
      int <- attr(x$LM$Terms, "intercept")
      Pint <- as.matrix(crossprod(Pcov, rep(int, n)))
      Un <- qr.Q(qr(Pint))
      glsfitted <- fastFit(Uf, PY, n, p)
      glsmeans <- fastFit(Un, PY, n, p)
      
      Sf <- crossprod(glsfitted) - crossprod(glsmeans)
      Sr <- crossprod((PY - glsfitted))
      Sy <- crossprod(PY) - crossprod(glsmeans)
      
      Cf <- Sf/(n - 1)
      Cr <- Sr/(n - 1)
      Cy <- Sy/(n - 1)
      
      pca.fitted <- pca.residuals <- pca.total <- list()
      
      svd.Y <- svd(Cy)
      keep <- which(zapsmall(svd.Y$d) > 0)
      pca.total$sdev <- sqrt(svd.Y$d[keep])
      pca.total$rotation <- as.matrix(svd.Y$v[, keep])
      pca.total$x <- (PY - glsmeans) %*% as.matrix(svd.Y$v[, keep])
      
      svd.f <- svd(Cf)
      keep <- which(zapsmall(svd.f$d) > 0)
      if(length(keep) == 0) {
        pca.fitted$sdev <- pca.fitted$rotation <- pca.fitted$x <- 0
      } else {
        pca.fitted$sdev <- sqrt(svd.f$d[keep])
        pca.fitted$rotation <- as.matrix(svd.f$v[, keep])
        pca.fitted$x <- glsfitted %*% as.matrix(svd.f$v[, keep])
      }
      
      svd.r <- svd(Cr)
      keep <- which(zapsmall(svd.r$d) > 0)
      pca.residuals$sdev <- sqrt(svd.r$d[keep])
      pca.residuals$rotation <- as.matrix(svd.r$v[, keep])
      pca.residuals$x <- (PY - glsfitted) %*% as.matrix(svd.r$v[, keep])
      
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
    eigs.f[1:rank.f] <- d.f[1:rank.f]
    eigs.r[1:rank.r] <- d.r[1:rank.r]
    eigs.t[1:rank.t] <- d.t[1:rank.t]
    eigs <- as.table(zapsmall(rbind(eigs.f, eigs.r, eigs.t)))
    rownames(eigs) <- c("Fitted", "Residuals", "Total")
    colnames(eigs) <- paste("PC", 1:rank.t, sep="")
    
    rfit <-refit(x)
    RR <- rfit$wResiduals.reduced
    RF <- rfit$wResiduals.full
    SSCP <- lapply(1:length(RF), function(j) crossprod(RR[[j]] - RF[[j]]))
    names(SSCP) <- LM$term.labels
    SSCP <- c(SSCP, list(Residuals = as.matrix(crossprod(RF[[length(RF)]]))))
    out <- list(table = tab, SSCP = SSCP, n = n, p = p, p.prime = p.prime, k = k, 
                perms = perms, dv = dv, SS = SS, SS.type = SS.type, redundancy = redund,
                eigenvalues = eigs, gls = x$LM$gls)
    class(out) <- "summary.lm.rrpp"
  }
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
  cat(paste("\nData space dimensions:", x$p.prime, dv))
  if(!is.null(x$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
  cat(paste("\nNumber of permutations:", x$perms))
  cat("\n\nFull Model Analysis of Variance\n\n")
  print(x$table)
  cat("\n\nRedundancy Analysis (PCA on fitted values and residuals)\n\n")
  if(x$gls) {
    cat("\nGLS mean used rather than center of gravity.  Projection is not orthogonal.\n\n")
  }
  print(x$redundancy)
  cat("\nEigenvalues\n\n")
  print(x$eigenvalues)
  cat("\n")
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
  cat(paste("\nData space dimensions:", x$p.prime))
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
#' @param PC Logical argument for whether to use predicted values 
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
#' indicates whether prediction line 
#' (PredLine) or regression score (RegScore) plotting is performed.  
#' For explanation of prediction line,
#' see Adams and Nistri (2010).  For explanation of regression score, see 
#' Drake and Klingenberg (2008).
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#' transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#' in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
plot.lm.rrpp <- function(x, type = c("diagnostics", "regression",
                                      "PC"), predictor = NULL,
                          reg.type = c("PredLine", "RegScore"), ...){
  plot.args <- list(...)
  r <- as.matrix(x$LM$wResiduals)
  f <- as.matrix(x$LM$wFitted)
  if(x$LM$gls) {
    r <- as.matrix(x$LM$gls.residuals)
    f <- as.matrix(x$LM$gls.fitted)
  }
  type <- match.arg(type)
  if(is.na(match(type, c("diagnostics", "regression", "PC")))) 
    type <- "diagnostics"
  PL <- Reg.proj <- PC.points <- NULL
  if(type == "diagnostics") {
    plot.args <- NULL
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
    PL <- prcomp(f)$x[,1]
    reg.type <- match.arg(reg.type)
    if(is.na(match(reg.type, c("PredLine", "RegScore")))) 
      if(is.null(predictor))
        stop("This plot type is not available without a predictor.")
    n <- NROW(r); p <- NCOL(r)
    if(!is.vector(predictor)) stop("Predictor must be a vector")
    if(length(predictor) != n) 
      stop("Observations in predictor must equal observations if procD.lm fit")
    
    plot.args$x <- predictor
    
    plot.args$ylab <- "Regression Score"
    if(is.null(plot.args$xlab)) plot.args$xlab <- deparse(substitute(predictor))

    xc <- predictor
    X <- cbind(xc, x$LM$X) * sqrt(x$LM$weights)
    b <- as.matrix(lm.fit(X, f)$coefficients)[1, ]
    Reg.proj <- center(x$LM$Y) %*% b %*% sqrt(solve(crossprod(b)))
    plot.args$y <- Reg.proj
    if(reg.type == "RegScore") {
      do.call(plot, plot.args)
    } else {
      plot.args$y <- PL
      plot.args$ylab <- "PC 1 for fitted values"
      do.call(plot, plot.args)
    }
  }
  if(type == "PC"){
    pca <- prcomp(f)
    eigs <- pca$rotation
    P <- center(x$LM$Y)%*%eigs
    v <- pca$sdev^2
    ev <- round(v[1:2]/sum(v)*100, 2)
    plot.args <- if(NCOL(P) > 1) list(x = P[,1], y = P[,2],  
                        xlab = paste("PC 1 for fitted values: ",ev[1],"%", sep = "") ,
                        ylab = paste("PC 2 for fitted values: ",ev[2],"%", sep = "") ,
                        ...) else list(x = 1:length(P), y = P, xlab = "Index", ylab = "PC 1 for fitted values: 100%",...)
    do.call(plot, plot.args)
    PC.points <- P
    rownames(P) <- rownames(x$LM$data)
  }
  out <- list(PredLine = PL, RegScore = Reg.proj, PC.points = PC.points, plot.args = plot.args)
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
#' @param x plot object (from \code{\link{predict.lm.rrpp}})
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
#' @param ... Other arguments passed onto pairwise
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
#' @param test.type Whether distances or vector correlations between vectors or variances (dispersion of residuals)
#' should be used in the test.
#' @param angle.type If test.type = "VC", whether angle results are expressed in radians or degrees.
#' @param confidence Confidence level to use for upper confidence limit; default = 0.95 (alpha = 0.05)
#' @param show.vectors Logical value to indicate whether vectors should be printed.
#' @param ... Other arguments passed onto pairwise
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.pairwise <- function(object, stat.table = TRUE, 
                             test.type = c("dist", "VC", "var"),
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
      if(stat.table) {
        tab <- makePWCorTable(L)
        if(angle.type == "deg") {
          options(warn = -1)
          tab$angle <- tab$angle*180/pi
          tab[,3] <- tab[,3]*180/pi
          options(warn = 0)
        }
      }
    }
  }
  
  if(type == "slopes") {

    if(test.type == "dist") {
      L <- d.summary.from.list(x$slopes.dist)
      if(stat.table) tab <- makePWDTable(L)
    }
    
    if(test.type == "VC") {
      L <- r.summary.from.list(x$slopes.vec.cor) 
      if(stat.table) {
        tab <- makePWCorTable(L)
        if(angle.type == "deg") {
          options(warn = -1)
          tab$angle <- tab$angle*180/pi
          tab[,3] <- tab[,3]*180/pi
          options(warn = -1)
        }
      }
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
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits between variances\n")
      print(L$CL)
      cat("\nPairwise effect sizes (Z) between variances\n")
      print(L$Z)
      cat("\nPairwise P-values between variances\n")
      print(L$P)
    }
  }
  
  if(type == "means") {
    
    cat("LS means:\n")
    if(x$show.vectors) print(x$x$LS.means[[1]]) else cat("Vectors hidden (use show.vectors = TRUE to view)\n")
    
    if(test.type == "dist") {
      if(stat.table) {
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
      if(stat.table) {
        cat("\nPairwise statistics based on mean vector correlations\n")
        print(tab)
      } else {
        cat("\nPairwise vector correlations between mean vectors\n")
        print(L$r)
        cat("\nPairwise angles between mean vectors\n")
        if(x$angle.type == "deg") print(L$angle*180/pi) else print(L$angle)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits for angles between mean vectors\n")
        if(x$angle.type == "deg") print(L$aCL*180/pi) else print(L$aCL)
        cat("\nPairwise effect sizes (Z) for angles between mean vectors\n")
        print(L$Z)
        cat("\nPairwise P-values for angles between mean vectors\n")
        print(L$P)
      }
    }
  }
  
  if(type == "slopes") {
    cat("Slopes (vectors of variate change per one unit of covariate change, by group):\n")
    if(x$show.vectors) print(x$x$slopes[[1]]) else cat("Vectors hidden (use show.vectors = TRUE to view)\n")
    
    if(test.type == "dist") {
      cat("\nSlope vector lengths\n")
      print(x$x$slopes.length[[1]])
      if(stat.table) {
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
      cat("\nPairwise statistics based on slopes vector correlations (r) and angles, acos(r)")
      cat("\nThe null hypothesis is that r = 1 (parallel vectors).")
      cat("\nThis null hypothesis is better treated as the angle between vectors = 0\n")
      if(stat.table) print(tab)
        else {
        cat("\nPairwise vector correlations between slope vectors\n")
        print(L$r)
        cat("\nPairwise angles between slope vectors\n")
        if(x$angle.type == "deg") print(L$angle*180/pi) else print(L$angle)
        cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits for angles between mean vectors\n")
        if(x$angle.type == "deg") print(L$aCL*180/pi) else print(L$aCL)
        cat("\nPairwise effect sizes (Z) for angles between slope vectors\n")
        print(L$Z)
        cat("\nPairwise P-values for angles between slope vectors\n")
        print(L$P)
      }
    }
  }
  invisible(x)
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{classify}}
#' @param ... Other arguments passed onto classify
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.classify <- function(x ,...){
  cat("\nGroups and group sizes\n")
  print(x$group.n)
  cat("\nPC means\n")
  print(x$means)
  cat("\nClassification for", length(x$class), "observations\n")
  cat("\nUse summary() to produce a table of posterior classification probabilities\n")
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{classify}}
#' @param ... Other arguments passed onto classify
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.classify <- function(object ,...){
  x <- object
  cat("\nGroups and group sizes\n")
  print(x$group.n)
  cat("\nPC means\n")
  print(x$means)
  cat("\nClassification for", length(x$class), "observations\n")
  cat("\nGeneralized (Mahalanobis) squared distances\n")
  print(x$Mah.dist.sq)
  cat("\nPrior probabilities\n")
  print(x$prior)
  cat("\nPosterior proabilities\n")
  print(x$posterior)
  cat("\n")
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{model.comparison}}
#' @param ... Other arguments passed onto model.comparison
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
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.model.comparison <- function(object ,...){
  x <- object
  type <- colnames(x$table)[[1]]
  if(type == "cov.trace") type.name <- "traces of model covariance matrices"
  if(type == "logLik") type.name <- "model log-likelihoods"
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
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.model.comparison <- function(x, ...){
  type <- colnames(x$table)[[1]]
  if(type == "cov.trace") type.name <- "Trace"
  if(type == "logLik") type.name <- "-2 * log-likelihood"
  tab <- x$table
  nms <- x$names
  x <- tab[,2]
  y <- tab[,1]
  if(type == "logLik") y <- -2 * y
  xrange <- range(x)
  dx <- xrange[2] - xrange[1]
  xlim <- xrange + c(-0.1*dx, 0.1*dx)
  
  yrange <- range(y)
  dy <- yrange[2] - yrange[1]
  ylim <- yrange + c(-0.1*dy, 0.1*dy)
  
  plot(x, y, xlab = "Parameter penalty", 
       ylab = type.name, xlim = xlim, ylim = ylim, ...)
  
  f <- lm(y ~ x)
  abline(f, lty = 3, lwd = 0.8, col = "red")
  text(x, y, nms, pos = 1, cex = 0.4)
}

#' Print/Summary Function for RRPP
#'
#' @param object Object from \code{\link{lm.rrpp}}, updated with \code{\link{manova.update}}
#' @param test Type of multivariate test statistic to use.
#' @param ... Other arguments passed onto manova.lm.rrpp
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
  perm.method <- object$PermInfo$perm.method
  if(perm.method == "RRPP") RRPP = TRUE else RRPP = FALSE
  ind <- object$PermInfo$perm.schedule
  trms <- object$LM$term.labels
  k <- length(trms)
  df <- object$ANOVA$df
  df.model <- sum(df[1:k])
  df <- c(df[1:k], df.model, df[k+1])
  names(df) <- c(trms, "Full.Model", "Residuals")
  
  MANOVA <- object$MANOVA
  if(MANOVA$verbose) {
    
    getEigs <- function(EH, EH.rank){
      r <- min(dim(na.omit(EH)))
      EH <- EH[1:r, 1:r]
      Re(eigen(EH, symmetric = FALSE, only.values = TRUE)$values)[1:EH.rank]
    }
    
    eigs <- lapply(1:(k+1), function(j){
      rh <- MANOVA$invR.H[[j]]
      rh.rank <- qr(rh[[1]])$rank
      lapply(rh, function(x) getEigs(x, rh.rank))
    })
  } else eigs <- MANOVA$eigs
  
  error <- MANOVA$error
  
  
  stats <- as.data.frame(matrix(NA, nrow = k + 2, ncol = 5, byrow = FALSE,
                                dimnames <- list(names(df), 
                                                 c("Df", "Rand", test, "Z", "Pr"))))
  stats$Df <- df
  if(!is.null(error)) stats$Rand[1:(k+1)] <- c(error, "Residuals") else stats$Rand[1:(k+1)] <- rep("Residuals", k+1)
  
  if(test == "Pillai") {
    
    rand.stats <- sapply(1:(k+1), function(j){
      y <- eigs[[j]]
      sapply(y, pillai)
    })
    colnames(rand.stats) <- c(trms, "Full.Model")
    test.stats <- rand.stats[1,]
    Z <- apply(log(rand.stats), 2, effect.size)
    P <- apply(rand.stats, 2, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Pillai[1:(k+1)] <- test.stats
  }
  else if(test == "Hotelling.Lawley") {
    
    rand.stats <- sapply(1:(k+1), function(j){
      y <- eigs[[j]]
      sapply(y, hot.law)
    })
    colnames(rand.stats) <- c(trms, "Full.Model")
    test.stats <- rand.stats[1,]
    Z <- apply(log(rand.stats), 2, effect.size)
    P <- apply(rand.stats, 2, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Hotelling.Lawley[1:(k+1)] <- test.stats
  }
  
  else if(test == "Wilks"){
    
    rand.stats <- sapply(1:(k+1), function(j){
      y <- eigs[[j]]
      sapply(y, wilks)
    })
    colnames(rand.stats) <- c(trms, "Full.Model")
    test.stats <- rand.stats[1,]
    Z <- apply(log(rand.stats), 2, effect.size)
    P <- apply(1 - rand.stats, 2, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Wilks[1:(k+1)] <- test.stats
  }
  else {
    rand.stats <- sapply(1:(k+1), function(j){
      y <- eigs[[j]]
      sapply(y, max)
    })
    colnames(rand.stats) <- c(trms, "Full.Model")
    test.stats <- rand.stats[1,]
    Z <- apply(log(rand.stats), 2, effect.size)
    P <- apply(rand.stats, 2, pval)
    stats$Z[1:(k+1)] <- Z
    stats$Pr[1:(k+1)] <- P
    stats$Roy[1:(k+1)] <- test.stats
  }
  
  if(test == "Wilks") names(stats)[[length(stats)]] <- paste("Pr(<", test, ")", sep = "") else
    names(stats)[[length(stats)]] <- paste("Pr(>", test, ")", sep = "")
  
  out <- list(stats.table = stats, rand.stats = rand.stats, stat.type = test,
              n = n, p = p, p.prime = p.prime, e.rank = MANOVA$e.rank, 
              manova.pc.dims = MANOVA$manova.pc.dims, PCA = MANOVA$PCA,
              SS.type = object$ANOVA$SS.type, perms = object$PermInfo$perms)
  class(out) <- "summary.manova.lm.rrpp"
  out
}

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{summary.manova.lm.rrpp}}
#' @param ... Other arguments passed onto summary.manova.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.summary.manova.lm.rrpp <- function(x, ...){
  
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of observations:", x$n))
  cat("\nNumber of dependent variables:", x$p)
  cat("\nData space dimensions:", x$p.prime)
  cat("\nResidual covariance matrix rank:", x$e.rank)
  pc.max <- max(x$manova.pc.dims)
  
  if(pc.max < x$e.rank) {
    cat("\n   Data reduced to", pc.max, "PCs, as required or prescribed.")
    PCA <- x$PCA
    d2 <- PCA$sdev^2
    d.p <- sum(d2[1:pc.max])/sum(d2)
    cat("\n  ", round(d.p*100, 1), "% of overall variation explained by these PCs.")
    cat("\n   See $MANOVA$PCA from manova.lm.rrpp object for more information.")
  }
    
  if(!is.null(x$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
  cat(paste("\nNumber of permutations:", x$perms), "\n\n")
  
  tab <- as.matrix(x$stats.table)
  print.table(tab, na.print = "")
  invisible(x)
}


#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{trajectory.analysis}}
#' @param ... Other arguments passed onto 
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
#' @param stat.table Logical argument for whether results should be returned in one table 
#' (if TRUE) or separate pairwise tables (if FALSE)
#' @param attribute Whether magnitude differences (MD, absolute difference in trajectory path lengths), 
#' trajectory correlations (TC), or trajectory shape differences (SD) are summarized.
#' @param angle.type If attribute = "TC", whether angle results are expressed in radians or degrees.
#' @param confidence Confidence level to use for upper confidence limit; default = 0.95 (alpha = 0.05)
#' @param show.trajectories Logical value to indicate whether trajectories should be printed.
#' @param ... Other arguments passed onto trajectory.analysis
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
    L <- r.summary.from.list(TC, confidence = confidence)
    tab <- makePWCorTable(L)
    if(angle.type == "deg") {
      options(warn = -1)
      tab$angle <- tab$angle*180/pi
      tab[,3] <- tab[,3]*180/pi
      options(warn = -1)
    }
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
    if(x$show.trajectories) print(x$x$trajectories[[1]]) else cat("Trajectories hidden (use show.trajectories = TRUE to view)\n")
    
    cat("\nObserved path distances by group\n\n")
    print(x$x$PD[[1]])
    
    if(stat.table) {
      cat("\nPairwise absolute differences in path distances, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise absolute differences in path distancess\n")
      print(L$D)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits, absolute differences in path distancess\n")
      print(L$CL)
      cat("\nPairwise effect sizes (Z) for absolute differences in path distancess\n")
      print(L$Z)
      cat("\nPairwise P-values for absolute differences in path distances\n")
      print(L$P)
    }
  }
  
  if(attribute == "TC") {
    
    cat("Trajectories:\n")
    if(x$show.trajectories) print(x$x$trajectories[[1]]) else cat("Trajectories hidden (use show.trajectories = TRUE to view)\n")
    
    if(stat.table) {
      cat("\nPairwise correlations between trajectories, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise statistics based on trajectory vector correlations (r) and angles, acos(r)")
      cat("\nThe null hypothesis is that r = 1 (parallel vectors).")
      cat("\nThis null hypothesis is better treated as the angle between vectors = 0\n")
      
      cat("\nPairwise vector correlations between trajectories\n")
      print(L$r)
      cat("\nPairwise angles between trajectories\n")
      if(x$angle.type == "deg") print(L$angle*180/pi) else print(L$angle)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits for angles\n")
      if(x$angle.type == "deg") print(L$aCL*180/pi) else print(L$aCL)
      cat("\nPairwise effect sizes (Z) for angles\n")
      print(L$Z)
      cat("\nPairwise P-values for angles\n")
      print(L$P)
    }
  }
  
  if(type == "trajectories"  && attribute == "SD") {
    
    cat("Trajectories:\n")
    if(x$show.trajectories) print(x$x$trajectories[[1]]) else cat("Trajectories hidden (use show.trajectories = TRUE to view)\n")

    if(stat.table) {
      cat("\nPairwise trajectory shape differences, plus statistics\n")
      print(tab)
    } else {
      cat("\nPairwise trajectory shape differences\n")
      print(L$D)
      cat("\nPairwise", paste(L$confidence*100, "%", sep=""), "upper confidence limits, shape differnces\n")
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
#'  \code{\link{lm.rrpp}} that are passed onto \code{\link{trajectory.analysis}}, and projects
#'  data onto them.  This function is a set.up, and \code{\link{add.trajectories}} is needed to 
#'  add trajectories to the plot.  By having two stages of control, the plotting functions are more 
#'  flexible.  This function also returns plotting information that can be valuable for making
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
#'  \item{trajectories}{pca Observed trajectories projected onto principal components.}
#'  
#' @seealso 
#' \code{\link{plot.default}} and \code{\link{par}}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Adams, D. C., and M. M. Cerney. 2007. Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @references Collyer, M. L., and D. C. Adams. 2013. Phenotypic trajectory analysis: comparison of shape change patterns 
#' in evolution and ecology. Hystrix 24: 75-83.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
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
    f <- if(x$fit$LM$gls) x$fit$LM$gls.fitted else x$fit$LM$wFitted
    pca <- prcomp(f)
    rot <- pca$rotation
    Y <- x$fit$LM$Y
    Y.cent <- colMeans(Y)
    props <- pca$sdev^2 / sum(pca$sdev^2)
    pc.points <- center(Y) %*% rot
    trajectories <- x$trajectories[[1]]
    if(is.matrix(trajectories)) trajectories <- list(trajectories)
    traj.c <- matrix(Y.cent, NROW(trajectories[[1]]), NCOL(trajectories[[1]]), byrow = TRUE)
    trajectories <- lapply(trajectories, function(x) (x - traj.c) %*% rot)
  }
  
  if(x$type == "single.factor") {
    f <- if(x$fit$LM$gls) x$fit$LM$gls.fitted else x$fit$LM$wFitted
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
    xlabel <- paste("PC 1 for fitted values: ", round(props[1] *100, 2), "%", sep = "")
  if(is.null(dots$ylab))
    ylabel <- paste("PC 2 for fitted values: ", round(props[2] *100, 2), "%", sep = "")
  
  if(!is.null(dots$xlab) && !is.null(dots$ylab)) plot(pc.points[,1], pc.points[,2], asp = 1, ...)
  if(!is.null(dots$xlab) && is.null(dots$ylab)) plot(pc.points[,1], pc.points[,2], asp = 1, ylab = ylabel, ...)
  if(is.null(dots$xlab) && !is.null(dots$ylab)) plot(pc.points[,1], pc.points[,2], asp = 1, xlab = xlabel, ...)
  if(is.null(dots$xlab) && is.null(dots$ylab))
    plot(pc.points[,1], pc.points[,2], asp = 1, xlab = xlabel, ylab = ylabel, ...)
  out <- list(pca = pca, pc.points = pc.points, trajectoy.analysis = x, trajectories = trajectories)
  invisible(out)
}

#' Plot Function for RRPP
#' 
#'  Function adds trajectories to a principal component plot
#'
#'  The function adds trajectories to a plot made by \code{\link{plot.trajectory.analysis}}.
#'  This function has a restricted set of plot parameters based on the number of trajectories
#'  to be added to the plot.
#'  
#' @param TP plot object (from \code{\link{plot.trajectory.analysis}})
#' @param traj.pch Plotting "character" for trajectory points.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for pch.
#' @param traj.col The color of trajectory lines.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for col.
#' @param traj.lty Trajectory line type.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for lty.
#' @param traj.lwd Trajectory line width.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for lwd.
#' @param traj.cex Trajectory point character expansion.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for cex.
#' @param traj.bg Trajectory point background.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for bg.
#' @param start.bg Trajectory point background, just the start points.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for bg.  Green start points are the default.
#' @param end.bg Trajectory point background, just the end points.  Can be a single value or vector 
#' of length equal to the number of trajectories.  See \code{\link{par}} and its description 
#' for bg.  Red end points are the default.
#' 
#' @seealso 
#' \code{\link{plot.default}} and \code{\link{par}}
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @references Adams, D. C., and M. M. Cerney. 2007. Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @references Collyer, M. L., and D. C. Adams. 2013. Phenotypic trajectory analysis: comparison of shape change patterns 
#' in evolution and ecology. Hystrix 24: 75-83.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
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
    stop("For add.trajectories, traj.pch must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.pch) == 1) traj.pch <- rep(traj.pch, nt)
  
  if(length(traj.col) != 1 && length(traj.col) != nt)
    stop("For add.trajectories, traj.col must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.col) == 1) traj.col <- rep(traj.col, nt)
  
  if(length(traj.lty) != 1 && length(traj.lty) != nt)
    stop("For add.trajectories, traj.lty must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.lty) == 1) traj.lty <- rep(traj.lty, nt)
  
  if(length(traj.lwd) != 1 && length(traj.lwd) != nt)
    stop("For add.trajectories, traj.lwd must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.lwd) == 1) traj.lwd <- rep(traj.lwd, nt)
  
  if(length(traj.cex) != 1 && length(traj.cex) != nt)
    stop("For add.trajectories, traj.cex must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.cex) == 1) traj.cex <- rep(traj.cex, nt)
  
  if(length(traj.bg) != 1 && length(traj.bg) != nt)
    stop("For add.trajectories, traj.bg must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(traj.bg) == 1) traj.bg <- rep(traj.bg, nt)
  
  if(length(start.bg) != 1 && length(start.bg) != nt)
    stop("For add.trajectories, start.bg must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(start.bg) == 1) start.bg<- rep(start.bg, nt)
  
  if(length(end.bg) != 1 && length(end.bg) != nt)
    stop("For add.trajectories, end.bg must be equal in length to the number of trajectories or just one value\n",
         call. = FALSE) else if(length(end.bg) == 1) end.bg<- rep(end.bg, nt)     
  
  for(i in 1:nt){
    x <- traj[[i]][,1]
    y <- traj[[i]][,2]
    lines(x, y, col = traj.col[i], lwd = traj.lwd[i], lty = traj.lty[i])
    points(x, y, col = 1, pch = traj.pch[i], lwd = traj.lwd[i], cex = traj.cex[i], bg = traj.bg[i])
    points(x[1], y[1], col = 1, pch = traj.pch[i], lwd = traj.lwd[i], cex = traj.cex[i], bg = start.bg[i])
    points(x[np], y[np], col = 1, pch = traj.pch[i], lwd = traj.lwd[i], cex = traj.cex[i], bg = end.bg[i])
  }
  
}

#' Plot Function for RRPP
#' 
#' @param x An object of class \code{\link{ordinate}}
#' @param axis1 A value indicating which component should be displayed as the X-axis (default = C1)
#' @param axis2 A value indicating which component should be displayed as the Y-axis (default = C2)
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' @return An object of class "plot.ordinate" is a list with components
#'  that can be used in other plot functions, such as the type of plot, points, 
#'  a group factor, and other information depending on the plot parameters used.
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.ordinate <- function(x, axis1 = 1, axis2 = 2, ...) {
  options(warn = -1)
  if(NCOL(x$x) == 1) stop("Only one component  No plotting capability with this function.\n", 
                          call. = FALSE)
  v <- x$d/sum(x$d)
  plot.args <- list(x = x$x[, axis1], y = x$x[, axis2],  ...)
  xlabel <- paste("C ", axis1, ": ", round(v[axis1] * 100, 2), "%", sep = "")
  ylabel <- paste("C ", axis2, ": ", round(v[axis2] * 100, 2), "%", sep = "")
  if(is.null(plot.args$xlab)) plot.args$xlab <- xlabel
  if(is.null(plot.args$ylab)) plot.args$ylab <- ylabel
  pcdata <- as.matrix(x$x[, c(axis1, axis2)])
  if(!is.null(plot.args$axes)) axes <- plot.args$axes else axes <- TRUE
  if(!is.logical(axes)) axes <- as.logical(axes)
  plot.args$xlim <- 1.05*range(plot.args$x)
  plot.args$ylim <- 1.05*range(plot.args$y)
  if(is.null(plot.args$asp)) plot.args$asp <- 1
  if(is.null(plot.args$phylo)) plot.args$phylo <- FALSE
  
  do.call(plot, plot.args)
  
  if(axes){
    abline(h = 0, lty=2, ...)
    abline(v = 0, lty=2, ...)
  }
  
  options(warn = 0)
  out <- list(points = pcdata,   
              call = match.call())
  
  out$plot.args <- plot.args
  class(out) <- "plot.ordinate"
  invisible(out)
  
}

