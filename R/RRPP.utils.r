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
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.lm.rrpp <- function(object, ...){
  cat("\nLinear Model fit with lm.rrpp\n")
  x <- object
  LM <- x$LM
  PI <- x$PermInfo
  AN <- x$ANOVA
  if(!is.null(x$LM$dist.coefficients)) dv <- "(dimensions of data after PCoA of distance matrix)" else
    dv <- " "
  cat(paste("\nNumber of observations:", LM$n))
  cat(paste("\nNumber of dependent variables:", LM$p, dv))
  if(!is.null(AN$SS.type)) cat(paste("\nSums of Squares and Cross-products: Type", AN$SS.type))
  cat(paste("\nNumber of permutations:", PI$perms))
  cat("\n\nFull Model Analysis of Variance\n\n")
  SS <- AN$SS
  k <- length(LM$term.labels)
  SSE <- unlist(SS[k+1,])
  SST <- unlist(SS[k+2,])
  SSM <- SST - SSE
  dfM <- LM$QR$rank
  df <- AN$df
  dfe <- df[k+1]
  Fs <- (SSM/df)/(SSE/dfe)
  P <- pval(log(Fs))
  Z <- effect.size(Fs)
  Rsq <- SSM[1]/SST[1]
  tab <- data.frame(dfM = dfM, dfe = dfe,
                    SSM = SSM[1], RSS = SSE[1], Rsq = Rsq,
                    F = Fs[1], Z = Z, P = P)
  dimnames(tab)[[2]] <- c("Df",
                     "Residual Df",
                     " SS",
                     "Residual SS",
                     "Rsq",
                     "F",
                     "Z (from F)",
                     "Pr(>F)")
  dimnames(tab)[[1]] <- deparse(x$call[[2]])
  print(tab)

  invisible(object)
}

## coef.lm.rrpp

#' Print/Summary Function for RRPP
#'
#' @param x Object from \code{\link{coef.lm.rrpp}}
#' @param ... Other arguments passed onto coef.lm.rrpp
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.coef.lm.rrpp <- function(x,...){
  cat("\nLinear Model fit with lm.rrpp\n")
  cat(paste("\nNumber of observations:", x$n))
  cat(paste("\nNumber of dependent variables:", x$p))
  cat(paste("\nSums of Squares and Cross-products: Type", x$SS.type))
  cat(paste("\nNumber of permutations:", x$nperms))
  cat("\n\nObserved coefficients\n\n")
  print(x$coef.obs)
  rrpp.type <- ifelse(x$RRPP, "RRPP", "FRPP")
  cat("\n\nStatistics (distances) of coefficients with ")
  cat(x$confidence*100, "percent confidence intervals,") 
  cat("\neffect sizes, and probabilities of exceeding observed values based on\n") 
  cat(x$nperms, "random permutations using", rrpp.type, "\n\n")
  print(x$stat.tab)
  cat("\n\n")
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



#' compare.models
#'
#' Print/Summary Function for RRPP
#'
#' @param x print/summary object (from \code{\link{lm.rrpp}})
#' @param confidence The confidence level for confidence intervals on model selection statistics
#' @param method One of AIC, adjusted R-squared (aRsq), first PC R-squared inprovement (pcRsq),
#' likelihood ratio tests (LRT), or all of these statistics for model comaprisons printed in the summary.
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.compare.models <- function(x, confidence = 0.95, 
      method = c("AIC", "aRsq", "pcRsq", "LRT", "all"), ...){
  method = match.arg(method)
  n <- x$n; p <- x$p; k <- x$k
  Dev <- x$Deviance
  Dev.obs <- unlist(lapply(1:k, function(j) (Dev[[j]])[[1]]))
  aRsq <- x$aRsq
  aRsq.obs <- aRsq[1,]
  aRsq.mean <-colMeans(aRsq)
  etamax <- max(aRsq.mean)
  delaRsq <- etamax - aRsq.mean
  aRsq.c <- etamax - aRsq
  eLCL <- apply(aRsq.c, 2, function(x) quantile(x, (1-confidence)/2))
  eUCL <- apply(aRsq.c, 2, function(x) quantile(x, confidence + (1-confidence)/2))
  Rsq.mean <- zapsmall(aRsq.mean)
  etab <- data.frame(Deviance = Dev.obs,
                     aRsq.mean = aRsq.mean, delRsq = delaRsq, 
                     LCL = eLCL, UCL = eUCL, CI.Range = eUCL - eLCL)
  
  pcRsq <- x$pcRsq
  pcRsq.obs <- pcRsq[1,]
  pcRsq.mean <-colMeans(pcRsq)
  etamax <- max(pcRsq.mean)
  delpcRsq <- etamax - pcRsq.mean
  pcRsq.c <- etamax - pcRsq
  eLCL <- apply(pcRsq.c, 2, function(x) quantile(x, (1-confidence)/2))
  eUCL <- apply(pcRsq.c, 2, function(x) quantile(x, confidence + (1-confidence)/2))
  Rsq.mean <- zapsmall(pcRsq.mean)
  ptab <- data.frame(Deviance = Dev.obs,
                     pcRsq.mean = pcRsq.mean, delRsq = delpcRsq, 
                     LCL = eLCL, UCL = eUCL, CI.Range = eUCL - eLCL)
  
  AIC <- x$AIC
  AIC.obs <- AIC[1,]
  AIC.mean <-colMeans(AIC)
  aicmin <- min(AIC.mean)
  delAIC <- AIC.mean - aicmin
  AIC.c <- AIC - aicmin
  aicLCL <- apply(AIC.c, 2, function(x) quantile(x, (1-confidence)/2))
  aicUCL <- apply(AIC.c, 2, function(x) quantile(x, confidence + (1-confidence)/2))
  AIC.mean <- zapsmall(AIC.mean)
  logL.obs <- x$logL[1,]
  AICw <- round(exp(-0.5*delAIC/p)/sum(exp(-0.5*delAIC/p)), 4)
  atab <- data.frame(Deviance = Dev.obs,
                     logL = logL.obs, AIC = AIC.mean, 
                     delAIC = delAIC, AICw = AICw, 
                     LCL = aicLCL, UCL = aicUCL, CI.Range = aicUCL - aicLCL)
  
  ltab <- x$LRT.table
  
  rownames(atab) <- rownames(etab) <- rownames(ptab) <- rownames(ltab)
  
  if(method == "AIC") {
    cat("\n\nAIC Table\n",
        "based on bootstrapped residuals\n",
        "(using residuals from each model)\n\n")
    print(atab)
  }
  
  if(method == "aRsq") {
    cat("\n\nAdjusted R-squared Table\n",
        "based on bootstrapped residuals\n",
        "(using residuals from each model)\n\n")
    print(etab)
  }
  
  if(method == "pcRsq") {
    cat("\n\nAdjusted R-squared Table\n",
        "based on bootstrapped residuals\n",
        "(using residuals from each model)\n\n")
    print(ptab)
  }
  
  if(method == "LRT") {
    cat("\n\nNon-parametric Likelihood Ratio Tests\n",
        "based on a randomized residual permutation procedure\n",
        "(using residuals from the reference model)\n\n")
    print(ltab)
  }

  if(method == "all") {
    cat("\n\nAIC Table\n",
        "based on bootstrapped residuals\n",
        "(using residuals from each model)\n\n")
    print(atab)
    
    cat("\n\nAdjusted R-squared Table\n",
        "based on bootstrapped residuals\n",
        "(using residuals from each model)\n\n")
    print(etab)
    
    cat("\n\nAdjusted R-squared Table\n",
        "based on bootstrapped residuals\n",
        "(using residuals from each model)\n\n")
    print(ptab)
    
    cat("\n\nNon-parametric Likelihood Ratio Tests\n",
        "based on a randomized residual permutation procedure\n",
        "(using residuals from the reference model)\n\n")
    print(ltab)
  }
}

#' Print/Summary Function for RRPP
#'
#' @param object print/summary object (from \code{\link{lm.rrpp}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.compare.models <- function(object, ...){
  print.compare.models(object, ...)
}

