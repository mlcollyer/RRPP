#' ANOVA for lm.rrpp model fits
#'
#' @description Computes an analysis of variance table using distributions of random
#' statistics from \code{\link{lm.rrpp}}.
#'
#' @param object Object from \code{\link{lm.rrpp}}
#' @param ... Other arguments passed to ANOVA
#' @param effect.type One of "F", "cohenf", "SS", "MS", "Rsq" to choose from which distribution
#' of statistics to calculate effect sizes (Z).  See \code{\link{lm.rrpp}}.
#' @param error An optional character string to define MS error term for calculation of F values.
#' See \code{\link{lm.rrpp}} for examples.
#' @export
#' @author Michael Collyer
#' @keywords utilities
anova.lm.rrpp <- function(object, ...,
                          effect.type = c("F", "cohenf", "SS", "MS", "Rsq"),
                          error = NULL) {
  x <- object$ANOVA
  df <- x$df
  k <- length(df)-2
  SS <- x$SS
  MS <- x$MS
  if(!is.null(error)) {
    if(!inherits(error, "character")) stop("The error description is illogical.  It should be a string of character values matching ANOVA terms.")
    kk <- length(error)
    if(kk != k) stop("The error description should match in length the number of ANOVA terms (not including Residuals)")
    trms <- rownames(MS)
    MSEmatch <- match(error, trms)
    if(any(is.na(MSEmatch))) stop("At least one of the error terms is not an ANOVA term")
  } else MSEmatch <- NULL
  if(length(SS) > 1) {
    Fs <- x$Fs
    if(!is.null(MSEmatch)){
      Fs[1:k,] <- MS[1:k,]/MS[MSEmatch,]
      F.effect.adj <- apply(Fs[1:k,], 1, effect.size)
    }
    effect.type <- match.arg(effect.type)
    if(effect.type == "F") {
      if(!is.null(error)) Z <- F.effect.adj else Z <- x$F.effect
    }
    if(object$LM$gls) {
      est <- "Generalized Least-Squares (via OLS projection)"
      if(effect.type == "SS") {
        cat("\nWarning: calculating effect size on SS is illogical with GLS.
            Effect type has been changed to F distributions.\n\n")
        effect.type = "F"
      }
      if(effect.type == "MS") {
        cat("\nWarning: calculating effect size on MS is illogical with GLS.
            Effect type has been changed to F distributions.\n\n")
        effect.type = "F"
      }
    }
    else
      est <- "Ordinary Least Squares"
    if(effect.type == "F") Z <- Fs
    if(effect.type == "SS") Z <- x$SS
    if(effect.type == "MS") Z <- x$MS
    if(effect.type == "Rsq") Z <- x$Rsq
    if(effect.type == "cohenf") Z <- x$cohenf
    if(effect.type == "Rsq") effect.type = "R-squared"
    if(effect.type == "cohenf") effect.type = "Cohen's f-squared"
    if(!is.matrix(Z)) Z <- matrix(Z, 1, length(Z))
    Fs <- Fs[,1]
    SS <- SS[,1]
    MS <- MS[,1]
    MS[length(MS)] <- NA
    Rsq <- x$Rsq
    cohenf <- x$cohenf
    P.val <- apply(Z, 1, pval)
    Z <- apply(Z, 1, effect.size)
    P.val[-(1:k)] <- NA
    Z[-(1:k)] <- NA
    Rsq <- Rsq[,1]
    Rsq[length(Rsq)] <- NA
    tab <- data.frame(Df=df, SS=SS, MS = MS, Rsq = Rsq, F=Fs, Z=Z, P.val=P.val)
    colnames(tab)[NCOL(tab)] <- paste("Pr(>", effect.type, ")", sep="")
    class(tab) = c("anova", class(tab))
    if(object$LM$gls)
      est <- "Generalized Least-Squares (via OLS projection)" else
        est <- "Ordinary Least Squares"
    perms <- object$PermInfo$perms
    pm <- object$PermInfo$perm.method
    if(pm == "RRPP") pm <- "Randomization of null model residuals" else
      pm <- ("Randomization of raw values (residuals of means)")
    SS.type <- x$SS.type
    cat("\nAnalysis of Variance, using Residual Randomization\n")
    cat(paste("Permutation procedure:", pm, "\n"))
    cat(paste("Number of permutations:", perms, "\n"))
    cat(paste("Estimation method:", est, "\n"))
    cat(paste("Sums of Squares and Cross-products: Type", SS.type, "\n"))
    cat(paste("Effect sizes (Z) based on", effect.type, "distributions\n\n"))
    print(tab)
    cat("\nCall: ")
    cat(deparse(object$call), fill=TRUE)
      } else {
        cat("\nANOVA Table\n\n")
        Residuals <- c(df, SS, MS)
        names(Residuals) <- c("Df", "SS", "MS")
        tab <- as.data.frame(Residuals)
        class(tab) = c("anova", class(tab))
        print(tab)
        cat("\nCall: ")
        cat(deparse(object$call), fill=TRUE)
      }
}