#' Evaluation of measurement error for two or more multivariate measurements, 
#' for common research subjects
#'
#' Function performs analyses concerned with the repeatability (reliability) of multivariate data 
#' (measurements) collected from the same research subjects.  Although there is no
#' requirement for repeated measurements on all research subjects, the analysis assumes
#' that multiple observations are made. 
#' 
#' This function performs analyses as described in Collyer and Adams (2024)
#'  to assess systematic and random components of 
#' measurement error (ME).  It basically performs ANOVA with RRPP,
#' but with different restricted randomization strategies.  The reliability of research subject variation 
#' can be considered by restricting randomization within replicates; the consistency of replicate measures
#' can be considered by restricting randomization within subjects.
#' Inter-subject variation remains constant across all random permutations within subjects and 
#' inter-replicate variation remains constant across all random permutations within replicates.  Type II
#' sums of squares and cross-products (SSCP) are calculated to assure conditional estimation.
#' 
#' The results include univariate-like statistics based on dispersion of values and
#' eigenanalysis performed on a signal to noise matrix product of SSCP matrices 
#' (sensu Bookstein and Mitteroecker, 2014) 
#' including the inverse of the random component of ME and the systematic
#' component of ME.  The multivariate test is a form of multivariate ANOVA (MANOVA), using
#' RRPP to generate sampling distributions of the major eigenvalue (Roy's maximum root).
#' Likelihood-ratio tests can also be performed using \code{\link{lr_test}}.
#' 
#' Intraclass correlation coefficients (ICC) can also be 
#' calculated (using \code{\link{ICCstats}}), 
#' both based on dispersion of values and 
#' covariance matrices, as descriptive statistics.  
#' Details are provided in \code{\link{ICCstats}}.
#'
#' 
#' @param data A required data frame, either of class \code{\link{data.frame}} or
#' class \code{\link{rrpp.data.frame}}.  This function cannot be used without a
#' data frame.  All arguments for data and variables are names that must exist in the data frame.
#' @param Y A name for a matrix (n x p) of data for n observations and p variables that can be found
#' in the data frame.  For example, Y = "morphData".
#' @param subjects A name for a vector or factor of research subjects, found within the data frame
#'(each subject should occur twice or more).  The length of the vector in the data frame must equal the 
#' number of observations and will be coerced into a factor.  For example, subjects = "ID".
#' @param replicates A name for a vector or factor for replicate measurements for research subjects, found
#' within the data frame.  The length of the vector  in the data frame must equal the number of observations and will 
#' be coerced into a factor.  For example, replicates = "Rep".
#' @param groups An optional name for a vector in the data frame, coercible to factor, to be included 
#' in the linear model (as an interaction with replicates).  This would be of interest if one 
#' were concerned with systematic ME occurring perhaps differently among certain strata within the data.  
#' For example, systematic ME because of an observer bias might only be observed with females or males,
#' in which case the argument might be: groups = "Sex".
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random 
#' permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found 
#' for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  
#' One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param multivariate Logical value for whether to include multivariate analyses.  Intraclass correlation 
#' matrices and relative eigenanalysis are based on products of sums of squares and cross-products (SSCP)
#' matrices, some of which must be inverted and potentially require
#'  significant computation time.  If FALSE, only statistics based on dispersion of values are calculated.
#' @param use.PCs A logical argument for whether to use the principal components of the data.  
#' This might be helpful for relative eigenanalysis, and if p > n, 
#' in which case inverting singular covariance matrices would not be possible.
#' @param tol A value indicating the magnitude below which 
#' components should be omitted, if use.PCs is TRUE. (Components are omitted if their 
#' standard deviations are less than or equal to tol times the 
#' standard deviation of the first component.)  See \code{\link{ordinate}} for more details.
#' @param Parallel The same argument as in \code{\link{lm.rrpp}} to govern parallel processing (
#' either a logical value -- TRUE or FALSE -- or the number of threaded cores to use).  See \code{\link{lm.rrpp}} 
#' for additional details.
#' @param turbo Logical value for whether to suppress coefficient estimation in RRPP iteration,
#' thus turbo-charging RRPP.
#' @param print.progress A logical value to indicate whether a progress 
#' bar should be printed to the screen.
#' @param verbose A logical value to indicate if all the output from an
#' \code{\link{lm.rrpp}} analysis should be retained.  If FALSE, only the needed
#' output for summaries and plotting is retained.
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @return Objects of class "measurement.error" return the same objects
#' as a \code{\link{lm.rrpp}} fit, plus a list of the following:
#'  \item{AOV}{Analysis of variance to test for systematic error, based on dispersion of values.}
#'  \item{mAOV}{Multivariate AOV based on product of the inverse of the random component (SSCP) of ME
#'  times the systematic component of ME.}
#'  \item{SSCP}{The sums of squares and cross-products matrices for model effects.}
#'  \item{SSCP.ME.product}{The products of the inverse of the random ME SSCP and the SSCP matrices
#'  for systematic ME,.  These are the same matrix products used for eigenanalysis.  
#'  This is the observed matrix.}
#'  \item{SSCP.ME.product.std}{A list of the symmetric forms of standardized SSCP.ME.products 
#'  that yield orthogonal eigenvectors.}

#' @references Collyer, M.L. and D.C. Adams.  2024. Interrogating Random and Systematic Measurement Error 
#' in Morphometric Data. Evolutionary Biology, 51, 179–20.
#' @references Bookstein, F.L., & Mitteroecker, P. (2014). Comparing covariance matrices by relative eigenanalysis, 
#' with applications to organismal biology. Evolutionary Biology, 41(2), 336-350.

#' @seealso \code{\link{lm.rrpp.ws}}, \code{\link{manova.update}}, 
#' \code{\link{lr_test}}

#' @examples
#' \dontrun{
#' # Measurement error analysis on simulated data of fish shapes
#' 
#' data(fishy)
#' 
#' # Example two digitization replicates of the same research subjects
#' rep1 <- matrix(fishy$coords[1,], 11, 2, byrow = TRUE)
#' rep2 <- matrix(fishy$coords[61,], 11, 2, byrow = TRUE)
#' plot(rep1, pch = 16, col = gray(0.5, alpha = 0.5), cex = 2, asp = 1)
#' points(rep2, pch = 16, col = gray(0.2, alpha = 0.5), cex = 2, asp = 1)
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
#' }
#' 
measurement.error <- function(data, 
                              Y, 
                              subjects, 
                              replicates, 
                              groups = NULL,
                              iter = 999, 
                              seed = NULL,
                              multivariate = FALSE,
                              use.PCs = TRUE, 
                              tol = 0.001, 
                              Parallel = FALSE,
                              turbo = TRUE,
                              print.progress = FALSE,
                              verbose = FALSE) {
  if(!inherits(data, "rrpp.data.frame")){
    if(!inherits(data, "data.frame") && 
       !inherits(data, "list"))
      stop(paste("\nThe data argument is neither an rrpp.data.frame",
                 "object, a data.frame object, nor a list.\n",
                 "Please see function details.\n", sep = " "), 
           call. = FALSE)
  }
  
  Y <- as.character(Y)
  Yslot <- which(names(data) %in% Y)
  if(length(Yslot) == 0)
    stop(paste("\nThe Y argument must be a character",
               "value, like 'myData',",
               "which can be found in the RRPP data frame.\n", 
               sep = " "), call. = FALSE)
  
  Y <- data[Yslot][[1]]
  if(use.PCs){
    Y <- ordinate(Y, tol = tol)$x
  }
  
  subjects <- as.character(subjects)
  Sslot <- which(names(data) %in% subjects)
  if(length(Sslot) == 0)
    stop(paste("\nThe subjects argument must be a character",
               "value, like 'mySubjects',",
               "which can be found in the RRPP data frame.\n", 
               sep = " "), call. = FALSE)
  
  subjects <- data[Sslot][[1]]
  
  replicates <- as.character(replicates)
  Rslot <- which(names(data) %in% replicates)
  if(length(Rslot) == 0)
    stop(paste("\nThe replicates argument must be a character",
               "value, like 'myReps',",
               "which can be found in the RRPP data frame.\n", 
               sep = " "), call. = FALSE)
  
  replicates <- data[Rslot][[1]]
  
  useGroups <- !is.null(groups)
  
  if(useGroups){
    
    groups <- as.character(groups)
    Gslot <- which(names(data) %in% groups)
    if(length(Gslot) == 0)
      stop(paste("\nThe groups argument must be a character",
                 "value, like 'myGroups',",
                 "which can be found in the RRPP data frame.\n", 
                 sep = " "), call. = FALSE)
    
    groups <- data[Gslot][[1]]
  }
  
  form <- if(useGroups) Y ~ subjects + groups * replicates else
    Y ~ subjects + replicates
  
  dat <- rrpp.data.frame(
    Y = Y,
    subjects = as.fctor(subjects),
    replicates = as.factor(replicates)
  )
  if(useGroups) dat$groups <- as.factor(groups)
  
  lm.rrpp.args <- list(
    f1 = form,
    data = dat,
    seed = seed,
    Parallel = Parallel,
    turbo = turbo,
    print.progress = print.progress,
    verbose = verbose,
    subjects = "subjects",
    iter = iter,
    Cov = NULL
  )
  
  fit <- suppressWarnings(do.call(lm.rrpp.ws, lm.rrpp.args))
  fit$call <- form
  data <- dat <- lm.rrpp.args <- NULL
  
  n <- fit$LM$n
  p <- fit$LM$p.prime
  
  if(multivariate && fit$LM$p > 1)
    fit <- manova.update(fit, tol = tol,
                         print.progress = print.progress,
                         verbose = verbose) 
  
  Df <- fit$ANOVA$df
  ANOVA <- getANOVAStats(fit, stat = "all")
  SNR <- ANOVA$SS / ANOVA$RSS
  AOVo <- AOV <- anova(fit, print.progress = FALSE)$table
  
  rownames(AOV)[which(
    rownames(AOV) == "subjects")] <- "Subjects"
  rownames(AOV)[which(
    rownames(AOV) == "replicates")] <- "Systematic ME"
  rownames(AOV)[which(
    rownames(AOV) == "groups")] <- "Groups"
  rownames(AOV)[which(
    rownames(AOV) == "groups:replicates")] <- 
    "Systematic ME:Groups"
  rownames(AOV)[which(
    rownames(AOV) == "Residuals")] <- "Random ME"
  
  if(multivariate){
    suppressWarnings(mAOV <- summary(fit)$stats.table)
    mAOV <- mAOV[rownames(mAOV) %in% rownames(AOVo),]
    mAOV <- mAOV[-NROW(mAOV),]
    rownames(mAOV)[which(
      rownames(mAOV) == "subjects")] <- "Subjects / Random ME"
    rownames(mAOV)[which(
      rownames(mAOV) == "replicates")] <- "Systematic ME / Random ME"
    rownames(mAOV)[which(
      rownames(mAOV) == "groups")] <- "Groups / Random ME"
    rownames(mAOV)[which(
      rownames(mAOV) == "groups:replicates")] <- 
      "Systematic ME:Groups / Random ME"
    mAOV <- mAOV[, -(1:2)]
  } else mAOV <- NULL
  
  rm(AOVo)
  rm(SNR)
  
  Xsub <- model.matrix(~ subjects)
  R <- LM.fit(Xsub, Y)$residuals
  rm(Xsub)
  SS.tot <- sum(R^2)
  rm(R)
  
  SNR <- ANOVA$SS / fit$ANOVA$RSS
  rmove <- which(Df == 0)
  if(length(rmove) > 0)
    SNR <- SNR[-rmove, ]
  
  AOV$EtaSq.ME <- NA
  AOV$EtaSq.ME[which(rownames(AOV) == "Systematic ME")] <- 
    AOV$SS[which(rownames(AOV) == "Systematic ME")] / SS.tot
  AOV$EtaSq.ME[which(rownames(AOV) == "Systematic ME:Groups")] <- 
    AOV$SS[which(rownames(AOV) == "Systematic ME:Groups")] / SS.tot
  AOV$EtaSq.ME[which(rownames(AOV) == "Random ME")] <- 
    AOV$SS[which(rownames(AOV) == "Random ME")] / SS.tot
  
  indx <- length(na.omit(AOV$F))    
  AOV$F[1:indx] <- AOV$SS[1:indx] / AOV$SS[indx + 1]
  names(AOV)[which(names(AOV) == "F")] <- "SNR"
  AOV <- AOV[c(1:4, 8, 5:7 )]
  AOV[1:indx, 7] <- apply(SNR, 1, effect.size)
  AOV[1:indx, 8] <- apply(SNR, 1, pval)
  names(AOV)[8] <- "Pr(>SNR)"
  
  out <- fit
  out$ANOVA <- ANOVA
  out$ANOVA$SS.type <- "Within-subject II"
  out$AOV <- AOV
  out$mAOV <- mAOV
  
  rm(AOV)
  rm(mAOV)
  rm(ANOVA)
  
  ### SSCP products
  
  if(p > 1){
    
    fitb <- fit
    class(fitb) <- "lm.rrpp"
    SSCP <- summary(fitb, formula = FALSE)$SSCP
    rm(fitb)
    if(length(rmove) > 0)
      SSCP <- SSCP[-rmove]
    
    SSCP.ME.product <- fast.solve(SSCP$Residuals) %*% 
      SSCP$replicates
    sscp.sqrt <- Cov.proj(SSCP$Residuals, symmetric = TRUE)
    SSCP.ME.product.std<- as.matrix(t(sscp.sqrt) %*% 
                                      SSCP$replicates %*% sscp.sqrt)
    
  } else  SSCP.ME.product <- SSCP.ME.product.std <- NULL
  
  out$SSCP <- SSCP
  rm(SSCP)
  out$SSCP.ME.product = SSCP.ME.product
  rm(SSCP.ME.product)
  out$SSCP.ME.product.std = SSCP.ME.product.std
  rm(SSCP.ME.product.std)
  rm(fit)
  
  if(!verbose) {
    out$ANOVA$MS <- out$ANOVA$Rsq <- out$ANOVA$Fs <-
      out$ANOVA$cohenf <- NULL
    out$MANOVA$SSCP <- out$MANOVA$invR.h <- NULL
    out$LM$QR <- NULL
    out$Models <- NULL
    out$PermInfo$perm.schedule <- NULL
  } 
  
  class(out) <- c("measurement.error", class(out))
  out
  
}
    
  
  