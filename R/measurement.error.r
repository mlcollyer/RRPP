#' Evaluation of measurement error for two or more multivariate measurements, 
#' for common research subjects
#'
#' Function performs analyses concerned with the repeatability (reliability) of multivariate data 
#' (measurements) collected from the same research subjects.  Although there is no
#' requirement for repeated measurements on all research subjects, the analysis assumes
#' that multiple observations are made. 
#' 
#' This function performs analyses as described in XXX to assess systematic and random components of 
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
#' component of ME.  The multivriate test is a permutational form of multivariate ANOVA (MANOVA).
#' 
#' The first eigenvalue (Roy's maximum root) is used as a test statistic 
#' for tests using MANOVA.
#' Intraclass correlation coefficients (ICC) are also calculated, both based on dispersion of values and 
#' covariance matrices, as descriptive statistics.  Multivariate generalizations of the statistics
#' described by Liljequist et al. (2019) are also used, along with eigenanalysis.  
#' Three statistics describe the ICC for the population,
#' agreement of measurements among subjects, and consistency between measurements.  The last statistic does not 
#' necessarily measure the sameness between measurements but the consistency of change between measurements,
#' which might be indicative of a systematic measurement error.  If groups are used, these three statistics are 
#' repeated, using the SSCP for groups-adjusted data.  This approach accounts for group differences,
#' which would avoid large subject variation compared to measurement error inflating ICC values.  If there are 
#' inherently disparate groups from which subjects are sampled, this approach can elucidate better agreement and 
#' consistency in light of group differences.
#'
#'   
#' @param Y A matrix (n x p) of data for n observations and p variables.
#' @param subj A vector or factor of research subjects (each subject should occur twice or more).  
#' The length of the vector must equal the number of observations and will be coerced into a factor.
#' @param reps A vector or factor for replicate measurements for research subjects.  
#' The length of the vector must equal the number of observations and will be coerced into a factor.
#' @param groups An optional vector, coercible to factor, to be included in the linear model
#' (as an interaction with replicates)..
#' This would be of interest if one were concerned with systematic ME occurring perhaps differently among 
#' certain strata within the data.  For example, systematic ME because of an observer bias might
#' only be observed with females or males.  
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
#' components should be omitted., if use.PCs is TRUE. (Components are omitted if their 
#' standard deviations are less than or equal to tol times the 
#' standard deviation of the first component.)  See \code{\link{ordinate}} for more details.
#' @param Parallel The same argument as in \code{\link{lm.rrpp}} to govern parallel processing (
#' either a logical vale -- TRUE or FALSE -- or the number of threaded cores to use).  See \code{\link{lm.rrpp}} 
#' for additional details.
#' @param print.progress A logical value to indicate whether a progress 
#' bar should be printed to the screen.
#' @param verbose A logical value to indicate if all the output from an
#' \code{\link{lm.rrpp}} analysis should be retained.  If FALSE, only the needed
#' output for summaries and plotting is retained.
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @return Objects of class "measurement.error" return a list of the following:
#'  \item{AOV}{Analysis of variance to test for systematic error, based on dispersion of values.}
#'  \item{mAOV}{Multivariate AOV based on product of the inverse of the random component (SSCP) of ME
#'  times the systematic component of ME.}
#'  \item{icc}{The intraclass correlation coefficient (ICC) based on the dispersion of values.}
#'  \item{mult.icc.eigs}{The eigenvalues of ICC matrices, culled to principal dimensions with positive eigenvalues.}
#'  \item{SSCP}{The sums of squares and cross-products matrices for model effects.}
#'  \item{SSCP.ME.product}{The products of the inverse of the random ME SSCP and the SSCP matrices
#'  for systematic ME,.  These are the same matrix products used for eigenanalysis.  
#'  This is the observed matrix.}
#'  \item{SSCP.ME.product.std}{A list of the symmetric forms of standradized SSCP.ME.products 
#'  that yield orthogonal eigenvectors.}
#'  \item{all.stats}{All SS, MS, eigen values, etc., from the RRPP analyses performed.  This is the same
#'  as the output found in an \code{\link{lm.rrpp}} object, updated with \code{\link{manova.update}}.
#'  This object only contains the many RRPP ANOVA and MANOVA statistics if verbose = TRUE.}

#' @references  yet to be determined.
#' @references Bookstein, F. L., & Mitteroecker, P. (2014). Comparing covariance matrices by relative eigenanalysis, 
#' with applications to organismal biology. Evolutionary biology, 41(2), 336-350.
#' @references Liljequist, D., Elfving, B., & Skavberg Roaldsen, K. (2019). Intraclass correlation–A discussion 
#' and demonstration of basic features. PloS one, 14(7), e0219854.
#' @seealso \code{\link{lm.rrpp}}, \code{\link{manova.update}}

#' @examples
#' # TBD

measurement.error <- function(Y, 
                              subj, 
                              reps, 
                              groups = NULL,
                              iter = 999, 
                              seed = NULL,
                              multivariate = FALSE,
                              use.PCs = TRUE, 
                              tol = 0.001, 
                              Parallel = TRUE,
                              print.progress = FALSE,
                              verbose = FALSE) {

  wrn <- options()$warn
  options(warn = -1)
  
  Y <- if(use.PCs) as.matrix(ordinate(as.matrix(Y), tol = tol)$x) else 
    as.matrix(Y)
  
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  
  if(p == 1) multivariate = FALSE
  
  if(length(subj) != n)
    stop("The number of observations do not match for data and subjects.\n",
         call. = FALSE)
  
  if(length(reps) != n)
    stop("The number of observations do not match for data and reps\n",
         call. = FALSE)
  
  subj <- as.factor(subj)
  reps <- as.factor(reps)
  ns <- nlevels(subj)
  nr <- nlevels(reps)
  levls <- levels(subj)
  levlr <- levels(rep)
  
  subj.check <- by(subj, subj, length)
  if(all(subj.check <= 1))
    stop("No replication of observations for any subjects.\n",
         call. = FALSE)
  
  rep.check <- by(reps, reps, length)
  if(all(rep.check <= 1))
    stop("No subjects have two or more reps.\n",
         call. = FALSE)
  
  dat <- rrpp.data.frame(Y = Y, subj = subj, reps = reps)
  
  if(print.progress) {
    cat("\nPerforming within-subject RRPP.\n")
  }
  
  if(!is.null(groups)) {
    if(length(groups) != n)
      stop("The number of observations do not match for data and groups factor\n",
           call. = FALSE)
    groups <- as.factor(groups)
    
    Xg <- model.matrix(~ groups)
    Xnull <- matrix(1, length(subj))
    dat$Ya <- LM.fit(Xg, Y)$residuals + LM.fit(Xnull, Y)$fitted.values
  }
  
  form <- if(!is.null(groups)) Y ~ subj + reps * groups else
    Y ~ subj + reps
  
  # within-subject analysis
  
  suppressWarnings(fit <- lm.rrpp(form, 
                 SS.type = "II", block = subj,
                 data = dat, iter = iter, 
                 Parallel = Parallel,
                 print.progress = print.progress, turbo = TRUE,
                 seed = seed, verbose = verbose))
  fit$call[[2]] <- form
  
  if(print.progress) {
    cat("\nPerforming within-replicate RRPP.\n")
  }
  
  # within replicate analysis
  
  suppressWarnings(fit2 <- lm.rrpp(form, 
                 SS.type = "II", block = reps,
                 data = dat, iter = iter, 
                 Parallel = Parallel,
                 print.progress = print.progress, turbo = TRUE,
                 seed = seed))
  
  ANOVA <- getANOVAStats(fit, stat = "all")
  ANOVA2 <- getANOVAStats(fit2, stat = "all")
  ANOVA$SS[1,] <- ANOVA2$SS[1,]
  ANOVA$MS[1,] <- ANOVA2$MS[1,]
  ANOVA$Rsq[1,] <- ANOVA2$Rsq[1,]
  ANOVA$Fs[1,] <- ANOVA2$Fs[1,]
  ANOVA$cohenf[1,] <-  ANOVA2$cohenf[1,]
  
  fit$ANOVA$RSS[1,] <- fit2$ANOVA$RSS[1,]
  fit$ANOVA$TSS[1,] <- fit2$ANOVA$TSS[1,]
  fit$ANOVA$RSS.model[1,] <- fit2$ANOVA$RSS.model[1,]
  fit$ANOVA$SS <- ANOVA$SS
  fit$ANOVA$MS <- ANOVA$MS
  fit$ANOVA$Rsq <- ANOVA$Rsq
  fit$ANOVA$Fs <- ANOVA$Fs
  fit$ANOVA$cohenf <-  ANOVA$cohenf
  ANOVA <- ANOVA2 <- NULL
  
  mAOV <- micc <- NULL
  
  if(multivariate) {
    if(print.progress) {
      cat("\nFor the random linear model fits performed in the previous step.\n")
      cat("The inverse of the SSCP matrix for random measurement error will be multiplied\n")
      cat("by the SSCP matrix for subjects and replicates, for performing realtive eigenanalysis.\n")
    }
    
    suppressWarnings(fitm <- if(p > 1)
      manova.update(fit, print.progress = print.progress) else NULL)
    
    suppressWarnings(fitm2 <- if(p > 1) 
                       manova.update(fit2, print.progress = print.progress) else NULL)
    if(p > 1) {
      for(i in 1:fitm$PermInfo$perms) {
        fitm$MANOVA$invR.H[[i]][[1]] <- fitm2$MANOVA$invR.H[[i]][[1]]
        fitm$MANOVA$eigs[[i]][1,] <- fitm2$MANOVA$eigs[[i]][1,] 
      }
      rm(fitm2)
    }
  }
  
  rm(fit2)
  
  
  # ICC Dispersion
  
  fita <- NULL
  
  if(!is.null(groups)) {
    fitf <- fit
    fit <- fitr <- lm.rrpp(Y ~ subj + reps, data = dat, 
                   SS.type = "II", iter = 0, print.progress = 0)
    
    ANOVA <- fit$ANOVA
    ANOVA2 <- getANOVAStats(fit, stat = "all")
    ANOVA[c("MS", "Rsq", "Fs", "cohenf")] <-
      ANOVA2[c("MS", "Rsq", "Fs", "cohenf")]
    
    MS <- ANOVA$MS[, 1]
    Df <- ANOVA$df
    MSBS <- MS[1]
    MSBM <- MS[2]
    RSS <- ANOVA$RSS[1]
    MSWS <- (MS[2] * Df[2] + RSS) / sum(Df[c(2:3)])
    MSE <- RSS / Df[3]
    icc1 <- (MSBS - MSWS) / (MSBS + (nr - 1) * MSWS)
    icca1 <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE + 
                               nr / ns * (MSBM - MSE))
    iccc1 <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE)
    
    fit <- fita <- lm.rrpp(Ya ~ subj + reps, data = dat, 
                   SS.type = "II", iter = 0, print.progress = 0)
    
    ANOVA <- fit$ANOVA
    ANOVA2 <- getANOVAStats(fit, stat = "all")
    ANOVA[c("MS", "Rsq", "Fs", "cohenf")] <-
      ANOVA2[c("MS", "Rsq", "Fs", "cohenf")]
    
    MS <- ANOVA$MS[, 1]
    Df <- ANOVA$df
    MSBS <- MS[1]
    MSBM <- MS[2]
    RSS <- ANOVA$RSS[1]
    MSWS <- (MS[2] * Df[2] + RSS) / sum(Df[c(2:3)])
    MSE <- RSS / Df[3]
    icc2 <- (MSBS - MSWS) / (MSBS + (nr - 1) * MSWS)
    icca2 <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE + 
                               nr / ns * (MSBM - MSE))
    iccc2 <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE)
    fit <- fitf
    if(!multivariate) rm(fitr, fitf)
    
  } else {
    
    ANOVA <- fit$ANOVA
    ANOVA2 <- getANOVAStats(fit, stat = "all")
    ANOVA[c("MS", "Rsq", "Fs", "cohenf")] <-
      ANOVA2[c("MS", "Rsq", "Fs", "cohenf")]
    
    MS <- ANOVA$MS[, 1]
    Df <- ANOVA$df
    MSBS <- MS[1]
    MSBM <- MS[2]
    RSS <- ANOVA$RSS[1]
    MSWS <- (MS[2] * Df[2] + RSS) / sum(Df[c(2:3)])
    MSE <- RSS / Df[3]
    icc1 <- (MSBS - MSWS) / (MSBS + (nr - 1) * MSWS)
    icca1 <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE + 
                               nr / ns * (MSBM - MSE))
    iccc1 <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE)
    icc2 <- icca2 <- iccc2 <- NULL
  }
  
  icc <- list(icc = icc1, icc.a = icca1, icc.c = iccc1,
              icc.g = icc2, icc.a.g = icca2, icc.c.g = iccc2)
  icc1 <- icc.a <- icc.c <- icc.g <- icc.a.g <- icc.c.g <- NULL
  # ICC (multivariate) need to generalized
  
  if(multivariate){
    if(!is.null(groups)) {
      fit <- fitr
      S <- summary(fit)
      SSCP <- S$SSCP
      RSSCP <- as.matrix(SSCP[[length(SSCP)]])
      SSCPBS <- as.matrix(SSCP[[1]])
      SSCPBM <- as.matrix(SSCP[[2]])
      SSCPWS <- RSSCP + SSCPBM
      Df <- fit$ANOVA$df
      MSCPBS <- SSCPBS / Df[1]
      MSCPBM <- SSCPBM / Df[2]
      MSCPWS <- SSCPWS / sum(Df[2:3])
      MSCPE <- RSSCP / Df[3]
      
      micc1 <- fast.solve(MSCPBS + (nr - 1) * MSCPWS) %*% (MSCPBS - MSCPWS) 
      micca1 <- fast.solve(MSCPBS + (nr - 1) * MSCPE + 
                             nr / ns * (MSCPBM - MSCPE)) %*% (MSCPBS - MSCPE) 
      miccc1 <- solve(MSCPBS + (nr - 1) * MSCPE) %*% (MSCPBS - MSCPE) 
      
      fit <- fita
      S <- summary(fit)
      SSCP <- S$SSCP
      RSSCP <- as.matrix(SSCP[[length(SSCP)]])
      SSCPBS <- as.matrix(SSCP[[1]])
      SSCPBM <- as.matrix(SSCP[[2]])
      SSCPWS <- RSSCP + SSCPBM
      Df <- fit$ANOVA$df
      MSCPBS <- SSCPBS / Df[1]
      MSCPBM <- SSCPBM / Df[2]
      MSCPWS <- SSCPWS / sum(Df[2:3])
      MSCPE <- RSSCP / Df[3]

      micc2 <- fast.solve(MSCPBS + (nr - 1) * MSCPWS) %*% (MSCPBS - MSCPWS) 
      micca2 <- fast.solve(MSCPBS + (nr - 1) * MSCPE + 
                             nr / ns * (MSCPBM - MSCPE)) %*% (MSCPBS - MSCPE) 
      miccc2 <- fast.solve(MSCPBS + (nr - 1) * MSCPE) %*% (MSCPBS - MSCPE) 
      
      fit <- fitf
      rm(fitr, fitf, fita)
      
    } else {
      
      S <- summary(fit)
      SSCP <- S$SSCP
      RSSCP <- as.matrix(SSCP[[length(SSCP)]])
      SSCPBS <- as.matrix(SSCP[[1]])
      SSCPBM <- as.matrix(SSCP[[2]])
      SSCPWS <- RSSCP + SSCPBM
      Df <- fit$ANOVA$df
      MSCPBS <- SSCPBS / Df[1]
      MSCPBM <- SSCPBM / Df[2]
      MSCPWS <- SSCPWS / sum(Df[2:3])
      MSCPE <- RSSCP / Df[3]
      
      micc1 <- fast.solve(MSCPBS + (nr - 1) * MSCPWS) %*% (MSCPBS - MSCPWS) 
      micca1 <- fast.solve(MSCPBS + (nr - 1) * MSCPE + 
                             nr / ns * (MSCPBM - MSCPE)) %*% (MSCPBS - MSCPE) 
      miccc1 <- fast.solve(MSCPBS + (nr - 1) * MSCPE) %*% (MSCPBS - MSCPE) 
      
      micc2 <- micca2 <- miccc2 <- NULL
    }
    
    micc <- list(icc = micc1, icc.a = micca1, icc.c = miccc1,
                 icc.g = micc2, icc.a.g = micca2, icc.c.g = miccc2)
    micc1 <- micca1 <- micc1 <- micc2 <- micca2 <- micc2 <- NULL
  } else micc <- NULL
  
  # clear memory
  rm(dat)
  
  # ANOVA 
  
  suppressWarnings(AOV <- anova(fit, effect.type = "Rsq")$table)
  if(multivariate) {
    suppressWarnings(mAOV <- summary(fitm)$stats.table)
    mAOV <- mAOV[rownames(mAOV) %in% rownames(AOV),]
    mAOV <- mAOV[-NROW(mAOV),]
  }

  rownames(AOV)[which(rownames(AOV) == "subj")] <- "Subjects"
  rownames(AOV)[which(rownames(AOV) == "reps")] <- "Systematic ME"
  rownames(AOV)[which(rownames(AOV) == "groups")] <- "Groups"
  rownames(AOV)[which(rownames(AOV) == "reps:groups")] <- "Systematic ME:Groups"
  rownames(AOV)[which(rownames(AOV) == "Residuals")] <- "Random ME"
  
  if(multivariate) {
    
    rownames(mAOV)[which(rownames(mAOV) == "subj")] <- "Subjects / Random ME"
    rownames(mAOV)[which(rownames(mAOV) == "reps")] <- "Systematic ME / Random ME"
    rownames(mAOV)[which(rownames(mAOV) == "groups")] <- "Groups / Random ME"
    rownames(mAOV)[which(rownames(mAOV) == "reps:groups")] <- "Systematic ME:Groups / Random ME"
    mAOV <- mAOV[, -(1:2)]
    
  }
  
  AOV$EtaSq.ME <- NA
 
  Xsub <- model.matrix(~subj)
  R <- LM.fit(Xsub, Y)$residuals
  SS.tot <- sum(R^2)
  
  AOV$EtaSq.ME[which(rownames(AOV) == "Systematic ME")] <- 
    AOV$SS[which(rownames(AOV) == "Systematic ME")] / SS.tot
  AOV$EtaSq.ME[which(rownames(AOV) == "Systematic ME:Groups")] <- 
    AOV$SS[which(rownames(AOV) == "Systematic ME:Groups")] / SS.tot
  AOV$EtaSq.ME[which(rownames(AOV) == "Random ME")] <- 
    AOV$SS[which(rownames(AOV) == "Random ME")] / SS.tot
    
  indx <- length(na.omit(AOV$F))    
  AOV$F[1:indx] <- AOV$F[1:indx] * AOV$Df[1:indx] / AOV$Df[indx + 1]
  names(AOV)[which(names(AOV) == "F")] <- "SNR"
  AOV <- AOV[c(1:4, 6, 7, 8, 5)]
  
  # SSCP products
  
  S <- summary(fit)
  
  
  # for relative EVs
  if(p > 1) {
    
    SSCP.ME.product <- fast.solve(S$SSCP$Residuals) %*% 
      S$SSCP[[2]]
    sscp.sqrt <- Cov.proj(S$SSCP$Residuals, symmetric = TRUE)
    SSCP.ME.product.std<- as.matrix(t(sscp.sqrt) %*% 
                                                S$SSCP[[2]] %*% sscp.sqrt)

    } else  SSCP.ME.product <- SSCP.ME.product.std <- NULL
  
  micc.eigs <- NULL
  if(!is.null(micc)){
    micc.eigs <- lapply(1:length(micc), function(j){
      x <- micc[[j]]
      if(!is.null(x)){
        res <- Re(eigen(as.matrix(x), only.values = TRUE)$values)
      } else res <- NULL
      res
    })
    names(micc.eigs) <- names(micc)
  } 
  
  options(warn = wrn)
  
  out <- list(AOV = AOV, mAOV = mAOV, 
              icc = icc, 
              mult.icc.eigs = micc.eigs,
              SSCP = S$SSCP,
              SSCP.ME.product = SSCP.ME.product,
              SSCP.ME.product.std = SSCP.ME.product.std)
  rm(S)
  out.fit = if(multivariate) fitm else fit
  fit <- fitm <- icc <- mult.icc.eigs <- SSCP <-
    SSCP.ME.product <- SSCP.ME.product.std <- NULL
  
  
  if(!verbose) {
    out.fit$ANOVA <- out.fit$MANOVA <- NULL
    out.fit$LM$QR <- NULL
    out.fit$Models <- NULL
    out.fit$Models$full <- NULL
    out.fit$PermInfo$perm.schedule <- NULL
  } 
  
  out$LM <- out.fit$LM
  out.fit$LM <- NULL
  out$ANOVA <- out.fit$ANOVA
  out.fit$ANOVA <- NULL
  if(!is.null(out.fit$MANOVA))
    out$MANOVA <- out.fit$MANOVA
  out.fit$MANOVA <- NULL
  out$PermInfo <- out.fit$PermInfo
  out.fit$PermInfo <- NULL
  out$Models <- out.fit$Models
  out.fit$Models <- NULL
  rm(out.fit)
  out$verbose = verbose
  class(out) <- "measurement.error"
  out
}
    
  
  