#' Intraclass correlation statistics from an lm.rrpp model fits
#'
#' Function performs analyses concerned with the repeatability (reliability) of multivariate data 
#' (measurements) collected from the same research subjects.  Although there is no
#' requirement for repeated measurements on all research subjects, the analysis assumes
#' that multiple observations are made. 
#' 
#' Function uses ANOVA statistics or SSCP matrices to find the ratio of among-subject
#' to within-subject variance.  The former is a dispersion-based approach and the latter
#' is a multivariate generalization of the ICC statistic (as a matrix product).  The
#' multivariate generalizations of the statistics
#' described by Liljequist et al. (2019) are used to find matrix products, 
#' from which eigenanalysis is performed, providing ICC statistics by eigenvectors.  
#' 
#' Three statistics describe the ICC for the population,
#' agreement of measurements among subjects, and consistency between measurements.  
#' The last statistic does not necessarily measure the sameness 
#' between measurements but the consistency of change between measurements,
#' which might be indicative of a systematic measurement error.  
#' If groups are used, these three statistics are 
#' repeated, using the SSCP for groups-adjusted data.  
#' This approach accounts for group differences,
#' which would avoid large subject variation compared to measurement error 
#' inflating ICC values.  If there are 
#' inherently disparate groups from which subjects are sampled, 
#' this approach can elucidate better agreement and 
#' consistency in light of group differences.
#' 
#' This function is most useful for analyses performed with 
#' \code{\link{measurement.error}}, but any \code{\link{lm.rrpp}} fit can be used,
#' so long as research subjects can be defined. 
#' 
#' It is essential that all arguments are terms that can be found in the model frame
#' of the model fit, as provoke by ANOVA.  Using anova(fit) will elucidate the row
#' names of the ANOVA that could be used.
#' 
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @return Objects of class "ICCstats" return the following:
#'  \item{ICC_disp}{The intraclass correlation coefficient (ICC) based on the dispersion of values.}
#'  \item{ICC_mult}{The eigenvalues of ICC matrices} 
#' @param fit The \code{\link{lm.rrpp}}, previously evaluated.
#' @param subjects A single character value indicating which term in an ANOVA table
#' corresponds to research subjects.
#' @param with_in One or more character values indicating which terms in an ANOVA table
#' are measured within subjects (replications, plus maybe interactions).  If NULL,
#' the only replication within-subject will be considered as residuals.
#' @param groups An optional character value to indicate if a factor in the 
#' model frame of the \code{\link{lm.rrpp}} fit that could account for subject variation.
#' Using this argument might minimize the importance of subject variation, if subjects
#' have disparate values that could inflate ICC.  Note that this name could be different
#' than what is shown in the ANOVA table, if \code{\link{measurement.error}}
#' was used.  Use names(fit$LM$data), substituting fit with the name assigned
#' to the \code{\link{measurement.error}} object, to know the groups factor, 
#' if used.
#' @param multivariate Logical value for whether to include to calculate ICC matrix 
#' generalizations and perform eigenanalysis.
#' @param print.AOV Logical value for whether to include ANOVA table as screen 
#' output, when calculating ISS statistics.
#'
#' Note that this function can return ICC statistics, even if they do not make sense.
#' It is possible to generate ICC stats with any ANOVA table, with at least one
#' term.  
#'
#' @references Liljequist, D., Elfving, B., & Skavberg Roaldsen, K. (2019). Intraclass correlation–A discussion 
#' and demonstration of basic features. PloS one, 14(7), e0219854.
#' 
#' @examples
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
#' ICCstats(ME2, subjects = "Subjects", 
#'   with_in = c("Systematic ME", "Systematic ME:Groups"), 
#'   groups = "groups")

ICCstats <- function(fit, 
                     subjects = NULL,
                     with_in = NULL,
                     groups = NULL,
                     multivariate = FALSE,
                     print.AOV = TRUE){
  if(!inherits(fit, "lm.rrpp"))
    stop(paste("\nThe fit argument must use an object of",
    " class lm.rrpp.", sep = ""), call. = FALSE)
  
  AOV <- anova(fit)$table
  if(print.AOV) print(AOV)
  
  if(is.null(subjects) || all(!rownames(AOV) %in% as.character(subjects)))
    stop(paste("\n\nThe ANOVA (if printed above) should be able to help you choose the",
               "correct subjects name (row name in table)\n", sep = " "), 
         call. = FALSE)
  
  if(!is.numeric(subjects)) 
    sub_no = which(rownames(AOV) %in% as.character(subjects))
  
  if(is.null(with_in)) with_in <- "Residuals"
  
  if(!is.numeric(with_in)) 
    with_in_no = which(rownames(AOV) %in% as.character(with_in))
  if(length(with_in_no) == 0) {
    with_in <- "Random ME"
    with_in_no = which(rownames(AOV) %in% with_in)
  }
  
  
  if(!is.null(groups)){
    gfac <- unlist(fit$LM$data[names(fit$LM$data) %in% groups])
    if(is.null(gfac))
      stop("\nThe groups variable is not found in the fit model frame.\n",
            call. = FALSE)
    Y <- fit$LM$Y
    df <- rrpp.data.frame(Y = Y, gfac = gfac)
    
    R <- resid(lm.rrpp(Y ~ gfac, 
               Cov = fit$LM$Cov, weights = fit$LM$weights,
               iter = 1, data = df,
               print.progress =  FALSE))
    form <- formula(fit$LM$Terms)
    fitb <- fit
    fitb$LM$Y <- fitb$LM$data$Y <- R
    if(fit$ANOVA$SS.type == "Within-subject II") {
      fit <- suppressWarnings(
        lm.rrpp(form, 
                Cov = fitb$LM$Cov, weights = fitb$LM$weights,
                iter = 1, data = fitb$LM$data, 
                print.progress = FALSE, 
                subjects = fit$LM$data[names(fit$LM$data) == fit$subjects.var]))
    } else {
      fit <- suppressWarnings(lm.rrpp(form,  
                                      Cov = fitb$LM$Cov, 
                                      weights = fitb$LM$weights,
                                      iter = 1, data = fitb$LM$data, 
                                      SS.type = fit$ANOVA$SS.type,
                                        print.progress =  FALSE))
    }
    
    rm(fitb)
    AOV2 <- anova(fit)$table
    rownames(AOV2) <- rownames(AOV)
    AOV <- AOV2
    rm(AOV2)
  }
  
  resid_no<- which(rownames(AOV) == "Residuals")
  if(length(resid_no) == 0)
    resid_no<- which(rownames(AOV) == "Random ME")
  
  type <- 1
  if(!identical(with_in_no, resid_no)) type <- 2
  
  sub.line <- AOV[sub_no,]
  w.lines <- AOV[unique(c(with_in_no, resid_no)),]
  res.line <- AOV[resid_no,]
  
  MSBS <- sub.line[3]
  if(type == 2)
    MSWS <- sum(w.lines[,2]) / sum(w.lines[,1]) else
      MSWS <- res.line[3]
  
  if(type == 2){
    MSBM <- w.lines[1, 3]
    MSE <- res.line[3]
  }
  ns <- sub.line[1] + 1
  nr <- w.lines[1, 1] + 1
  
  icc <- (MSBS - MSWS) / (MSBS + (nr - 1) * MSWS)
  if(type == 2){
    icca <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE + 
                              nr / ns * (MSBM - MSE))
    iccc <- (MSBS - MSE) / (MSBS + (nr - 1) * MSE)
  } else icca <- iccc <- NULL
  
  icc.out <- unlist(c(icc , icca, iccc))
  names(icc.out) <- c("ICC", "ICC_a", "ICC_c")[1:length(icc.out)]
  
  if(multivariate){
    
    fitb <- fit
    class(fitb) <- "lm.rrpp"
    SSCP <- summary(fitb, formula = FALSE)$SSCP
    rm(fitb)
    
    Df <- fit$ANOVA$df
    rmove <- which(Df == 0)
    if(length(rmove) > 0)
      SSCP <- SSCP[-rmove]
    
    if(!identical(names(SSCP), rownames(AOV))){
      if(length(SSCP) != (nrow(AOV) - 1))
        stop("\nIncommensurate ANOVA and MANOVA terms\n",
             call. = FALSE)
      names(SSCP) <- rownames(AOV)[1:length(SSCP)]
    }
    
    sub_no <- which(names(SSCP) %in% subjects)
    with_in_no <- which(names(SSCP) %in% with_in)
    resid_no <- length(SSCP)

    RSSCP <- as.matrix(SSCP[[resid_no]])
    SSCPBS <- as.matrix(SSCP[[sub_no]])
    if(type == 2){
      SSCPBM <- lapply(with_in_no, function(x) SSCP[[x]])
      SSCPBM <- Reduce("+", SSCPBM)
      SSCPWS <- RSSCP + SSCPBM
    } else SSCPWS <- RSSCP
    
    Df <- AOV$Df
    
    MSCPBS <- SSCPBS / Df[sub_no]
    
    if(type == 2) {
      MSCPBM <- SSCPBM / sum(Df[with_in_no])
      MSCPWS <- SSCPWS / sum(Df[c(with_in_no, resid_no)])
      MSCPE <- RSSCP / Df[resid_no]
    } else MSCPWS <- SSCPWS / Df[resid_no]
  
    micc <- fast.solve(MSCPBS + (nr - 1) * MSCPWS) %*% (MSCPBS - MSCPWS) 
    
    if(type == 2){
      micca <- fast.solve(MSCPBS + ((nr - 1) * MSCPE) + 
                            as.numeric(nr / ns) * 
                            (MSCPBM - MSCPE)) %*% (MSCPBS - MSCPE) 
      miccc <- solve(MSCPBS + (nr - 1) * MSCPE) %*% (MSCPBS - MSCPE) 
    } else micca <- miccc <- NULL
    
    
    ev_icc <- eigen(micc)$values
    ev_icc <- Re(ev_icc)
    ev_icc <- ev_icc[ev_icc > 0]
    
    if(type == 2){
      ev_icc_a <-eigen(micca)$values
      ev_icc_a <- Re(ev_icc_a)
      ev_icc_a <- ev_icc_a[ev_icc_a > 0]
      ev_icc_c <-eigen(miccc)$values
      ev_icc_c <- Re(ev_icc_c)
      ev_icc_c <- ev_icc_c[ev_icc_c > 0]
    } else ev_icc_a <- ev_icc_c <- NA
    
    k <- max(c(length(ev_icc), 
               length(ev_icc_a), 
               length(ev_icc_c)), na.rm = TRUE)
    nrm <- if(type == 1) 1 else 3
    iccm.out <- matrix(NA, nrm, k)
    
    iccm.out[1, 1:length(ev_icc)] <- ev_icc
    if(type == 2){
      iccm.out[2, 1:length(ev_icc_a)] <- ev_icc_a
      iccm.out[3, 1:length(ev_icc_c)] <- ev_icc_c
    }
    rownames(iccm.out) <- names(icc.out)
    colnames(iccm.out) <- paste("EV", 1:k, sep = "")
  } else iccm.out <- NULL
  
  out <- list(ICC_disp = icc.out, ICC_mult = iccm.out)
  class(out) <- "ICCstats"
  out
}
  