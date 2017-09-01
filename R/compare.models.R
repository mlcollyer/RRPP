#' Compare Linear Models via Randomized Residual Permutation Procedure
#'
#' Function performs various "model selection" procedures and a non-parametric likelihood
#' ratio test for separate models previously evaluated with \code{lm.rrpp}
#'
#' The function allows input of an unlimited number of \code{\link{lm.rrpp}} objects, to evaluate 
#' model performance in a comaprative manner.  Emphasis is placed on "model performance", 
#' as some of the analyses performed use RRPP to generate residual sum of squares (RSS) for model 
#' evaluation, rather than using parametric model selection.  
#' 
#' There are four important a priori distinctions that must be made to understand how this function operates.
#'
#' 1. Although the function is non-parametric in nature, a parametric model section criterion, 
#' Akaike's Information Criteion (AIC) is also provided for descriptive purposes.  The method of 
#' Bedrick and Tsai (1994) is used to estimate AIC.  A small-sample correction is not used, as it 
#' would preclude use of high-dimensional data.  This limitation is more reason to use a non-parametric 
#' method in such cases.  Additionally, error-covariance matrices are projected into a subspace of 
#' appropriate dimensionality before calculating a determinant.  This subspace is defined in terms of 
#' dimensions by the eigen values of a singular value decomposition of the covariance matrix.  To assure 
#' the result is a positive-definite matrix, the subspace is limited to dimensions associated with eigenvalues 
#' comprising 99.9 percent of the sum of all eigenvalues.  This process allows the estimation of the 
#' determinant for calculating the log-likelihood of the model in order to calculate AIC.
#'
#' 2. Adjusted multivariate coefficients of determination (aRsq) provide a non-parametric alternative 
#' to AIC that mimics AIC in purpose.  Like AIC, the standard coefficient of determinaton (Rsq) is 
#' penalized by the number of parameters used to estimate model sum of squares (SS), limiting the selection of overfitted
#' models when comparing among multiple candidate models.  Unlike AIC, aRsq is not influenced by the
#' number of data dimensions.  (SS and RSS are calculated as the trace of model and model residual sums of
#' squares and cross-products matrices, respectively, which is the same as the sum of SS over all variables.
#' Because Rsq is a ratio, the number of dimensions does not directly influence the scale of this statistic.
#' However, it is possible that the number and types of variables used can influence the upper limit of Rsq
#' and aRsq values.  This property is of less concern though for comparative procedures, as it would apply
#' to all models compared.)  Also, unlike AIC, a higher value rather than a lower value suggests a 
#' "better" model, in comparisons.  Finally, aRsq as a statistic can be appreciated as the amount of 
#' variation explained by the model, adjusted for the number of parameters.  In this sense, it has more
#' absolute meaning than AIC.

#'
#' 3. Irrespective of parameter-adjustement, aRsq as a sole statistic will have a tendency to 
#' favor over-fitted models comapred to AIC.  However, it is not the observed value of aRsq but the 
#' distribution of aRsq over many RRPP permutations that is of interest.  The distribution of aRsq values
#' across RRPP permutations is approximately normally distributed (Collyer and Adams, insipient).
#' This function provides the mean aRsq value from RRPP permutations, delta aRsq values after subtracting 
#' all aRsq values from the maximum value, plus confidence intervals for delta aRsq values.  Any 
#' confidence interval bracketing 0 might suggest a viable model.  One can choose the simplest among these 
#' as a preferred model.  One can also choose the model with the smallest confidence range as a preferred model, 
#' as it most consistently produces favorable aRsq results.  This approach has some similarity to 
#' Deviance Information Criteria (DIC; see Spigelhalter et al. 2002; Gelman et al. 2004; Rafferty et al. 2007).  
#' DIC evaluates model performance (by measuring Deviance) rather than penalizing log-likelihood over many MCMC
#' permutations.  RRPP is analagous to MCMC and Rsq is the ratio of model SS to model deviance for continuous
#' variables with normally distributed residuals, for models fit with least squares.  However, deviance does not
#' have an applicable multivairate generalization (RSS could be used, but would scale by the number of variables)
#' and issues with DIC favoring overfitted models are known (because adding superfluous model parameters can only 
#' improve model performance).  The approach used here with aRsq captures the spirit of DIC but does not employ
#' some of the deficiencies DIC would impose with highly multivariate data.
#' 
#' 4. The difference in model RSS across RRPP permutations provides a distribution of test-statistics, from which 
#' a non-parametric analog of a likelihood ratio test can be performed.  For these tests, the probability of these 
#' statistics (SS between models, which is the difference in RSS between models) can be inferred from the many 
#' outcomes generated by RRPP.  By standardizing this distribution, the observed SS is also a standard deviate 
#' (effect size) for the difference between a reference and target model.  If large, it suggests the target 
#' model is better than the reference model for explaining the data (and the probability of finding such large SS 
#' by chance can be used as a p-value for the null hypthesis that the two models have the same deviance).  
#' A test like this normally expects the reference model 
#' to be nested within the target model, though using non-nested models in this permutational strategy does not preclude 
#' generation of results.
#'
#' This function requires input as: model fit 1, model fit 2, model fit 3, ...; the first model will be considered a 
#' reference model for hypothesis testing or just one of the models for AIC or aRsq comparisons.
#'
#' \subsection{Additional Computational Notes and Suggestions}{ 
#' *** The method of calculating AIC weights is generally exp(-delta AICi/2) / sum(exp(-delta AICi/2)), where i refers 
#' to the ith model in consideration and delta means the difference between model i AIC and the minimum AIC for all 
#' candidate models.  This method assumes a univariate dependent variable, and weights calculated this way are 
#' dimension-dependent. Rather, calculating AIC weights as exp(-delta AICi/(2p)) / sum(exp(-delta AICi/(2p))), where 
#' p is the rank of the data space (generally, the number of variables, but could be less if there are variable redundancies),
#' generalizes the univariate method (which remains unchanged for univariate data) and removes the association between 
#' weighting and data dimensionality.  Thus, interpretations of weights can be made consistently across different data sets.
#' 
#' *** Likewise, the AIC parameter penalty is calculated as (2 times) p*k + 0.5p(p+1), where k is the rank of the model design 
#' matrix (typically the number of model parameters), and p is the rank of the data space.  The first part of this equation is 
#' the number of coefficients for the Beta matrix; the second part is the number of unique elements in the error covariance matrix.  
#' If p = 1, this equation reduces to k + 1, where k contains parameters for a mean (intecerpt), slopes, and the + 1 is the for the 
#' error variance.  Adding a noise variable to p increases to 2k + 3; adding another increases to 3k + 6; adding another increases 
#' to 4k + 10; and so on, illustrating the impact highly multivariate data can have on AIC weights estimated in a univariate context.
#' 
#' *** These dimensionality issues are not found with confidence intervals for aRsq values.
#' 
#' *** Results will include "Deviance", which is the sum of variable deviances (RSS) across all dependent
#' variables.  This is not intended to introduce a multivariate continuous variable definition of deviance, but rather
#' offer results analagous to univariate model selection tables, for descriptive purpose.
#' 
#' *** This function will allow one to compare covariance structures.  Caution should be used in
#' interpreting AIC scores in these cases.  The number of parameters used to estimate covariance matrices
#' will have large influence on AIC scores.  (An argument for the number of covariance parameters can be
#' provided in \code{\link{lm.rrpp}} in the cases that correaltion structure is desired.)  The non-parametric 
#' aRsq results will be more valuable in these cases, and should be preferred to AIC score interpretation.  
#' Additionally, one should realize that the model comparison results indicate which models consistently 
#' produce the least error, but when comparing covariance structure, this is not necessarily key to choosing the most
#' appropriate model.  For example, a covariance matrix for phylogenetic relatedness might yield a model
#' that performs poorly comapred to one with no covariance matrix, but this does not mean the phylogenetic relatedness
#' should be ignored.  (Rather it probably signals the non-independence of observations if ignored.)
#' }
#' 
#' @param refModel A class \code{lm.rrpp} object with more than one random permutation.
#' @param ... Any additional class \code{lm.rrpp} objects intended for comparison.  There are two methods of
#' comparison: information criteria (AIC, aRsq) or non-parametric testing of model RSS.  The former does
#' not require that the \code{refModel} be chose carefully, but for direct testing, the refModel will be the sole model
#' for hypothesis test statistics.  It is assumed that the same seed was used in \code{\link{lm.rrpp}} for the random permutations,
#' and that the same number of iterations were used.  Failure to ensure this prerequisite could result in spurious outcomes.
#' See \code{\link{lm.rrpp}} for further details.
#' @param confidence The level for confidence intervals on Delta-aRsq measures.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.
#' This is helpful for long-running analyses.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{compare.models} is a list containing the following
#' \item{AIC.table}{A table of model comparison statistics using AIC.}
#' \item{aRsq.table}{A table of model comparison statistics using aRsq}
#' \item{LRT.table}{A table for LRT test stats for models compared to the reference model.  logLR is
#' calculated as twice the difference in log-likelihoods between null and comapred models.  Note that 
#' statistics will be generated, even if the comparisons are not nested or do not make much sense.}
#' \item{logL}{Model by model log-likelihoods.}
#' \item{deviances}{Model by model random deviances, based on RRPP iterations described by the model.}
#' \item{aRsq}{Model by model random aRsq values, based on RRPP iterations described by the model.}
#' \item{LRT.SS}{The sum of squares for every model compared to the reference model, as a statistic for a 
#' likelihood ratio test.  These values are differences in deviances for all RRPP permutations.
#' If a complex model is chosen as the reference model comapre to subsequent models, SS can be negative, 
#' and the test will not make much sense.}
#'
#' @references Bedrick, E.J., and C.L. Tsai.  1994. Model selection for multivariate regression in small samples.
#' Biometrics 50:226–231.
#' @references Collyer, M.L. and D.C. Adams. 2018?  Isipient world class manuscript.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., Rubin, D.B. 2004. Bayesian Data Analysis: Second Edition.
#' Texts in Statistical Science, CRC Press.
#' @references Raftery, A.E., M.A. Newton, J.M. Satagopan, and P.N. Krivitsky. 2006. Estimating 
#' the integrated likelihood via posterior simulation using the harmonic mean identity.
#' In. Bernardo et al. (eds) Bayesian Statistics. Oxford University Press.
#' @references Spiegelhalter, D.J.,  Best, N.G., Carlin, B.P., van der Linde, A. 2002. Bayesian measures of model complexity and fit
#' (with discussion). Journal of the Royal Statistical Society, Series B. 64:583–639.
#' @examples
#' ### Example with pupfish shape data, to understand if population source, sex, and/or
#' # fish size are important sources of variation.  Variable interactions are also considered.
#' 
#' data(pupfish)
#' # cs = centroid size, sex, population are possible independent variables;
#' # coords are Procrustes residuals from a generalized Procrustes analysis (GPA)
#'
#' rdf <- rrpp.data.frame(Y = pupfish$coords, cs = pupfish$CS, 
#' pop = pupfish$Pop, sex = pupfish$Sex)
#' fit0 <- lm.rrpp(Y ~ 1, data = rdf)
#' fit1 <- lm.rrpp(Y ~ cs, data = rdf)
#' fit2 <- lm.rrpp(Y ~ sex, data = rdf)
#' fit3 <- lm.rrpp(Y ~ pop, data = rdf)
#' fit4 <- lm.rrpp(Y ~ cs * sex, data = rdf)
#' fit5 <- lm.rrpp(Y ~ cs * pop, data = rdf)
#' fit6 <- lm.rrpp(Y ~ cs + sex * pop, data = rdf)
#' fit7 <- lm.rrpp(Y ~ cs * sex * pop, data = rdf)
#'
#' CM <- compare.models(fit0, fit1, fit2, fit3, 
#' fit4, fit5, fit6, fit7)
#' summary(CM)
#' 
#' # Different outcomes based on AIC and aRsq
#' # try again to perform a likelihood ratio test between relevant models
#' 
#' CM2 <- compare.models(fit5, fit6, fit7)
#' summary(CM2)
#' 
#' # Improvement over fit 5 is "not significant" (at alpha = 0.05), 
#' # even if fit 6 and fit 7 perform better.  

# NEED BETTER EXAMPLES

compare.models <- function(refModel, ..., confidence = 0.95, print.progress = TRUE){
  dots <- list(...)
  if(length(dots) == 0) stop("This function requires multiple lm.rrpp objects")

  classes <- unlist(Map(function(j) inherits(j, "lm.rrpp"), dots))

  if(any(!classes)) stop("Not all objects are class lm.rrpp.")

  if(!inherits(refModel, "lm.rrpp")) stop("The reference model is not class lm.rrpp")

  dots <- c(list(refModel), dots)
  nlist <- sapply(dots, function(j) j$LM$n)
  plist <- sapply(dots, function(j) j$LM$p)

  if(length(unique(nlist)) > 1) stop("Unequal numbers of observations among models.")

  if(length(unique(plist)) > 1) stop("Unequal numbers of dependent variables among models.")

  n <- unique(nlist)
  p <- unique(plist)
  k <- length(dots)
  ranks <- sapply(dots, function(j) j$LM$QR$rank)
  rank.check <- sapply(1:k, function(j) is.null(ranks[[j]]))
  ranks[rank.check == TRUE] <- 1
  if(print.progress){
    cat("\n\nCalculating model deviances for the number of 
        lm.rrpp permutations for", k, "models\n")
    pb <- txtProgressBar(min = 0, max = k+1, initial = 0, style=3)
  }
  deviances <- lapply(1:k, function(j){
    x <- dots[[j]]
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    if(!is.null(x$LM$Cov)) DevianceGLS(x) else
      DevianceOLS(x)
  })
  
  step <- k+1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  
  Dev <- lapply(1:k, function(j) deviances[[j]]$Deviance)
  Dev.obs <- unlist(lapply(1:k, function(j) (Dev[[j]])[[1]]))
  
  # AIC
  AIC <- sapply(1:k, function(j) deviances[[j]]$AIC)
  logL <-sapply(1:k, function(j) deviances[[j]]$logL)
  delAIC <- AIC - min(AIC)
  AICw <- round(exp(-0.5*delAIC/p)/sum(exp(-0.5*delAIC/p)), 4)
  atab <- data.frame(Deviance = Dev.obs,
                    logL = logL, AIC = AIC, 
                    delAIC = delAIC, AICw = AICw)
  dotnms <- substitute(list(...))[-1]
  dotnms <- c(deparse(substitute(refModel)), sapply(substitute(dotnms), deparse))
  colnames(atab)[1] <- "Deviance (RSS)"
  rownames(atab) <- names(Dev) <- dotnms
  rownames(atab)[1] <- paste(rownames(atab)[1], "(ref)")
  
  # eta-sq
  Rsq <- sapply(1:k, function(j) deviances[[j]]$Eta.sq.a)
  Rsq.obs <- Rsq[1,]
  Rsq.mean <-colMeans(Rsq)
  etamax <- max(Rsq.mean)
  delRsq <- etamax - Rsq.mean
  Rsq.c <- etamax - Rsq
  eLCL <- apply(Rsq.c, 2, function(x) quantile(x, (1-confidence)/2))
  eUCL <- apply(Rsq.c, 2, function(x) quantile(x, confidence + (1-confidence)/2))
  Rsq.mean <- zapsmall(Rsq.mean)
  etab <- data.frame(Deviance = Dev.obs,
                     aRsq.mean = Rsq.mean, delRsq = delRsq, 
                     LCL = eLCL, UCL = eUCL, CI.Range = eUCL - eLCL)
  
  rownames(etab) <- rownames(atab)
  colnames(etab)[1] <- "Deviance (RSS)"
  colnames(etab)[4] <- paste((1-confidence)/2*100, "% LC", sep="")
  colnames(etab)[5] <- paste((confidence + (1-confidence)/2)*100, "% UC", sep="")

  # LRT
  ind.list <- lapply(1:k, function(j) dots[[j]]$PermInfo$perm.schedule)
  ind.check <- sapply(2:k, function(j) identical(ind.list[[j]], ind.list[[j - 1]]))
  ind <- ind.list[[1]]
  perms <- length(ind)
  if(any(!ind.check)) {
    cat("\nWarning: Unequal permutation schedules, so an LRT cannot be perfomed\n")
    ltab <- NULL
  } else {
    
    if(print.progress){
      cat("\n\nPerforming", k-1, "non-parametric likelihood ratio tests, each with",
          perms, "permutations\n")
      pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
    }
    
    refFit <- refit(refModel)
    refFit.set <- remodel.set(refFit)
    if(refModel$LM$gls) {
      P <- refModel$LM$Pcov
      X <- refModel$LM$X * sqrt(refModel$LM$weights)
      refFit.set$U <- P%*%qr.Q(qr(crossprod(P, X)))
    }
    
    SS <- lapply(1:perms, function(j){
      refFit.set$ind.i <- ind[[j]]
      y <- do.call(remodel.Y, refFit.set)
      SSn <- SS.mean(y, n)
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      SSm <- sapply(2:k, function(jj){
        f <- dots[[jj]]
        if(f$LM$gls){
          x <- f$LM$X * sqrt(f$LM$weights)
          P <- f$LM$Pcov
          U <- P%*%qr.Q(qr(crossprod(P, x)))
        } else U <- qr.Q(f$LM$wQR)
        sum(crossprod(U, y)^2) - SSn
      })
    })
    step <- perms + 1
    if(print.progress) {
      setTxtProgressBar(pb,step)
      close(pb)
    }
    SS <- cbind(1, t(matrix(unlist(SS), k-1, perms)))
    pvals <- apply(SS, 2, pval)
    Z <- apply(log(SS), 2, effect.size)
    logLR <- 2*(logL - rep(logL[1],k))
    ltab <- data.frame(Deviance = Dev.obs,
                       logLR = logLR, SS = SS[1,], 
                       Z = Z, P = pvals)
    ltab[1,-1] <- NA
    colnames(ltab)[NCOL(ltab)] <- "Pr(>SS)"
    colnames(ltab)[1] <- "Deviance (RSS)"
    rownames(ltab) <- colnames(SS) <- rownames(atab)
    class(ltab) <- c("anova", class(ltab))
  }
  out <- list(AIC.table = atab, 
              aRsq.table = etab,
              LRT.table = ltab, 
              logL = logL, Deviance = Dev, 
              aRsq = Rsq, LRT.SS = SS[,-1])
  class(out) <- "compare.models"
  out
}

