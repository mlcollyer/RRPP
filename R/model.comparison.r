#' Model Comparisons, in terms of the log-likelihood, covariance trace,
#' or Z-score.
#'
#' Function calculates either log-likelihoods or traces of covariance 
#' matrices for comparison with 
#' respect to parameter penalties, or calculates Z-scores from RRPP, which
#' can be profiled across a gradient (predictor).
#' 
#'
#' The function calculates either log-likelihoods or traces of (residual) 
#' covariance matrices, plus parameter 
#' penalties, to assist in comparative model evaluation or selection.  
#' Because high-dimensional data often
#' produce singular or ill-conditioned residual covariance matrices, 
#' this function does one of two things: 1) uses 
#' the trace of a covariance matrix rather than its determinant; or 2) 
#' provides a ridge-regularization (Warton, 2008)
#' of the covariance matrix, only if it is determined that it is 
#' ill-conditioned.  Regardless of implementation,
#' covariance matrices are projected into a principal component (PC) 
#' space of appropriate dimensions.
#' 
#' The parameter penalty is based on that proposed by Bedrick and Tsai 
#' (1994), equal to 2(pk + p(p + 1)/2), where 
#' p is the appropriate dimension (not number of variables) of the 
#' covariance matrix.  The parameter, k,
#'  is the rank of the model design matrix.
#'  
#' In the case that "logLik" is chosen for the argument, type, 
#' AIC scores are calculated.  These scores
#' may not perfectly match other packages or software that 
#' calculate AIC for multivariate data, if ridge regularization
#' was used (and if other packages require p = the number of 
#' data variables).  When choosing logLik as the type of comparison,
#' it might be a good idea to adjust the tolerance or number of data 
#' principal components.  The default (NULL) values will
#' use all data dimensions to calculate log-likelihoods, which 
#' might cause problems if the number of variables exceeds the number 
#' of observations (producing singular residual covariance matrices).  
#' However, one should not reduce data dimensions haphazardly,
#' as this can lead to poor estimates of log-likelihood.  Furthermore, 
#' using the tolerance argument could result in different
#' numbers of principal components used for each model to calculate 
#' log-likelihoods, which might be a concern for comparing models.  
#' If both tol and pc.no arguments are used, the solution will use 
#' the fewest PCs produced by either argument.  Because the trace
#' of a covariance matrix is not sensitive to matrix singularity, 
#' no PC adjustment is used for the cov.trace argument.
#' 
#' This function can also calculate Z-scores from RRPP on model log-likelihoods, 
#' which can be compared directly or profiled along a gradient (predictor).  
#' This might be useful for for comparing generalized least-squares (GLS) models, 
#' for example, along a gradient of a parameter used to scale the covariance 
#' matrix for GLS estimation. See Collyer et al. 2022 for an example of using
#' RRPP on log-likelihoods with different covariance matrices.
#' 
#' Users can construct their own tables 
#' from the results but this function does not attempt to 
#' summarize results, as interpreting results requires 
#' some arbitrary decisions.  The \code{\link{anova}} function 
#' explicitly tests multiple models and can be used for nested 
#' model comparisons.
#' 
#' Results can also be plotted using the generic \code{\link{plot}} 
#' function.
#' 
#' Caution: For models with GLS estimation, the number of 
#' parameters used to estimate the covariance matrix
#' is not taken into consideration.  A generalized information 
#' criterion is currently in development.
#' 
#' 
#' 
#' @param ... Any number of lm.rrpp class objects for model fits 
#' to be compared.
#' @param type An argument to choose between log-likelihood, 
#' covariance trace, or Z results.  If Z is chosen, Z-scores are calculated, the same log-likelihoods
#' are calculated as with the log-likelihood type, but also in every RRPP permutation,
#' as describe for the initial mode fits, with choice of null model (below).
#' @param predictor An optional vector that can be used to profile the results based on type
#' across a range of numerical values described by the predictor.  A spline will also be fit, which
#' will reveal estimated values of the predictor that yield maximum and minimum values of model
#' comparison metric.
#' @param tol If type = logLik or Z, tol is a tolerance value between 0 and 1, 
#' indicating the magnitude below which 
#' components should be omitted (if standard deviations of 
#' components are less than the eigenvalue of the first 
#' component times the tolerance), for calculating the log-likelihood.
#' @param pc.no If type = logLik or Z, an optional value to indicate the 
#' number of principal components (maximum rank) to use 
#' for calculating the log-likelihood.
#' @param gls.null A logical value indicating whether GLS estimation should be
#' used with the null (intercept) model, for calculating Z scores via RRPP of 
#' log-likelihoods.  This should be FALSE if comparing different GLS estimations
#' of covariance matrices.  It should be TRUE if comparing different model fits
#' with the same GLS-estimated covariance matrix.
#' 
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{model.comparison} is a data 
#' frame with either log-likelihoods
#' or covariance traces, plus parameter penalties.  AIC scores 
#' might be include, if applicable
#' @references Bedrick, E.J., and C.L. Tsai. 1994. Model selection 
#' for multivariate regression in small samples. 
#' Biometrics, 226-231.
#' @references Warton, D.I., 2008. Penalized normal likelihood and 
#' ridge regularization of correlation and covariance matrices. 
#' Journal of the American Statistical Association. 103: 340-349.
#' @references Collyer,  M.L., E.K. Baken, & D.C. Adams.  A standardized 
#' effect size for evaluating and comparing the strength of phylogenetic 
#' signal. Methods in Ecology and Evolution. 13: 367–382.
#' @examples 
#' \dontrun{
#' data(Pupfish)
#' Pupfish$logSize <- log(Pupfish$CS)
#' fit1 <- lm.rrpp(coords ~ logSize, data = Pupfish, iter = 0, 
#' print.progress = FALSE)
#' fit2 <- lm.rrpp(coords ~ Pop, data = Pupfish, iter = 0, 
#' print.progress = FALSE)
#' fit3 <- lm.rrpp(coords ~ Sex, data = Pupfish, iter = 0, 
#' print.progress = FALSE)
#' fit4 <- lm.rrpp(coords ~ logSize + Sex, data = Pupfish, iter = 0, 
#' print.progress = FALSE)
#' fit5 <- lm.rrpp(coords ~ logSize + Pop, data = Pupfish, iter = 0, 
#' print.progress = FALSE)
#' fit6 <- lm.rrpp(coords ~ logSize + Sex * Pop, data = Pupfish, iter = 0, 
#' print.progress = FALSE)
#' 
#' modComp1 <- model.comparison(fit1, fit2, fit3, fit4, fit5, 
#' fit6, type = "cov.trace")
#' modComp2 <- model.comparison(fit1, fit2, fit3, fit4, fit5, 
#' fit6, type = "logLik", tol = 0.01)
#' 
#' summary(modComp1)
#' summary(modComp2)
#' 
#' par(mfcol = c(1,2))
#' plot(modComp1)
#' plot(modComp2)
#' 
#' # Comparing fits with covariance matrices
#' # an example for scaling a phylogenetic covariance matrix with
#' # the scaling parameter, lambda
#' 
#' data("PlethMorph")
#' Cov <- PlethMorph$PhyCov
#' lambda <- seq(0, 1, 0.1)
#' 
#' Cov1 <- scaleCov(Cov, scale. = lambda[1])
#' Cov2 <- scaleCov(Cov, scale. = lambda[2])
#' Cov3 <- scaleCov(Cov, scale. = lambda[3])
#' Cov4 <- scaleCov(Cov, scale. = lambda[4])
#' Cov5 <- scaleCov(Cov, scale. = lambda[5])
#' Cov6 <- scaleCov(Cov, scale. = lambda[6])
#' Cov7 <- scaleCov(Cov, scale. = lambda[7])
#' Cov8 <- scaleCov(Cov, scale. = lambda[8])
#' Cov9 <- scaleCov(Cov, scale. = lambda[9])
#' Cov10 <- scaleCov(Cov, scale. = lambda[10])
#' Cov11 <- scaleCov(Cov, scale. = lambda[11])
#' 
#' 
#' fit1 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov1)
#' fit2 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov2)
#' fit3 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov3)
#' fit4 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov4)
#' fit5 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov5)
#' fit6 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov6)
#' fit7 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov7)
#' fit8 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov8)
#' fit9 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov9)
#' fit10 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov10)
#' fit11 <- lm.rrpp(SVL ~ 1, data = PlethMorph, Cov = Cov11)
#' 
#' par(mfrow = c(1,1))
#' 
#' MC1 <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6,
#' fit7, fit8, fit9, fit10, fit11,
#' type = "logLik")
#' MC1
#' plot(MC1)
#' 
#' MC2 <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6,
#' fit7, fit8, fit9, fit10, fit11,
#' type = "logLik", predictor = lambda)
#' MC2
#' plot(MC2)
#' 
#' 
#' MC3 <- model.comparison(fit1, fit2, fit3, fit4, fit5, fit6,
#' fit7, fit8, fit9, fit10, fit11,
#' type = "Z", predictor = lambda)
#' MC3
#' plot(MC3)
#' }

model.comparison<- function(..., type = c("cov.trace", "logLik", "Z"), 
                            predictor = NULL, tol = NULL, pc.no = NULL,
                            gls.null = FALSE) {
  
  dots <- list(...)
  check <- unlist(lapply(dots, inherits, "lm.rrpp"))
  if(any(!check)) stop("\nObjects must be lm.rrpp fits\n.")
  dot.names <- lapply(dots, function(x) x$LM$Terms[[3]])
  if(length(unique(dot.names)) < length(dot.names)) {
    cat("Models do not have unique term combinations.\n")
    cat("Labelling models based on order presented.\n")
    dot.names <- paste("Mod", 1:length(dots), sep = ".")
  }
  type = match.arg(type)
  
  if(type == "logLik") {
    if(!is.null(tol)) {
      if(tol > 1 || tol < 0) 
        stop("tol must be a value between 0 and 1\n", call. = FALSE)
    }
    if(!is.null(pc.no)) {
      p.list <- sapply(dots, function(x) x$LM$p.prime)
      rank.list <- ifelse(p.list < pc.no, p.list, pc.no)
      rank.list <- as.data.frame(matrix(rank.list, nrow = 1))
      rownames(rank.list) <- "pc.no"
      colnames(rank.list) <- dot.names
      
      if(any(pc.no > p.list)) {
        
        warning(
          paste(
            "\nThis is not an error!  It is a friendly warning.\n",
            "\nThe number of PCs requested exceeds possible dimensions for some model fits.",
            "\nPC number will be dropped initially to", rank.list, 
            "and might be reduced further, if the tolerance was also adjusted.\n",
            "\nUse options(warn = -1) to turn off these warnings. \n\n", sep = " "),
          noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
        
        pc.no <- p.list
      }
    }
    
    ll.args <- list(fit = dots[[1]], tol = tol, pc.no = pc.no)
    temp <- sapply(1:length(dots), function(j){
      ll.args$fit <- dots[[j]]
      do.call(logL, ll.args)
    }) 
    res <- unlist(temp[1,])
    rank <- unlist(temp[2,])
  } else res <- sapply(dots, cov.trace) 
  
  if(type == "Z") {
    res <- sapply(1:length(dots), function(j){
      f <- dots[[j]]
      z <- logLik(f, Z = TRUE, tol = tol, 
                  pc.no = pc.no, gls.null = gls.null)$Z
      if(is.na(z)) z <- 0
      z      
    })
  }
  
  par.pen <- function(f){
    p <- f$LM$p.prime
    k <- f$LM$QR$rank
    2 * (p * k + 0.5 * p* (p + 1))
  }
  pp <- sapply(dots, par.pen)
  
  if(type == "logLik") {
    out <- as.data.frame(cbind(res, rank, pp, -2*res + pp))
    names(out) <- c("logLik", "residual.pc.no", "penalty", "AIC")
  } else if(type == "cov.trace") {
    out <- as.data.frame(cbind(res, pp))
    names(out) <- c("cov.trace", "penalty")
  } else {
    out <- data.frame(Z = res)
  }
  
  if(!is.null(predictor)) {
    out$predictor <- predictor
    names(out)[which(names(out) == "predictor")] <- deparse(substitute(predictor))
  }
   
  rownames(out) <- unlist(dot.names)
  out <- list(table = out, names = dot.names)
  attr(out, "class") = "model.comparison"
  out
  
}
