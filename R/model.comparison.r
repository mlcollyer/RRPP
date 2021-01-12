#' Model Comparisons, in terms of the log-likelihood or covariance trace
#'
#' Function calculates either log-likelihoods or traces of covariance 
#' matrices for comparison with 
#' respect to parameter penalties.
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
#' @param ... Any number of lm.rrpp class objects for model fits 
#' to be compared.
#' @param type An argument to choose between log-likelihood or 
#' covariance trace results
#' @param tol If type = logLik, tol is a tolerance value between 0 and 1, 
#' indicating the magnitude below which 
#' components should be omitted (if standard deviations of 
#' components are less than the eigenvalue of the first 
#' component times the tolerance), for calculating the log-likelihood.
#' @param pc.no If type = logLik, an optional value to indicate the 
#' number of principal components (maximum rank) to use 
#' for calculating the log-likelihood.
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
#' @examples 
#' 
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
model.comparison<- function(..., type = c("cov.trace", "logLik"), 
                            tol = NULL, pc.no = NULL) {
  
  dots <- list(...)
  check <- unlist(lapply(dots, inherits, "lm.rrpp"))
  if(any(!check)) stop("\nObjects must be lm.rrpp fits\n.")
  dot.names <- lapply(dots, function(x) x$LM$Terms[[3]])
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
        cat("Warning: The number of PCs requested exceeds 
            possible dimensions for some model fits\n")
        cat("PC number will be dropped initially to\n\n")
        print(rank.list)
        cat("\nand might be reduced further,
            if the tolerance was also adjusted. \n\n")
        pc.no <- p.list
      }
    }
    
    ll.args <- list(fit = dots[[1]], tol = tol, pc.no = pc.no[[1]])
    temp <- sapply(1:length(dots), function(j){
      ll.args$fit <- dots[[j]]
      ll.args$pc.no <- pc.no[[j]]
      do.call(logL, ll.args)
    }) 
    res <- unlist(temp[1,])
    rank <- unlist(temp[2,])
  } else res <- sapply(dots, cov.trace) 
  
  par.pen <- function(f){
    p <- f$LM$p.prime
    k <- f$LM$QR$rank
    2 * (p * k + 0.5 * p* (p + 1))
  }
  pp <- sapply(dots, par.pen)
  
  if(type == "logLik") {
    out <- cbind(res, rank, pp, -2*res + pp)
    colnames(out) <- c("logLik", "residual.pc.no", "penalty", "AIC")
  } else {
    out <- cbind(res, pp)
    colnames(out) <- c("cov.trace", "penalty")
  }
  rownames(out) <- dot.names
  out <- list(table = out, names = dot.names)
  attr(out, "class") = "model.comparison"
  out
  
}
