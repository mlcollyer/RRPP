#' Linear Model Evaluation with RRPP performed within subjects
#'
#' Function performs a linear model fit over many random permutations of 
#' data, using a randomized residual permutation procedure restricted to
#' subjects.
#'
#' The function fits a linear model using ordinary least squares (OLS) or 
#' generalized least squares (GLS) estimation of coefficients over any 
#' number of random permutations of the data, but the permutations are mostly 
#' restricted to occur with subject blocks for any model terms other than subjects.  
#' All functionality should resemble that of \code{\link{lm.rrpp}}.  However, 
#' an argument for research subjects is also required.  The purpose of this function 
#' is to account for the non-independence among observations of research subjects 
#' (due to sampling within subjects), while also allowing for the non-independence 
#' among subjects to be considered.  
#' 
#' By comparison, the covariance matrix option in \code{\link{lm.rrpp}} must have a
#' one-to-one match to observations, which can be matched by the row names of the data.
#' In this function, the covariance matrix can be the same one used in \code{\link{lm.rrpp}}
#' but the number of observations can be greater.  For example, if subjects are
#' species or some other level of taxonomic organization, data can comprise measurements
#' on individuals.  Users have the option to expand the covariance matrix for subjects
#' or input one they have generated.
#' 
#' Irrespective of covariance matrix type, the row names of the data matrix must match the
#' subjects.  This step assures that the analysis can proceed in \code{\link{lm.rrpp}}.  It
#' is also best to make sure to use an \code{\link{rrpp.data.frame}}, so that the subjects
#' can be a name in that data frame.  For example, if research subjects are species and
#' data (observations) are collected from individuals within species, then a procedure like 
#' the following should produce results:
#' 
#' rownames(Y) <- species
#' rdf <- rrpp.data.frame(Y = Y, subjects = species, x = x)
#' fit <- lm.rrpp.ws(Y ~ species * x, data = rdf, Cov = myCov, ...)
#' 
#' where ... means other arguments.  The covariances in the the Covariance matrix can be 
#' sorted by the subjects factor but data will not be sorted.  Therefore, names matching
#' the subjects is essential.  Additionally, subjects must be a factor in the data frame
#' or a factor in the global environment.  It cannot be part of a list.  Something like
#' subjects - mylist$species will not work.  Assuring that data and subjects are in the 
#' same \code{\link{rrpp.data.frame}} object as data is the best way to avoid errors.
#' 
#' Most attributes for this analysis are explained with \code{\link{lm.rrpp}}.  
#' The notable different attributes for this function are that: (1) a covariance 
#' matrix for the non-independence of subjects can be either a symmetric matrix 
#' that matches in dimensions the number of subjects or the number of observations; 
#' (2) a parameter (delta) that can range between near 0 and 1 to calibrate the 
#' covariances between observations of different subjects; and (3) a 
#' parameter (gamma) that is either 1 (equal) or the square-root of the subject 
#' sample size (sample) to calibrate the covariances among observations 
#' within subjects.  If delta = 0, it is expected that the covariance between
#' individual observations, between subjects, is the same as expected from the
#' covariance matrix, as if observations were the single observations made on subjects.
#' As delta approaches 1, the observations become more independent, as if it is 
#' expected that the many observations would not be expected to be
#' as correlated as if from one observation.  Increasing delta might be useful, if, 
#' for example, many individuals are sampled within species, from different locations,
#' different age groups, etc.  Alternatively, the sample size (n_i) for subject i 
#' can also influence the trend of inter-subject covariances.  If more individual
#' observations are sampled, the correlation between subjects might be favored
#' to be smaller compared to fewer observations.  The covariances can be adjusted
#' to allow for greater independence among observations to be assumed for larger samples.
#'
#' A design matrix, \bold{X}, is constructed with 0s and 1s to indicate subjects association,
#' and it is used to expand the covariance matrix (\bold{C}) by \bold{XCt(X)}, where \bold{t(X)}
#' is the matrix transpose.  The parameters in \bold{X} are multiplied by exp(-delta * gamma)
#' to scale the covariances.  (If delta = 0 and gamma = 1, they are unscaled.)  
#' 
#' These options for scaling covariances could be important for data with 
#' hierarchical organization.
#' For example, data sampled from multiple species with expected covariances 
#' among species based on phylogenetic distances, might be expected to not covary as
#' strongly if sampling encounters other strata like population, sex, and age.
#' An a priori expectation is that covariances among observations would be expected to
#' be smaller than between species, if only one observation per species were made.
#' 
#' If one wishes to have better control over between-subject and within-subject
#' covariances, based on either a model or empirical knowledge, a covariance matrix should 
#' be generated prior to analysis.  One can input a covariance matrix with dimensions 
#' the same as \bold{XCt(X)}, if they prefer to define covariances in an alternative way.
#' A function to generate such matrices based on separate
#' inter-subject and intra-subject covariance matrices is forthcoming.
#' 
#' IMPORTANT.  It is assumed that either the levels of the covariance matrix (if 
#' subject by subject) match the subject levels in the subject argument, or that
#' the order of the covariance matrix (if observation by observation) matches the order
#' of the observations in the data.  
#' No attempt is made to reorder a covariance matrix by observations
#' and row-names of data are not used to re-order the covariance matrix.  If the
#' covariance matrix is small (same in dimension as the number of subject levels), the function
#' will compile a large covariance matrix that is correct in terms of order, but
#' this is based on the subjects argument.
#' 
#' The covariance matrix is important for describing the expected covariances
#' among observations, especially knowing observations between and within subjects 
#' are not independent.  However, the randomization of residuals in a permutation
#' procedure (RRPP) is also important for testing inter-subject and 
#' intra-subject effects.  There are two RRPP philosophies used.  If the
#' variable for subjects is part of the formula, the subject effect is evaluated with
#' type III sums of squares and cross-products (estimates SSCPs between a model with all
#' terms and a model lacking subject term), and RRPP performed for all residuals of the
#' reduced model.  Effects for all other terms are evaluated with type II SSCPs and RRPP
#' restricted to randomization of reduced model residuals, within subject blocks.  This
#' assures that subject effects are held constant across permutations, so that intra-subject
#' effects are not confounded by inter-subject effects.
#' 
#' 
#' More details will be made.
#' 
#' The \code{\link{lm.rrpp}} arguments not available for this function include: 
#' full.resid, block, and SS.type.  These arguments are fixed because of
#' the within-subject blocking for tests, plus the requirement for type II SS
#' for within-subject effects.
#' 
#' 
#' @param f1 A formula for the linear model (e.g., y~x1+x2).  
#' @param subjects A required factor (or vector coercible to factor) of 
#' subjects that might have multiple data observations.
#' @param iter Number of iterations for significance testing
#' @param turbo A logical value that if TRUE, suppresses coefficient estimation 
#' in every random permutation.  This will affect subsequent analyses that 
#' require random coefficients (see \code{\link{coef.lm.rrpp}})
#' but might be useful for large data sets for which only ANOVA is needed.
#' @param seed An optional argument for setting the seed for random 
#' permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found 
#' for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  
#' One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param int.first A logical value to indicate if interactions of first 
#' main effects should precede subsequent main effects
#' @param RRPP A logical value indicating whether residual randomization 
#' should be used for significance testing
#' @param Cov An optional argument for including a covariance matrix 
#' to address the non-independence
#' of error in the estimation of coefficients (via GLS).  If included, 
#' any weights are ignored.  This matrix must math in dimensions either
#' the number of subjects or the number of observations.
#' @param delta A within-subject scaling parameter for covariances, ranging from 
#' 0 to 1.  If delta = 0, a sight value (0.001) is added to assure variances of the 
#' covariance matrix are 0.1 percent larger than covariances.
#' @param gamma A sample-size scaling parameter that is adjusted to be 1 ("equal")
#' scaling or the square-root of the sample size for subject observations ("sample").
#' @param data A data frame for the function environment, see 
#' \code{\link{rrpp.data.frame}}
#' @param print.progress A logical value to indicate whether a progress 
#' bar should be printed to the screen.
#' This is helpful for long-running analyses.
#' @param Parallel Either a logical value to indicate whether parallel processing 
#' should be used, a numeric value to indicate the number of cores to use, or a predefined
#' socket cluster.  This argument defines parallel processing via the \code{parallel} library. 
#' If TRUE, this argument invokes forking or socket cluster assignment of all processor cores, 
#' except one.  If FALSE, only one core is used. A numeric value directs the number of cores to 
#' use, but one core will always be spared.  If a predefined socket cluster (Windows) is provided,
#' the cluster information will be passed to \code{parallel}.
#' @param ... Arguments typically used in \code{\link{lm}}, such as 
#' weights or offset, passed on to
#' \code{lm.rrpp} for estimation of coefficients.  If both weights and 
#' a covariance matrix are included,
#' weights are ignored (since inverses of weights are the diagonal elements 
#' of weight matrix, used in lieu
#' of a covariance matrix.)
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{lm.rrpp.ws} is a list containing the 
#' following
#' \item{call}{The matched call.}
#' \item{LM}{Linear Model objects, including data (Y), coefficients, 
#' design matrix (X), sample size
#' (n), number of dependent variables (p), dimension of data space (p.prime),
#' QR decomposition of the design matrix, fitted values, residuals,
#' weights, offset, model terms, data (model) frame, random coefficients 
#' (through permutations),
#' random vector distances for coefficients (through permutations), 
#' whether OLS or GLS was performed, 
#' and the mean for OLS and/or GLS methods. Note that the data returned 
#' resemble a model frame rather than 
#' a data frame; i.e., it contains the values used in analysis, which 
#' might have been transformed according to 
#' the formula.  The response variables are always labeled Y.1, Y.2, ..., 
#' in this frame.}
#' \item{ANOVA}{Analysis of variance objects, including the SS type, 
#' random SS outcomes, random MS outcomes,
#' random R-squared outcomes, random F outcomes, random Cohen's f-squared 
#' outcomes, P-values based on random F
#' outcomes, effect sizes for random outcomes, sample size (n), number of 
#' variables (p), and degrees of freedom for
#' model terms (df).  These objects are used to construct ANOVA tables.}
#' \item{PermInfo}{Permutation procedure information, including the number 
#' of permutations (perms), The method
#' of residual randomization (perm.method), and each permutation's sampling 
#' frame (perm.schedule), which
#' is a list of reordered sequences of 1:n, for how residuals were 
#' randomized.}
#' \item{Models}{Reduced and full model fits for every possible model 
#' combination, based on terms
#' of the entire model, plus the method of SS estimation.}
#' @references TBD
#' @seealso \code{\link{lm.rrpp}}; 
#' @references ter Braak, C.J.F. 1992. Permutation versus bootstrap significance tests in 
#' multiple regression and ANOVA. pp .79–86 In Bootstrapping and Related Techniques. eds K-H. Jockel, 
#' G. Rothe & W. Sendler.Springer-Verlag, Berlin.
#' \code{\link[stats]{lm}} for more on linear model fits.
#' @examples 
#' # TBD
lm.rrpp.ws <- function(f1, subjects, 
                       iter = 999, turbo = FALSE, 
                       seed = NULL, int.first = FALSE,
                       RRPP = TRUE, 
                       data = NULL, Cov = NULL,
                       delta = 0, 
                       gamma = c("sample", "equal"),
                       print.progress = FALSE, Parallel = FALSE, ...) {
  
  sub.var.name <- deparse(substitute(subjects))
  if(!is.null(data)) {
    subj.no <- which(names(data) %in% sub.var.name)
    subj <- data[[subj.no]]
    if(!is.null(subj)) subjects <- subj
  }
  
  L <- L.args <- c(as.list(environment()), list(...))
  names(L)[which(names(L) == "f1")] <- "formula"
  
  if(inherits(f1, "formula")) {
    exchange.args <- lm.args.from.formula(L)
    if(!is.null(exchange.args$D)) D <- exchange.args$D
    exchange.args <- exchange.args[c("Terms", "Y", "model")]
    Terms <- exchange.args$Terms
    Y <- exchange.args$Y
    
  } else stop("\n lm.rrpp.ws currently requires a formula rather than an lm object.\n",
              call. = FALSE)
  L <- NULL
  Ynames <- rownames(Y)

  n <- NROW(Y)
  if(is.null(L.args$data)) L.args$data <- rrpp.data.frame(Y = Y) else
    L.args$data$Y <- Y
  L.args$f1 <- update(L.args$f1, Y ~ .)
  Y <- NULL
  
  if(!identical(Ynames, as.character(subjects)))
    stop("The rownames of the data must match the subjects.\n", 
         call. = FALSE)
  
  STerm <- which(attr(Terms, "term.labels") == sub.var.name)
  subTest <- !(length(STerm) == 0)
  subjects <- as.factor(subjects)
  sub.lev <- levels(subjects)
  Xsub <- model.matrix(~ subjects + 0)
  colnames(Xsub) <- sub.lev
  
  if(!is.null(Cov)) {
    
    cov.lev <- levels(as.factor(dimnames(Cov)[[1]]))
    if(!all(sub.lev %in% cov.lev) || !all(cov.lev %in% sub.lev))
      stop("\nThere is a mismatch between subject levels and covariance levels\n
         Make sure there is direct correspondence between subjects in the model\n
         and subjects in the covariance matrix.  Check spelling of names.\n", 
           call. = FALSE)
    Xs <- Matrix(Xsub, sparse = TRUE)
    nn <- nrow(Cov)
    if(nn == nrow(Xsub)) cov.type <- 1 else 
      if(nn == length(sub.lev)) cov.type <- 2 else
        stop(paste("There is an irreconcilable covariance matrix, in terms of",
                   "\nsubject number or subject factor levels.  Please make sure that",
                   "\nthe number of Cov rows or columns match the number of subject levels",
                   "\nor the number of observations.", sep = " "), call. = FALSE)
    if(cov.type == 2) {
      Cov <- Cov[colnames(Xsub), colnames(Xsub)]
      Cov.names <- dimnames(Cov)[[1]]
      Cov <- Matrix(Cov, sparse = TRUE)
      if(is.vector(delta) && length(delta) > 1 &&
         length(delta) != length(sub.lev)) {
        stop(paste("Either a single value or vector equal in length to the number",
                   "of subject levels is required.\n", sep = " "), call. = FALSE)
        
      }
      if(!is.numeric(delta))
        stop("delta must be a numeric value or vector.\n", call. = FALSE)
      
      delta[which(delta <= 0)] <- 0.001
      delta[which(delta > 1)] <- 1
      
      n.list <- as.vector(by(subjects, subjects, length))
      gamma <- match.arg(gamma)
      gam <- if(gamma == "equal") 1 else 
        sqrt(n.list)
      
      Xr <- t(t(Xs) * exp(-delta * gam))
      CovEx <- Xr %*% Cov %*% t(Xr)
      D <- diag(Cov)
      D <- unlist(lapply(1:length(D), function(j) 
        rep(D[j], n.list[j])))
      diag(CovEx) <- D
      Xs <- Xr <- NULL
    } else {
      CovEx <- Cov
      cat("\nIt is assumed that the Covariance matrix is ordered the same as the data...\n")
    }
    
    Xsub <- NULL
    
  } else CovEx <- NULL
  
  L.args$Cov <- CovEx
  if(!is.null(L.args$Cov)) rownames(L.args$data$Y) <- rownames(L.args$Cov) <-
    colnames(L.args$Cov) <- 1:n
  L.args$SS.type <- "III"
  L.args$block <- NULL
  L.args$turbo <- TRUE
  CovEx <- Cov <- NULL
  
  
  if(subTest){
    
    if(print.progress)
      cat("\nPerforming among-subjects analysis...\n")
    
    fit.subjects <- suppressWarnings( do.call(lm.rrpp, L.args) )
    
    sub.rpl <- which(rownames(fit.subjects$ANOVA$SS) == sub.var.name)
    sSS <- fit.subjects$ANOVA$SS[sub.rpl, ]
    sMS <- fit.subjects$ANOVA$MS[sub.rpl, ]
    sRSS <- fit.subjects$ANOVA$RSS[sub.rpl, ]
    sTSS <- fit.subjects$ANOVA$TSS[sub.rpl, ]
    sRSS.model <- fit.subjects$ANOVA$RSS.model[sub.rpl, ]
    sRsq <- fit.subjects$ANOVA$Rsq[sub.rpl, ]
    sFs <- fit.subjects$ANOVA$Fs[sub.rpl, ]
    scohenf <- fit.subjects$ANOVA$cohenf[sub.rpl, ]
    sRed <- fit.subjects$Models$reduced[[sub.rpl]]
    sFull <- fit.subjects$Models$full[[sub.rpl]]
    fit.subjects <- NULL
  }
  
  L.args$SS.type <- "II"
  L.args$block <- subjects
  L.args$turbo <- turbo

  if(print.progress)
    cat("\nPerforming within-subjects analysis...\n")
  
  out <- suppressWarnings( do.call(lm.rrpp, L.args) )
  L.args <- NULL
  
  if(subTest){
    out$ANOVA$SS.type <- "Within-subject type II"
    out$ANOVA$SS[sub.rpl, ] <- sSS
    out$ANOVA$MS[sub.rpl, ] <- sMS
    out$ANOVA$RSS[sub.rpl, ] <- sRSS
    out$ANOVA$TSS[sub.rpl, ] <- sTSS
    out$ANOVA$RSS.model[sub.rpl, ] <- sRSS.model
    out$ANOVA$Rsq[sub.rpl, ] <- sRsq
    out$ANOVA$Fs[sub.rpl, ] <- sFs
    out$ANOVA$cohenf[sub.rpl, ] <- scohenf
    out$Models$reduced[[sub.rpl]] <- sRed
    out$Models$full[[sub.rpl]] <- sFull
  }
  
  rownames(out$LM$Y) <- Ynames
  if(!is.null(out$LM$Cov)){
    rownames(out$LM$Cov) <- colnames(out$LM$Cov) <- Ynames
  }

  out$call <- match.call()
  out$subjects <- subjects
  out$subjects.var <- sub.var.name
  if(!is.null(Cov)) dimnames(out$LM$Cov) <- list(subjects, subjects)
  class(out) <- c("lm.rrpp.ws", class(out))
  out

}

