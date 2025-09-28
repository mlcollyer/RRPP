#' Mixed Linear Model Evaluation with RRPP performed within subjects
#'
#' Function performs a mixed linear model fit over many random permutations of 
#' data, using a randomized residual permutation procedure restricted to
#' subjects.  This function is likely to evolve greater flexibility in the future.
#'
#' The function fits a mixed linear model using ordinary least squares (OLS) or 
#' generalized least squares (GLS) estimation of coefficients for fixed effects, and 
#' maximum likelihood (ML) or restricted ML (REML) for random effects, over any 
#' number of random permutations of the data, but the permutations are mostly 
#' restricted to occur with subject blocks for any model terms other than subjects.  
#' Most functionality should resemble that of \code{\link{lm.rrpp}} and \code{\link{lm.rrpp.ws}},  
#' with some caveats.  First, estimation of random effects uses the methodology of the 
#' \code{\link[lme4]{lmer}} function (Bates et al., 2015), 
#' for only one random term based on a subjects designation.  
#' Second, this means that a "random slopes" model involves an interaction between subjects and a numeric 
#' variable (not a categorical variable).  Third, although any number of fixed effects can be used,
#' only one fixed effect can interaction with subjects. Finally, the function is intended to allow
#' analysis of variance (ANOVA), like with \code{\link{lm.rrpp.ws}}, but not treated random 
#' effects as fixed.  
#' There is not a coefficients test associated eith this function, as with other \code{\link{lm.rrpp}}
#' that include only fixed effects.
#' 
#' Although this function allows covariance matrices to be used for generalized least squares (GLS)
#' estimation of fixed effects, it does not currently 
#' expand the covariance mnatrix based only on subjects, as does 
#' \code{\link{lm.rrpp.ws}}.  However, that function can be used to obtain a covariance matrix for this 
#' function.  The covariance structure (Lambda matrix) among subjects is estimated via ML or REML (Bates et al., 2015) and 
#' is not impacted by the expected covariances among fixed effects, however, the expected covariances
#' will influence interactions between fixed and random effects.  The Lambda matrix is estimated as a single matrix,
#' even for multivariate data.  This matrix is obtained by projecting data onto the first component of
#' a singular value decomposition of the product between the data and a hat matrix of the the fixed effects 
#' (basically projection following a principal component analysis of model fitted values).  This is a computationally
#' proficient step that is similar to averaging covariance parameters (theta) over all variables.
#' If independent Lambda matrices are desired, then this function can be used on, e.g., each principal
#' component of data, followed by combining the coefficients.  For ANOVA, the precise estimation of Lambda
#' is not as important for tests of fixed effects, as it means a linear transformation of coefficients
#' that would be consistent across random permutations of RRPP.  It might be more important for model 
#' comparison and selection, in which case whether to use LS, ML, REML, and whether to use separate Lambda matrices
#' is a much deeper consideration.
#' 
#' A data frame (preferably an \code{\link{rrpp.data.frame}} object) must be
#' used for this function.  Not all "downstream" functions (e.g., logLik, predict) are guaranteed 
#' to work as expected, as this function is still very much under development.
#' 
#' The \code{\link{lm.rrpp}} arguments not available for this function include: 
#' full.resid, block, and SS.type.  These arguments are fixed because of
#' the within-subject blocking for tests, plus the requirement for type II SS
#' for within-subject effects.
#' 
#' @param fixed A formula for the fixed effects of linear model (e.g., y~x1+x2).  
#' @param subjects A variable that can be found in the data frame indicating the research subjects
#' for the analysis.  This variable must be in the data frame.  Is can be either numeric 
#' (if its slot in the data frame is known) or a character, e.g., "sub_id".  It is imperative that
#' it is ordered the same as the data but that the data do not have row names the same as subjects.
#' For example, the subjects variable in the data frame might be sub_id: sub1, sub1, sub1, sub2,
#' sub2, sub2, ... and the row names of the data might be obs1, obs2, obs3, obs4, obs5, obs6, ...  
#' The data do not need to have row names but the subjects variable has to be provided.
#' @param type An indication whether to use random intercepts ("intercepts") or random intercepts 
#' and slopes ("slopes") for random effects associated with subjects.
#' @param slopeTerm An otional character value (e.g., "mass") indicating which term in the data frame is to be
#' used in the estimation of random slopes.  This term should also be in the fixed formula.
#' @param estimation One of "LS", "ML", or "REML", for guiding how random effects should be estimated.
#' @param iter Number of iterations for significance testing
#' @param data A data frame for the function environment, see 
#' \code{\link{rrpp.data.frame}}.  A data frame is required for this analysis.
#' @param ... Arguments typically used in \code{\link{lm}}, \code{\link{lm.rrpp}},
#' or \code{\link{lm.rrpp.ws}} such as weights or offset, seed for permutations, etc., passed on to
#' \code{lm.rrpp.ws} for estimation of coefficients.  
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{lmm.rrpp} is a list containing the 
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
#' in this frame.  Most important objects for mixed models are also found here. }
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
#' @seealso \code{\link{lm.rrpp}}; \code{\link{measurement.error}}
#' @references Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear
#' Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.
#' doi:10.18637/jss.v067.i01.
#' 
#' \subsection{Limitations}{
#' This function currently only works for "Random Intercepts" and "Random slopes" models,
#' where the intercepts or slopes of a single subject term (random effect) are free to vary.  
#' If a Random Slopes model is used, a single covariate (numeric variable) must be specified for slopes.  
#' Generally, this is term that is also among the fixed effects.  These limitations are imposed
#' to allow ML or REML estimation.  More complex models can be entertained with \code{lmm.rrpp.ws},
#' but using only "LS" estimation.  If "LS" estimation is used in this function, results should be the same
#' with a comparable formula, using \code{lmm.rrpp.ws}
#' 
#' If a ML or REML fit is near singular (one or more covariance parameters are 0), the analysis will be arrested
#' rather than continued with some random effects rendered as 0.  This should be an indication that the model
#' is over-fit with the type of estimation used. 
#' 
#' These limitations might be addressed with future versions of this function, as theoretical advances are made.
#'  }

#' @examples 
#' \dontrun{
#' # Random intercepts models
#' 
#' data(ReptAbund)
#' fit_ls <- lmm.rrpp(fixed = reptile_abundance ~ Year, 
#' subjects = "Site",      
#' type = "intercepts", 
#' estimation = "LS", 
#' data = ReptAbund)
#' 
#' fit_ml <- lmm.rrpp(fixed = reptile_abundance ~ Year, 
#' subjects = "Site",      
#' type = "intercepts", 
#' estimation = "ML", 
#' data = ReptAbund)
#'
#' fit_reml <- lmm.rrpp(fixed = reptile_abundance ~ Year, 
#' subjects = "Site",      
#' type = "intercepts", 
#' estimation = "REML", 
#' data = ReptAbund) 
#' 
#' anova(fit_ls, error = c("Site", "Residuals"))
#' anova(fit_ml, error = c("Site", "Residuals"))
#' anova(fit_reml, error = c("Site", "Residuals"))
#' 
#' coef(fit_ls, type = "list")
#' coef(fit_ml, type = "list")
#' coef(fit_reml, type = "list")
#' 
#' # Random slopes models
#' 
#' fit_ls <- lmm.rrpp(fixed = reptile_abundance ~ Year, 
#' subjects = "Site",      
#' type = "slopes",
#' slopeTerm = "Year", 
#' estimation = "LS", 
#' data = ReptAbund)
#' 
#' fit_ml <- lmm.rrpp(fixed = reptile_abundance ~ Year, 
#' subjects = "Site",      
#' type = "slopes",
#' slopeTerm = "Year", 
#' estimation = "ML", 
#' data = ReptAbund)
#'
#' # should not work!
#' 
#' fit_reml <- lmm.rrpp(fixed = reptile_abundance ~ Year, 
#' subjects = "Site",      
#' type = "slopes",
#' slopeTerm = "Year",  
#' estimation = "REML", 
#' data = ReptAbund) 
#' 
#' anova(fit_ls, error = c("Year:Site", "Year:Site", "Residuals"))
#' anova(fit_ml, error = c("Year:Site", "Year:Site", "Residuals"))
#' 
#' coef(fit_ls, type = "list")
#' coef(fit_ml, type = "list")
#'  }

lmm.rrpp <- function(fixed,
                     subjects,
                     type = c("intercepts", "slopes"),
                     slopeTerm = NULL,
                     estimation = c("REML", "ML", "LS"),
                     iter = 999,
                     data,
                     ...) {
  estimation <- match.arg(estimation)
  dots <- list(...)
  largs <- as.vector(names(formals(lm.rrpp.ws)))
  lm.rrpp.args <- dots[names(dots) %in% largs]
  names(lm.rrpp.args) <- lm.rrpp.args[largs %in% names(dots)]
  
  response <- as.character(fixed[[2]])
  trms <- attr(terms(fixed), "term.labels")
  form <- if(type == "intercepts") 
    reformulate(response = response, 
                termlabels = c(trms, subjects)) else
                  reformulate(response = response, 
                              termlabels = c(trms, subjects,
                                             paste(slopeTerm, subjects, sep = ":")))
  
  lm.rrpp.args$f1 <- form
  lm.rrpp.args$subjects <- subjects
  lm.rrpp.args$verbose <- FALSE
  lm.rrpp.args$iter <- 0
  lm.rrpp.args$data <- data
  if(isTRUE(dots$print.progress)){
    cat("\nInitial lm.rrpp.ws fit:\n ")
  }
  
  fit <- suppressWarnings(do.call(lm.rrpp.ws, 
                                  lm.rrpp.args))
  lm.rrpp.args$f1 <- fixed
  
  fixed.fit <- suppressWarnings(do.call(lm.rrpp.ws, 
                                        lm.rrpp.args))
  
  fit$call <- match.call()
  dat <- fit$LM$data
  dat$resp <- as.matrix(dat$Y[,1])
  
  rands <- if(type == "intercepts")
    paste("(1|", subjects, ")", sep = "") else
      paste("(", slopeTerm, "|", subjects, ")", sep = "")
  
  randform <- reformulate(response = response, 
                          termlabels = c(trms, rands))
  randform <- update(randform, resp ~ .)
  
  init.rand.fit <- suppressMessages(suppressWarnings(lmer(randform, 
                        REML = (estimation == "REML"),
                        data = dat)))
  
  Y <- dat$Y
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  
  if(fit$LM$gls) {
    if(!is.null(fit$LM$Cov)) {
      Pcov <- getModelCov(fit, "Pcov")
      Y <- Pcov %*% Y
    } else {
      w <- sqrt(fit$LM$weights)
      Y <- Y * w
    }
  }
  
  subjFactor <- unlist(dat[[which(names(dat) == subjects)]])
  ns <- nlevels(subjFactor)
  subbjLevels <- levels(subjFactor)
  
  if(type == "slopes") {
    slopeVar <- unlist(dat[[which(names(dat) == slopeTerm)]])
    if(is.factor(slopeVar) && estimation != "LS"){
      stop("\nThe slope term must be numeric for ML or REML fits...",
           "\nThis is a requirement of the lme4::lmer function dependency. ",
           "\nConsider using dummy (binary) variables for two-group ",
           "comparisons. \nSee the lme4::dummy function for details. \n",
           call. = FALSE)
    }
  }
  
  Xs <- suppressWarnings(
    getModels(fixed.fit, "qr"))
  k <- length(Xs$full)
  QR <- Xs$full[[k]]
  
  for(i in 1:2){
    reduced <- lapply(Xs$reduced, 
                      function(x) as.matrix(x$X))
    full <- lapply(Xs$full, 
                   function(x) as.matrix(x$X))
  }
  Xs$reduced <- reduced
  Xs$full <- full
  rm(reduced, full)
  
  Xs$reduced[[k + 1]] <- Xs$full[[k]]
  Xs$full[[k + 1]] <- Xs$full[[k]]
  names(Xs$reduced)[[k + 1]] <- 
    names(Xs$full)[[k + 1]] <- subjects
  if(type == "slopes"){
    Xs$reduced[[k + 2]] <- Xs$full[[k]]
    Xs$full[[k + 2]] <- Xs$full[[k]]
    names(Xs$reduced)[[k + 2]] <- 
      names(Xs$full)[[k + 2]] <- paste(
        slopeTerm, subjects, sep = ":")
  }
  
  Zall <- getME(init.rand.fit, 
                "Z")
  if(type == "slopes"){
    nzc <- ncol(Zall)
    Zslope <- Zall[, seq(2, nzc, 2)]
    Zsub <- Zall[, seq(1, nzc - 1, 2)]
  } else {
    Zsub <- Zall
    Zslope <- NULL
  }
  
  Snames <- rand.coef.names <- dimnames(Zsub)[[2]]
  if(type == "slopes"){
    rand.coef.names <-paste(
      rep(Snames, each = 2), 
      c("", paste(":", slopeTerm, sep = "")), sep = "")
    
  }
  
  
  if(estimation != "LS"){
    
    if(p > 1){
      H <- tcrossprod(QR$Q)
      yp <- ordinate(Y, A = H, rank. = 1)$x
      ft <- refit(init.rand.fit, yp)
    } else ft <- init.rand.fit
    
    if(isSingular(ft))
      stop("Model is near-singular.  Try a different estimation method or a simpler model.\n",
           call. = FALSE)
    
    Lambda <- getME(ft, "Lambda")
    
    if(!is.null(Zslope)){
      subs <- seq(1, nrow(Lambda) - 1, 2)
      Lambda_sub <- Lambda[subs, subs]
      Lambda_slopes <- Lambda[-subs, -subs]
    } else Lambda_sub <- Lambda
    
    Zsub <- Zsub %*% Lambda_sub
    Zall <- Zall %*% Lambda
    if(!is.null(Zslope)) Zslope <- Zslope %*% Lambda_slopes
    
  } else {
    Lambda <- Lambda_sub <- Lambda_slopes <- NULL
  }
  
  Hbs.reduced <- Hbs.full <- XZs.reduced <- XZs.full <- 
    Qs.reduced <- Qs.full <- list()
  
  for(i in 1:k){
    Hbs.reduced[[i]] <- getLMM_Hb(Xs$reduced[[i]], Zsub) 
    Hbs.full[[i]] <- getLMM_Hb(Xs$full[[i]], Zsub)
    XZs.reduced[[i]] <- cbind(Xs$reduced[[i]], Zsub)
    XZs.full[[i]] <- cbind(Xs$full[[i]], Zsub)
    Qs.reduced[[i]] <- getLMMQ(XZs.reduced[[i]], Hbs.reduced[[i]])
    Qs.full[[i]] <- getLMMQ(XZs.full[[i]], Hbs.full[[i]])
  }
  
  
  # If Intercepts model
  
  if(type == "intercepts"){
    Hbs.reduced[[k+1]] <- getLMM_Hb(Xs$full[[k]])
    XZs.reduced[[k+1]] <- Xs$full[[k]]
    Hbs.full[[k+1]] <- getLMM_Hb(Xs$full[[k]], Zsub)
    XZs.full[[k+1]] <- cbind(Xs$full[[k]], Zsub)
    Qs.reduced[[k+1]] <- getLMMQ(XZs.reduced[[k+1]], 
                                 Hbs.reduced[[k+1]])
    Qs.full[[k+1]] <- getLMMQ(XZs.full[[k+1]], 
                                 Hbs.full[[k+1]])
  }
  
  # If slopes model
  
  if(type == "slopes"){
    Hbs.reduced[[k+1]] <- getLMM_Hb(Xs$full[[k]], Zslope)
    XZs.reduced[[k+1]] <- cbind(Xs$full[[k]], Zslope)
    Hbs.full[[k+1]] <- getLMM_Hb(Xs$full[[k]], Zall)
    XZs.full[[k+1]] <- cbind(Xs$full[[k]], Zall)
    Hbs.reduced[[k+2]] <- getLMM_Hb(Xs$full[[k]], Zsub)
    XZs.reduced[[k+2]] <- cbind(Xs$full[[k]], Zsub)
    Hbs.full[[k+2]] <- getLMM_Hb(Xs$full[[k]], Zall)
    XZs.full[[k+2]] <- cbind(Xs$full[[k]], Zall)
    
    Qs.reduced[[k+1]] <- getLMMQ(XZs.reduced[[k+1]], 
                                 Hbs.reduced[[k+1]])
    Qs.full[[k+1]] <- getLMMQ(XZs.full[[k+1]], 
                              Hbs.full[[k+1]])
    Qs.reduced[[k+2]] <- getLMMQ(XZs.reduced[[k+2]], 
                                 Hbs.reduced[[k+2]])
    Qs.full[[k+2]] <- getLMMQ(XZs.full[[k+2]], 
                              Hbs.full[[k+2]])
  }
  

  Fitted <- lapply(Qs.reduced, function(q){
    fastFit(q, Y, n, p)
  })
  
  Resid <- lapply(Fitted, function(f) Y - f)
  
  
  ### start RRPP
  
  ind <- perm.index(n, iter, block = subjFactor,
                    seed = dots$seed)
  
  ind_s <- perm.index(n, iter, block = NULL,
                      seed = dots$seed)
  
  perms <- length(ind)
  kk <- length(XZs.full)
  
  XZ_full <- XZs.full[[kk]]
  XZ_null <- XZ_full[, 1] # need to consider no intercept
  Hb_full <- Hbs.full[[kk]]
  QRnull <- QRforX(matrix(1, n, 1))
  Hb_null <- getHb(QRnull)
  Qnull <- QRnull$Q
  QRfull <- QRforX(XZ_full %*% Hb_full)
  Qfull <- QRfull$Q[, 1:QRfull$rank]
  
  B <- as.matrix(Hb_full %*% Y)
  dimnames(B) <- list(c(colnames(fixed.fit$LM$X),
                        rand.coef.names),
                      colnames(Y))
  U <- as.matrix(B[-seq_len(ncol(fixed.fit$LM$X)), ])
  UtU_det <- try(det(crossprod(U)), silent = TRUE)
  if(inherits(UtU_det, "try-error"))
    UtU_det <- NA
  
  sub_no <- which(names(Xs$full) == subjects)
  
  Result <- as.list(array(NA, perms))
  names(Result) <- names(ind)
  
  ss <- function(qf, qr, y) {
    ssr <- sum(crossprod(qr, y)^2)
    ssf <- sum(crossprod(qf, y)^2)
    ssm <- sum(y^2) - sum(crossprod(Qfull, y)^2)
    c(ssf - ssr, ssm)
  } 
  
  trmnms <- names(Xs$full)
  
  for(i in 1:perms){
    s <- ind[[i]]
    Y_i <- lapply(1:kk, function(j){
      sj <- if(j == sub_no) ind_s[[i]] else s
      Fitted[[j]] + as.matrix(Resid[[j]][sj, ])
    })
    yy <- Y[s, ]
    
    res <- vapply(1:max(1, kk), function(j){
      ss(Qs.full[[j]], Qs.reduced[[j]], Y_i[[j]])
    }, numeric(2))
    
    dimnames(res) <- list(c("SS", "RSS"), trmnms)
    
    Falls <- Map(function(y){
      fastFit(Qfull, y, n, p)}, Y_i)
    Ralls <- Map(function(y, f){
      y - f}, Y_i, Falls)
    
    RSS <- sapply(Ralls, function(r) sum(r^2))
    
    Fm <- fastFit(Qfull, yy, n, p)
    Rm <- yy - Fm
    RSS.model <- sum(Rm^2)
    F0 <- fastFit(Qnull, yy, n, p)
    R0 <- yy - F0
    TSS <- sum(R0^2)
    Result[[i]] <- list(
      SS = res,
      RSS.model = RSS.model,
      TSS = TSS
    )
    
  }
  
  
  
  SS <- sapply(Result, function(x) x$SS[1,])
  RSS <- sapply(Result, function(x) x$SS[2,])
  RSS.model <- sapply(Result, function(x) x$RSS.model)
  RSS.model <- matrix(RSS.model, kk, perms, byrow = TRUE)
  TSS <- sapply(Result, function(x) x$TSS)
  TSS <- matrix(TSS, kk, perms, byrow = TRUE)
  dimnames(TSS) <- dimnames(RSS.model) <- dimnames(SS)
  Df <- fit$ANOVA$df
  MS <- SS / Df[1:kk]
  Rsq <- SS / TSS
  
  fit$LM$LMM <- TRUE
  fit$LM$X <- fixed.fit$LM$X
  fit$LM$Z <- getME(init.rand.fit, "Z")
  fit$LM$Lambda <- Lambda
  fit$LM$QR <- NULL
  fit$LM$coefficients <- B
  fit$LM$coef.fixed <- seq_len(ncol(fixed.fit$LM$X))
  fit$LM$coef.random <- seq_len(NROW(B))[-fit$LM$coef.fixed]
  if(fit$LM$gls) 
    fit$L$gls.fitted <- XZ_full %*% fit$LM$coefficients else 
      fit$LM$fitted <- XZ_full %*% fit$LM$coefficients
  if(fit$LM$gls) 
    fit$L$gls.residuals <- Y - fit$LM$fitted else 
      fit$LM$residuals <- Y - fit$LM$fitted
  fit$LM$fixed <- fixed
  fit$LM$estimation <- estimation
  fit$LM$ranef.type <- type
  fit$LM$ranef.slopeTerm <- slopeTerm
  fit$LM$lm_form <- form
  fit$LM$lmer_form <- randform
  fit$LM$cnms <- init.rand.fit@cnms
  fit$LM$UtU_det <- UtU_det
  
  fit$ANOVA$SS <- SS
  fit$ANOVA$RSS <- RSS
  fit$ANOVA$RSS.model <- RSS.model
  fit$ANOVA$TSS <- TSS
  fit$ANOVA$Rsq <- Rsq
  
  fit$PermInfo$perms <- perms
  
  class(fit) <- c("lmm.rrpp", class(fit))
  fit$call <- match.call()
  fit
}


getLMMQ <- function(XZ, Hb){
  QR <- QRforX(XZ %*% Hb)
  Q <- QR$Q[, 1:QR$rank]
  Q
}

# Not used but could in the future
getSubBlocks <- function(flmer){
  Zt <- getME(flmer, "Ztlist")[[1]]
  block <- as.factor(Zt@i + 1)
  block
}

# Not used but could in the future
indexZtoX <- function(X, Z){
  Z <- Z0 <- as.matrix(Z)
  kx <- ncol(X)
  kz <- ncol(Z)
  nx <- nrow(X)
  nz <- nrow(Z)
  if(kx != kz || nx != nz)
    stop("Dimensions of matrices do not match.\n",
         call. = FALSE)
  res <- array(NA, kx)
  z.index <- 1:kz
  
  for(i in seq_len(kx)){
    x <- X[, i]
    r <- sapply(1:kz, function(j){
      identical(x, Z[, j])
    })
    a <- which(r)
    res[i] <- z.index[a]
    Z <- as.matrix(Z[, -a])
    z.index <- z.index[-a]
    kz <- kz - 1
  }
  
  res
}