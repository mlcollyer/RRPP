#' @name RRPP-package
#' @docType package
#' @aliases RRPP
#' @title Linear Model Evaluation with Randomized Residual Permutation Procedures
#' @author Michael Collyer and Dean Adams
#' @return Key functions for this package:
#' \item{\code{\link{lm.rrpp}}}{Fits linear models, using RRPP.}
#' \item{\code{\link{anova.lm.rrpp}}}{ANOVA on linear models, using RRPP, plus model comparisons.}
#' \item{\code{\link{coef.lm.rrpp}}}{Extract coefficients or perform test on coefficients, using RRPP.}
#' \item{\code{\link{predict.lm.rrpp}}}{Predict values from lm.rrpp fits and generate bootstrapped confidence intervals.}
#' \item{\code{\link{pairwise}}}{Perform pairwise tests, based on lm.rrpp model fits.}
#' 
#' @description Functions in this package allow one to evaluate linear models with residual randomization.
#' The name, "RRPP", is an acronym for, "Randomization of Residuals in a Permutation Procedure."  Through
#' the various functions in this package, one can use randomization of residuals to generate empirical probability
#' distributions for linear model effects, for high-dimensional data or distance matrices.
#' 
#' An especially useful option of this package is to fit models with either ordinary or generalized
#' least squares estimation (OLS or GLS, respectively), using theoretic covariance matrices.  Mixed linear
#' effects can also be evaluated.
#' 
#' @import parallel
#' @import stats
#' @import graphics
#' @import utils
#' @export print.lm.rrpp
#' @export summary.lm.rrpp
#' @export print.summary.lm.rrpp
#' @export print.anova.lm.rrpp
#' @export summary.anova.lm.rrpp
#' @export print.coef.lm.rrpp
#' @export summary.coef.lm.rrpp
#' @export print.predict.lm.rrpp
#' @export summary.predict.lm.rrpp
#' @export plot.lm.rrpp
#' @export plot.predict.lm.rrpp
#' @export print.pairwise
#' @export summary.pairwise
#' @export print.summary.pairwise
#' @export print.trajectory.analysis
#' @export summary.trajectory.analysis
#' @export print.summary.trajectory.analysis
#' 
#' @section RRPP TOC:
#' RRPP-package
NULL

#' Landmarks on pupfish
#'
#' @name Pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprinodon pecosensis body shapes, with indication of Sex and
#' Population from which fish were sampled (Marsh or Sinkhole).
#' @details These data were previously aligned with GPA.  Centroid size (CS) is also provided.  
#' See the \pkg{geomorph} package for details.
#' 
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic
#' change for phenotypes described by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
NULL

#' Landmarks on pupfish heads
#'
#' @name PupfishHeads
#' @docType data
#' @author Michael Collyer
#' @description Landmark data from Cyprinodon pecosensis head shapes, with variables for 
#' sex, month and year sampled, locality, head size, and coordinates of landmarks for head shape,
#' per specimen.  These data are a subset of a larger data set.
#' @details The variable, "coords", are data that were previously aligned
#' with GPA.  The variable, "headSize", is the Centroid size of each vector of coordinates.
#' See the \pkg{geomorph} package for details.
#' @references Gilbert, M.C. 2016. Impacts of habitat fragmentation on the cranial morphology of a 
#' threatened desert fish (Cyprinodon pecosensis). Masters Thesis, Western Kentucky University.
NULL

#' Plethodon comparative morphological data 
#'
#' @name PlethMorph
#' @docType data
#' @author Michael Collyer and Dean Adams
#' @keywords datasets
#' @description Data for 37 species of plethodontid salamanders.  Variables include snout to vent length
#' (SVL) as species size, tail length, head length, snout to eye length, body width, forelimb length,
#' and hind limb length, all measured in mm.  A grouping variable is also included for functional guild size.  
#' The data set also includes a phylogenetic covariance matrix based on a Brownian model of evolution, to assist in 
#' generalized least squares (GLS) estimation.
#' @details The covariance matrix was estimated with the vcv.phylo function of the R package, ape, based on the tree
#' described in Adams and Collyer (2018).
#' @references Adams, D.C and Collyer, M.L. 2018. Multivariate phylogenetic anova: group-clade aggregation, biological 
#' challenges, and a refined permutation procedure. Evolution. In press.
NULL

#' Simulated motion paths
#'
#' @name motionpaths
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @keywords datasets
NULL

#####----------------------------------------------------------------------------------------------------

# HELP FUNCTIONs

#' Create a data frame for lm.rrpp analysis
#'
#' Create a data frame for lm.rrpp analysis, when covariance or distance matrices are used
#'
#' This function is not much different than \code{\link{data.frame}} but is more flexible to allow
#' distance matrices and covariance matrices to be included.  Essentially, this function creates a list,
#' much like an object of class \code{data.frame} is also a list.  However, \code{rrpp.data.frame} is
#' less concerned with coercing the list into a matrix and more concerned with matching the number of observations (n).
#' It is wise to use this function with any \code{lm.rrpp} analysis so that \code{\link{lm.rrpp}} does not have to search
#' the global environment for needed data.
#'
#' It is assumed that multiple data sets for the same subjects are in the same order.
#'
#' See \code{\link{lm.rrpp}} for examples.
#'
#' @param ... Components (objects) to combine in the data frame.
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @examples
#' # Why use a rrpp.data.frame?
#' y <- matrix(rnorm(30), 10, 3)
#' x <- rnorm(10)
#' df <- data.frame(x = x, y = y)
#' df
#' rdf <- rrpp.data.frame(x = x, y = y)
#' rdf # looks more like a list
#' 
#' is.list(df)
#' is.list(rdf)
#' 
#' d <- dist(y) # distance matrix as data
#' 
#' # One can try this but it will result in an error
#' # df <- data.frame(df, d = d) 
#' rdf <- rrpp.data.frame(rdf, d = d) # works
#' 
#' fit <- lm.rrpp(d ~ x, data = rdf)
#' summary(fit)

rrpp.data.frame<- function(...){
  dots <- list(...)
  if(length(dots) == 1 && is.data.frame(dots[[1]])) {
    dots <- dots[[1]]
    class(dots) <- "rrpp.data.frame"
  } else if(length(dots) == 1 && inherits(dots[[1]], "geomorph.data.frame")) {
    dots <- dots[[1]]
    cat("\nWarning: Some geomorph.data.frame objects might not be compatible with RRPP functions.")
    cat("\nIf any part of the geomorph.data.frame conatins a 3D array,")
    cat("\nconsider converting it to a matrix before attempting to make an rrpp.data.frame.")
    class(dots) <- "rrpp.data.frame"
  } else if(length(dots) == 1 && inherits(dots[[1]], "rrpp.data.frame")) {
    dots <- dots[[1]]
  } else {
    if(length(dots) > 1 && inherits(dots[[1]], "rrpp.data.frame")) {
      dots1 <- dots[[1]]
      dots2 <- dots[-1]
      dots <- c(dots1, dots2)
    }
    N <- length(dots)
    dots.ns <- array(NA,N)
    for(i in 1:N){
      if(is.array(dots[[i]])) {
        if(length(dim(dots[[i]])) == 3) dots.ns[i] <- dim(dots[[i]])[[3]]
        if(length(dim(dots[[i]])) == 2) dots.ns[i] <- dim(dots[[i]])[[2]]
        if(length(dim(dots[[i]])) == 1) dots.ns[i] <- dim(dots[[i]])[[1]]
      }
      if(is.matrix(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[1]]
      if(class(dots[[i]]) == "dist") dots.ns[i] <- attr(dots[[i]], "Size")
      if(is.data.frame(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[1]]
      if(is.vector(dots[[i]])) dots.ns[i] <- length(dots[[i]])
      if(is.factor(dots[[i]])) dots.ns[i] <- length(dots[[i]])
      if(is.logical(dots[[i]])) dots.ns[i] <- length(dots[[i]])
    }
    if(any(is.na(dots.ns))) stop("Some input is either dimensionless or inappropriate for data frames")
    if(length(unique(dots.ns)) > 1) stop("Inputs have different numbers of observations")
    class(dots) <- c("rrpp.data.frame")
  }
  dots
}

#####----------------------------------------------------------------------------------------------------

# SUPPORT FUNCTIONS


# makeDf
# makes any list into a data.frame/model.frame object
# used in lm.rrpp/rrpp.fit functions
makeDf <- function(Terms, Y, data = NULL) {
  if(length(attr(Terms, "term.labels")) == 0) {
    
    df <- list()
    df$Y <- as.matrix(Y)
    df <- as.data.frame(df)
    
  } else {
    
    ind.var.names <- rownames(attr(Terms, "factors"))[-1]
    df <- vector("list", length = length(ind.var.names))
    names(df) <- ind.var.names
    for(i in 1:length(ind.var.names)) {
      f <- as.formula(paste("~", ind.var.names[i]))
      temp <- if(is.null(data)) try(eval(f[[2]]), silent = TRUE) else 
        try(eval(f[[2]], as.environment(data), enclos = parent.frame()), silent = TRUE)
      if(inherits(temp, "try-error")) temp <- try(eval(f[[2]]), silent = TRUE)
      if(inherits(temp, "try-error"))
        stop("Cannot find data in global environment.\n",
                                           call. = FALSE) else
                                             df[[i]] <- temp
    }
    
    if(is.null(names(df))) names(df) <- ind.var.names
    
    df <- as.data.frame(df)
    
  }
  
  df$Y <- as.matrix(Y)
  df <- droplevels(df)
  df
}
  
# rrpp.fit + subfunctions
# lm-like fit modified for all submodels
# general workhorse for all 'rrpp.lm' functions
# used in all 'rrpp.lm' functions

lm.args.from.lm <- function(f1){
  Terms <- terms(f1)
  data <- model.frame(Terms, data = f1$model)
  weights <- f1$weights
  offset <- f1$offset
  form <- formula(f1)
  
  dep <- eval(form[[2]], data, parent.frame())
  if(inherits(dep, "dist")) {
    if(any(dep < 0)) stop("Distances in distance matrix cannot be less than 0")
    D <- dep
  } else if((is.matrix(dep) || is.data.frame(dep))
            && isSymmetric(unname(as.matrix(dep)))) {
    D <- as.dist(dep)
  } else D <- NULL
  
  if(!is.null(D)) Y <- pcoa(D) else Y <- as.matrix(dep)
  data$Y <- Y
  form <- update(form, Y ~ .)
  Terms <- terms(form)
  data <- model.frame(Terms, data = data)
  
  list(Terms = Terms, data = data,  weights = weights, offset = offset, 
       Y = Y, D = D)
}

lm.args.from.formula <- function(f1, data = NULL){
  Terms <- terms(f1)
  dep <- eval(f1[[2]], data, parent.frame())
  if(is.array(dep) && length(dim(dep)) > 2)
    stop("Data are arranged in an array rather than a matrix.  Please update data first.\n", 
         call. = FALSE)
  if(inherits(dep, "dist")) {
    if(any(dep < 0)) stop("Distances in distance matrix cannot be less than 0")
    D <- dep
  } else if((is.matrix(dep) || is.data.frame(dep))
            && isSymmetric(unname(as.matrix(dep)))) {
    D <- as.dist(dep)
  } else D <- NULL
  
  if(!is.null(D)) Y <- pcoa(D) else Y <- as.matrix(dep)
  data <- try(makeDf(Terms, Y, data), silent = TRUE)
  if(inherits(data, "try-error"))
    stop("It was not possible to find model terms in the global environment or the data frame used.\n",
         call. = FALSE)
  form <- update(f1, Y ~ .)
  Terms <- terms(form)
  list(Terms = Terms, data = data,  weights = NULL, offset = NULL, 
       Y = Y, D = D)
}

rrpp.fit <- function(Terms, data, weights = weights, offset = offset,
                     Y, D, 
                     keep.order, pca,
                     SS.type = "I"){
  
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","II", "III")))) SS.type <- "I"
  
  dims.Y <- dim(Y)
  n <- dims.Y[1]
  p <- dims.Y[2]
  Terms <- terms(Terms, keep.order = keep.order)
  form <- formula(Terms)
  form <- update(form, Y ~ .)
  
  if(pca) data$Y <- prcomp(Y, tol = sqrt(.Machine$double.eps))$x 
  
  pdf.args <- list(data=data,
                   x = model.matrix(Terms, data = data),
                   w = weights,
                   offset = offset,
                   terms = Terms,
                   formula = form,
                   SS.type = SS.type,
                   D = D)
  
  if(is.null(pdf.args$w))
    pdf.args$w <- rep(1, NROW(pdf.args$data))
  if(is.null(pdf.args$offset))
    pdf.args$offset <- rep(0, NROW(pdf.args$data))
  if(sum(attr(pdf.args$x, "assign")) == 0)
    out <- rrpp.fit.int(pdf.args) else
      out <- rrpp.fit.lm(pdf.args)
  out
}

# rrpp.lm.fit
# helper for rrpp.fit, resembles lm.fit

rrpp.lm.fit <- function(x, y, o = NULL){
  k <- ncol(x)
  p <- NCOL(y)
  if(!is.null(o)) y <- y - o
  n <- nrow(x)
  Q <- qr(x)
  out <- list()
  out$coefficients <- qr.coef(Q, y)
  out$residuals <- qr.resid(Q, y)
  out$rank <- Q$rank
  out$fitted.values <- qr.fitted(Q, y, Q$rank)
  out$assign <- attr(x, "assign")
  out$qr <- Q
  out$df.residual <- n - Q$rank
  out
}


# rrpp.fit.lm
# base for rrpp.fit
rrpp.fit.lm <- function(a){
  # set-up
  w <- a$w
  if(any(w <= 0)) stop("Weights must be positive")
  o <-a$offset
  dat <- a$data
  Y <- dat$Y
  D <- a$D
  X <- a$x
  n <- NROW(Y)
  SS.type <- a$SS.type
  dat <- a$data
  wY <- Y*sqrt(w); wX <- X*sqrt(w)
  # data and design matrix
  Terms <- a$terms
  X.k <- X.k.obs <- attr(X, "assign")
  X.n.k.obs <- length(X.k.obs)
  QRx <- qr(X)
  X.n.k <- QRx$rank
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0, X.k))
  
  if(X.n.k < X.n.k.obs) fix <- TRUE else fix <- FALSE
  
  if(fix) {
    Terms <- Terms[uk]
    cat("\nWarning: Because variables in the linear model are redundant,")
    cat("\nthe linear model design has been truncated (via QR decomposition).")
    cat("\nOriginal X columns:", X.n.k.obs)
    cat("\nFinal X columns (rank):", X.n.k)
    cat("\nCheck coefficients or degrees of freedom in ANOVA to see changes.\n\n")
  } 
  
  term.labels <- attr(Terms, "term.labels")
  k <- length(term.labels)


  if(SS.type == "III"){
    uk0 <- uk[-1]
    xk0 <- unique(X.k[-1])
    terms.r <- lapply(1:k, function(j) Terms[xk0 %in% uk0[-j]])
    terms.f <- lapply(1:k, function(.) Terms)
    Xrs <- lapply(2:length(uk), function(j)  X[, X.k %in% uk[-j]])
    Xfs <- lapply(2:length(uk), function(j)  X)
  }
  
  if(SS.type == "II") {
    uk0 <- uk[-1]
    xk0 <- unique(X.k[-1])
    fac <- crossprod(attr(Terms, "factor"))
    terms.r <- lapply(1:k, function(j){
      ind <- ifelse(fac[j,] < fac[j,j], 1, 0)
      ind <- as.logical(ind)
      Terms[xk0 %in% uk0[ind]]
    })
    terms.f <- lapply(1:k, function(j){
      ind <- ifelse(fac[j,] < fac[j,j], 1, 0)
      ind[j] <- 1
      ind <- as.logical(c(ind))
      Terms[xk0 %in% uk0[ind]]
    })
    Xrs <- lapply(1:NROW(fac), function(j){
      ind <- ifelse(fac[j,] < fac[j,j], 1, 0)
      ind <- as.logical(c(1,ind))
      X[, X.k %in% uk[ind]]
    })
    Xfs <- lapply(1:NROW(fac), function(j){
      ind <- ifelse(fac[j,] < fac[j,j], 1, 0)
      ind[j] <- 1
      ind <- as.logical(c(1,ind))
      X[, X.k %in% uk[ind]]
    })
  }
  
  if(SS.type == "I") {
    terms.f <- lapply(1:k, function(j) Terms[1:j])
    terms.r <- lapply(1:k, function(j) Terms[0:(j - 1)])
    Xs <- lapply(1:length(uk), function(j)  Xj <- X[, X.k %in% uk[1:j]])
    Xrs <- Xs[1:k]
    Xfs <- Xs[2:(k+1)]
  }
  
  names(terms.r) <- names(terms.f) <- names(Xrs) <- names(Xfs) <- term.labels
  
  # unweighted output
  QRs.reduced <- lapply(Xrs, function(x) qr(x))
  fits.reduced <- lapply(Xrs, function(x) rrpp.lm.fit(as.matrix(x),Y, o = o))
  fitted.reduced <- lapply(fits.reduced, function(x) as.matrix(x$fitted.values))
  residuals.reduced <- lapply(fits.reduced, function(x) as.matrix(x$residuals))
  coefficients.reduced <- lapply(fits.reduced, function(x) as.matrix(x$coefficients))
  QRs.full <- lapply(Xfs, function(x) qr(x))
  fits.full <- lapply(Xfs, function(x) rrpp.lm.fit(as.matrix(x),Y, o = o))
  fitted.full <- lapply(fits.full, function(x) as.matrix(x$fitted.values))
  residuals.full <- lapply(fits.full, function(x) as.matrix(x$residuals))
  coefficients.full <- lapply(fits.full, function(x) as.matrix(x$coefficients))
  # weighted output
  if(sum(w) == n) {
    wXrs <- Xrs
    wQRs.reduced <- QRs.reduced
    wFitted.reduced <- fitted.reduced
    wResiduals.reduced <- residuals.reduced
    wCoefficients.reduced <- coefficients.reduced
    wXfs <- Xfs
    wQRs.full <- QRs.full
    wFitted.full <- fitted.full
    wResiduals.full <- residuals.full
    wCoefficients.full <- coefficients.full
  } else{
    wXrs <- lapply(Xrs, function(x) x*sqrt(w))
    wQRs.reduced <- lapply(wXrs, function(x) qr(x))
    wfits.reduced <- lapply(Xrs, function(x) lm.wfit(as.matrix(x),Y, w, offset = o))
    wFitted.reduced <- lapply(wfits.reduced, function(x) as.matrix(x$fitted.values))
    wResiduals.reduced <- lapply(wfits.reduced, function(x) as.matrix(x$residuals))
    wCoefficients.reduced <- lapply(wfits.reduced, function(x) as.matrix(x$coefficients))
    wXfs <- lapply(Xfs, function(x) x*sqrt(w))
    wQRs.full <- lapply(wXfs, function(x) qr(x))
    wfits.full <- lapply(Xfs, function(x) lm.wfit(as.matrix(x),Y, w, offset = o))
    wFitted.full<- lapply(wfits.full, function(x) as.matrix(x$fitted.values))
    wResiduals.full <- lapply(wfits.full, function(x) as.matrix(x$residuals))
    wCoefficients.full <- lapply(wfits.full, function(x) as.matrix(x$coefficients))
  }
  # additional output
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  
  out <- list(Y=Y, wY=wY, X=X, Xrs=Xrs, Xfs=Xfs,
              wX=wX, wXrs=wXrs, wXfs=wXfs,
              QRs.reduced = QRs.reduced,
              QRs.full = QRs.full,
              wQRs.reduced = wQRs.reduced,
              wQRs.full = wQRs.full,
              fitted.reduced = fitted.reduced,
              fitted.full = fitted.full,
              wFitted.reduced =wFitted.reduced,
              wFitted.full = wFitted.full,
              residuals.reduced = residuals.reduced,
              residuals.full = residuals.full,
              wResiduals.reduced = wResiduals.reduced,
              wResiduals.full = wResiduals.full,
              coefficients.reduced = coefficients.reduced,
              coefficients.full = coefficients.full,
              wCoefficients.reduced = wCoefficients.reduced,
              wCoefficients.full = wCoefficients.full,
              weights = w, offset = o, data = dat,
              SS.type = SS.type, D = D,
              Terms = Terms, term.labels = term.labels,
              model.sets = model.sets)
  class(out) <- "rrpp.fit"
  invisible(out)
}

# rrpp.fit.int
# base for rrpp.fit
# same as rrpp.fit.lm, but special case where only an intercept is found
rrpp.fit.int <- function(a) {
  # set-up
  w <- a$w
  if(any(w <= 0)) stop("Weights must be positive")
  o <-a$offset
  dat <- a$data
  Y <- dat$Y
  D <- a$D
  X <- a$x
  n <- NROW(Y)
  SS.type <- a$SS.type
  dat <- a$data
  wY <- Y*sqrt(w); wX <- X*sqrt(w)
  Terms <- a$terms
  # unweighted output
  Xrs <- NULL
  QRs.reduced <- NULL
  fits.reduced <- NULL
  fitted.reduced <- NULL
  residuals.reduced <- NULL
  coefficients.reduced <- NULL
  Xfs <- list(X)
  QRs.full <- list(qr(X))
  fits.full <- rrpp.lm.fit(as.matrix(X),Y, o = o)
  fitted.full <- list(as.matrix(fits.full$fitted.values))
  residuals.full <- list(as.matrix(fits.full$residuals))
  coefficients.full <- list(as.matrix(fits.full$coefficients))
  # weighted output
  if(sum(w) == n){
    wQRs.reduced <- QRs.reduced
    wFitted.reduced <- fitted.reduced
    wResiduals.reduced <- residuals.reduced
    wCoefficients.reduced <- coefficients.reduced
    wXrs <- Xrs
    wFitted.reduced <- fitted.reduced
    wResiduals.reduced <- residuals.reduced
    wCoefficients.reduced <- coefficients.reduced
    wXfs <- Xfs
    wQRs.full <- QRs.full
    wFitted.full <- fitted.full
    wResiduals.full <- residuals.full
    wCoefficients.full <- coefficients.full
  } else{
    wXrs <- NULL
    wQRs.reduced <- NULL
    wFitted.reduced <- NULL
    wResiduals.reduced <- NULL
    wCoefficients.reduced <- NULL
    wXfs <- list(X*sqrt(w))
    wQRs.full <- list(qr(X*sqrt(w)))
    wfits.full <- lm.wfit(as.matrix(X),Y, w = w, offset = o)
    wFitted.full<- list(as.matrix(wfits.full$fitted.values))
    wResiduals.full <- list(as.matrix(wfits.full$residuals))
    wCoefficients.full <- list(as.matrix(wfits.full$coefficients))
  }
  term.labels <- attr(Terms, "term.labels")
  model.sets <- NULL
  out <- list(Y=Y, wY=wY, X=X, Xrs=Xrs, Xfs=Xfs,
              wX=wX, wXrs=wXrs, wXfs=wXfs,
              QRs.reduced = QRs.reduced,
              QRs.full = QRs.full,
              wQRs.reduced = wQRs.reduced,
              wQRs.full = wQRs.full,
              fitted.reduced = fitted.reduced,
              fitted.full = fitted.full,
              wFitted.reduced =wFitted.reduced,
              wFitted.full =wFitted.full,
              residuals.reduced = residuals.reduced,
              residuals.full = residuals.full,
              wResiduals.reduced = wResiduals.reduced,
              wResiduals.full = wResiduals.full,
              coefficients.reduced = coefficients.reduced,
              coefficients.full = coefficients.full,
              wCoefficients.reduced = wCoefficients.reduced,
              wCoefficients.full = wCoefficients.full,
              weights = w, offset = o, data = dat,
              SS.type = NULL, D = D,
              Terms = Terms, term.labels = term.labels,
              model.sets = model.sets)
  class(out) <- "procD.fit"
  invisible(out)
}

# rrpp.fit
# calls one of previous functions, depending on conditions


# rrpp + subfunctions
# The principal function for RRPP
# used in SS.iter

# set-up
rrpp.setup <- function(fit, ind){
  fitted <- fit$fitted.reduced
  residuals <- fit$residuals.reduced
  n <- NROW(fit$Y)
  perms <- length(ind)
  if(sum(fit$weights) == n) w <- NULL else w <- sqrt(fit$weights)
  if(sum(fit$offset) == 0) o <- NULL else o <- fit$offset
  lapply(1:perms, function(j){
    list(fitted = fitted, residuals = residuals,
         w = w, o = o, ind.i = ind[[j]])
  })
}

# RRPP for OLS
rrpp.basic <- function(fitted, residuals, ind.i){
  Map(function(f, r) f+r[ind.i,], fitted, residuals)
}

# RRPP for wLS
rrpp.w <- function(fitted, residuals, ind.i, w){
  Map(function(f, r) (f+r[ind.i,])*w, fitted, residuals)
}

# if there is an offset
rrpp.o <- function(fitted, residuals, ind.i, o){
  Map(function(f, r) (f+r[ind.i,]) - o, fitted, residuals)
}

# weighted LS with offset
rrpp.w.o <- function(fitted, residuals, ind.i, w, o){
  Map(function(f, r) (f+r[ind.i,])*w - o, fitted, residuals)
}

# calls other functions as needed
rrpp <- function(fitted, residuals, ind.i, w, o){
  if(!is.null(w) && !is.null(o)) r <- rrpp.w.o(fitted, residuals, ind.i, w, o)
  if(!is.null(w) && is.null(o)) r <- rrpp.w(fitted, residuals, ind.i, w)
  if(is.null(w) && !is.null(o)) r <- rrpp.o(fitted, residuals, ind.i, o)
  if(is.null(w) && is.null(o)) r <- rrpp.basic(fitted, residuals, ind.i)
  r
}


# droplevels.rrpp.data.frame
# same as droplevels for data.frame objects
# used in lm.rrpp

droplevels.rrpp.data.frame <- function (x, except = NULL, ...) {
  ix <- vapply(x, is.factor, NA)
  if (!is.null(except))
    ix[except] <- FALSE
  x[ix] <- lapply(x[ix], factor)
  x
}

# SS.iter
# workhorse for lm.rrpp
# used in lm.rrpp

SS.iter <- function(fit, ind, P = NULL, RRPP = TRUE, print.progress = TRUE) {
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  fitted <- fit$wFitted.reduced
  res <- fit$wResiduals.reduced
  Y <- as.matrix(fit$wY)
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  trms <- fit$term.labels
  k <- length(trms)
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  rrpp.args <- list(fitted = fitted, residuals = res,
                    ind.i = NULL, w = NULL, o = NULL)
  if(offset) rrpp.args$o <- o
  if(print.progress){
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  if(gls){
    Y <- crossprod(P, Y)
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Qr <- lapply(Xr, qr)
    Qf <- lapply(Xf, qr)
    Ur <- lapply(Qr, function(x) qr.Q(x))
    Uf <- lapply(Qf, function(x) qr.Q(x))
    Ufull <- Uf[[k]]
    int <- attr(fit$Terms, "intercept")
    Unull <- qr.Q(qr(crossprod(P, rep(int, n))))
    yh0 <- fastFit(Unull, Y, n, p)
    r0 <- Y - yh0
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      fitted <- Map(function(u) crossprod(tcrossprod(u), Y), Ur)
      res <- lapply(fitted, function(f) Y - f)
    }
    rrpp.args$fitted <- fitted
    rrpp.args$residuals <- res
    result <- lapply(1: perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- yh0 + r0[x,]
      yy <- sum(y^2)
      SSr <- sapply(1:k, function(j) sum(crossprod(Ur[[j]], Yi[[j]])^2))
      SSf <- sapply(1:k, function(j) sum(crossprod(Uf[[j]], Yi[[j]])^2))
      RSS <- sapply(1:k, function(j) sum(Yi[[j]]^2) - sum(crossprod(Ufull, Yi[[j]])^2))
      TSS <- yy - sum(crossprod(Unull, y)^2)
      TSS <- rep(TSS, k)
      RSS.model <- yy - sum(crossprod(Ufull, y)^2)
      list(SS = SSf - SSr, RSS=RSS, TSS=TSS, RSS.model = RSS.model)
    })
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
      rrpp.args$fitted <- fitted
      rrpp.args$residuals <- res
    }
    Ur <- lapply(fit$wQRs.reduced, function(x) qr.Q(x))
    Uf <- lapply(fit$wQRs.full, function(x) qr.Q(x))
    Ufull <- Uf[[k]]
    int <- attr(fit$Terms, "intercept")
    Unull <- qr.Q(qr(rep(int, n)))
    yh0 <- fastFit(Unull, Y, n, p)
    r0 <- Y - yh0
    result <- lapply(1: perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- yh0 + r0[x,]
      yy <- sum(y^2)
      SSr <- sapply(1:k, function(j) sum(crossprod(Ur[[j]], Yi[[j]])^2))
      SSf <- sapply(1:k, function(j) sum(crossprod(Uf[[j]], Yi[[j]])^2))
      RSS <- sapply(1:k, function(j) sum(Yi[[j]]^2) - sum(crossprod(Ufull, Yi[[j]])^2))
      TSS <- yy - sum(crossprod(Unull, y)^2)
      TSS <- rep(TSS, k)
      RSS.model <- yy - sum(crossprod(Ufull, y)^2)
      list(SS = SSf - SSr, RSS=RSS, TSS=TSS, RSS.model = RSS.model)
    })
  }
  SS <- sapply(1:perms, function(j) result[[j]]$SS)
  RSS <- sapply(1:perms, function(j) result[[j]]$RSS)
  TSS <- sapply(1:perms, function(j) result[[j]]$TSS)
  RSS.model <- sapply(1:perms, function(j) result[[j]]$RSS.model)
  if(k == 1) {
    SS <- matrix(SS, 1, length(SS))
    RSS <- matrix(RSS, 1, length(RSS))
    TSS <- matrix(TSS, 1, length(TSS))
  }
  rownames(SS) <- rownames(RSS) <- rownames(TSS) <- trms
  colnames(SS) <- colnames(RSS) <- colnames(TSS) <- names(RSS.model) <- 
    c("obs", paste("iter", 1:(perms-1), sep="."))
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
}

# SS.iterPP
# workhorse for lm.rrpp, same as SS.iter, but with parallel processing
# used in lm.rrpp

SS.iterPP <- function(fit, ind, P = NULL, RRPP = TRUE, print.progress = TRUE) {
  cl <- detectCores()-1
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  if(print.progress){
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    cat(paste("Progress bar not possible with parallel processing, but this shouldn't take long...\n"))
  }
  fitted <- fit$wFitted.reduced
  res <- fit$wResiduals.reduced
  Y <- as.matrix(fit$wY)
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  trms <- fit$term.labels
  k <- length(trms)
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  rrpp.args <- list(fitted = fitted, residuals = res,
                    ind.i = NULL, w = NULL, o = NULL)
  if(offset) rrpp.args$o <- o
  if(gls){
    Y <- crossprod(P, Y)
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Ur <- lapply(Xr, function(x) qr.Q(qr(x)))
    Uf <- lapply(Xf, function(x) qr.Q(qr(x)))
    Ufull <- Uf[[k]]
    int <- attr(fit$Terms, "intercept")
    Unull <- qr.Q(qr(crossprod(P, rep(int, n))))
    yh0 <- fastFit(Unull, Y, n, p)
    r0 <- Y - yh0
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      fitted <- Map(function(u) crossprod(tcrossprod(u), Y), Ur)
      res <- lapply(fitted, function(f) Y - f)
    }
    rrpp.args$fitted <- fitted
    rrpp.args$residuals <- res
    result <- mclapply(1: perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- yh0 + r0[x,]
      yy <- sum(y^2)
      SSr <- sapply(1:k, function(j) sum(crossprod(Ur[[j]], Yi[[j]])^2))
      SSf <- sapply(1:k, function(j) sum(crossprod(Uf[[j]], Yi[[j]])^2))
      RSS <- sapply(1:k, function(j) sum(Yi[[j]]^2) - sum(crossprod(Ufull, Yi[[j]])^2))
      TSS <- yy - sum(crossprod(Unull, y)^2)
      TSS <- rep(TSS, k)
      RSS.model <- yy - sum(crossprod(Ufull, y)^2)
      list(SS = SSf - SSr, RSS=RSS, TSS=TSS, RSS.model = RSS.model)
    }, mc.cores = cl)
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
      rrpp.args$fitted <- fitted
      rrpp.args$residuals <- res
    }
    Ur <- lapply(fit$wQRs.reduced, function(x) qr.Q(x))
    Uf <- lapply(fit$wQRs.full, function(x) qr.Q(x))
    Ufull <- Uf[[k]]
    int <- attr(fit$Terms, "intercept")
    Unull <- qr.Q(qr(rep(int, n)))
    yh0 <- fastFit(Unull, Y, n, p)
    r0 <- Y - yh0
    result <- mclapply(1: perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- yh0 + r0[x,]
      yy <- sum(y^2)
      SSr <- sapply(1:k, function(j) sum(crossprod(Ur[[j]], Yi[[j]])^2))
      SSf <- sapply(1:k, function(j) sum(crossprod(Uf[[j]], Yi[[j]])^2))
      RSS <- sapply(1:k, function(j) sum(Yi[[j]]^2) - sum(crossprod(Ufull, Yi[[j]])^2))
      TSS <- yy - sum(crossprod(Unull, y)^2)
      TSS <- rep(TSS, k)
      RSS.model <- yy - sum(crossprod(Ufull, y)^2)
      list(SS = SSf - SSr, RSS=RSS, TSS=TSS, RSS.model = RSS.model)
    }, mc.cores = cl)
  }
  SS <- sapply(1:perms, function(j) result[[j]]$SS)
  RSS <- sapply(1:perms, function(j) result[[j]]$RSS)
  TSS <- sapply(1:perms, function(j) result[[j]]$TSS)
  RSS.model <- sapply(1:perms, function(j) result[[j]]$RSS.model)
  if(k == 1) {
    SS <- matrix(SS, 1, length(SS))
    RSS <- matrix(SS, 1, length(RSS))
    TSS <- matrix(SS, 1, length(TSS))
  }
  rownames(SS) <- rownames(RSS) <- rownames(TSS) <- trms
  colnames(SS) <- colnames(RSS) <- colnames(TSS) <- names(RSS.model) <- 
    c("obs", paste("iter", 1:(perms-1), sep="."))

  list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
}

SS.iter.null <- function(fit, ind, P = NULL, RRPP=TRUE, print.progress = TRUE) {
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  if(print.progress){
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  fitted <- fit$wFitted.full
  res <- fit$wResiduals.full
  Y <- as.matrix(fit$wY)
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  k <- 1
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(gls){
    Y <- crossprod(P, Y)
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      int <- attr(fit$Terms, "intercept")
      U <- qr.Q(qr(crossprod(P, rep(int, n))))
      fitted <- crossprod(tcrossprod(U), Y)
      res <- Y - fitted
    }
    SS <- lapply(1:perms, function(j){
      x <-ind[[j]]
      y <- fitted + res[x,]; yy <- sum(y^2)
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      yy - sum(crossprod(U, y)^2)
    })
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    }
    int <- attr(fit$Terms, "intercept")
    U <- qr.Q(qr(rep(int, n)))
    SS <- lapply(1:perms, function(j){
      x <-ind[[j]]
      y <- Y[x,]; yy <- sum(y^2)
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      yy - sum(crossprod(U, y)^2)
    })
  }
  SS <- matrix(unlist(SS), 1, perms)
  rownames(SS) <- "Residuals"
  colnames(SS) <- c("obs", paste("iter", 1:(perms-1), sep=":"))
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  SS
}

# anova.parts
# construct an ANOVA tablefrom random SS output
# used in lm.rrpp

anova.parts <- function(fit, SS){
  SS.type <- fit$SS.type
  perms <- NCOL(SS)
  QRf <- fit$wQRs.full
  QRr <- fit$wQRs.reduced
  k <- length(QRf)
  dims <- dim(as.matrix(fit$Y))
  n <- dims[1]; p <- dims[2]
  df <- unlist(Map(function(qf, qr) qf$rank - qr$rank, QRf, QRr))
  dfe <- n - QRf[[k]]$rank
  
  RSS <- SS$RSS
  TSS <- SS$TSS
  RSS.model <- SS$RSS.model
  SS <- SS$SS
  Rsq <- SS/TSS
  MS <- SS/df
  RMS <- RSS/dfe
  Fs <- MS/RMS
  dft <- sum(df, dfe)
  df <- c(df, dfe, dft)
  if(SS.type == "III") {
    etas <- SS/TSS
    cohenf <- etas/(1-etas)
  } else {
    etas <- Rsq
    if(k == 1) unexp <- 1 - etas else unexp <- 1 - apply(etas, 2, cumsum)
    cohenf <- etas/unexp
  }
  rownames(Fs) <- rownames(cohenf) <- rownames(SS)
  out <- list(SS.type = SS.type, SS = SS, MS = MS, RSS = RSS,
              TSS = TSS, RSS.model = RSS.model, Rsq = Rsq,
              Fs = Fs, cohenf = cohenf,
              n = n, p = p, df=df
  )
  out
}

# SS.mean
# support function for calculating SS quickly in SS.iter
# used in lm.rrpp

SS.mean <- function(x, n) if(is.vector(x)) sum(x)^2/n else sum(colSums(x)^2)/n

# generate U matrices
# makes a U matrix based on model design
# not used bt was used in compare.models, dureing development.  Retained
# in case it is useful down the road.

cm.U <- function(x){
  X <- x$LM$X
  w <- sqrt(x$LM$weights)
  if(x$LM$ols) out <- qr.Q(qr(X*w))
  if(x$LM$gls){
    P <- x$LM$Pcov
    out <- crossprod(P, qr.Q(qr(crossprod(P,X*w))))
  }
  out
}

# beta.iter
# gets appropriate beta vectors for random permutations in lm.rrpp
# generates distances as statistics for summary

beta.iter <- function(fit, ind, P = NULL, RRPP = TRUE, print.progress = TRUE) {
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  fitted <- fit$wFitted.reduced
  res <- fit$wResiduals.reduced
  Y <- fit$wY
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  trms <- fit$term.labels
  k <- length(trms)
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  rrpp.args <- list(fitted = fitted, residuals = res,
                    ind.i = NULL, w = NULL, o = NULL)
  if(offset) rrpp.args$o <- o
  if(print.progress){
    cat(paste("\nCoefficients estimation:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  if(gls){
    Y <- crossprod(P, Y)
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Qr <- lapply(Xr, qr)
    Qf <- lapply(Xf, qr)
    Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
    Ur <- lapply(Qr, function(x) qr.Q(x))
  
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      fitted <- Map(function(u) crossprod(tcrossprod(u), Y), Ur)
      res <- lapply(fitted, function(f) Y - f)
    }
    rrpp.args$fitted <- fitted
    rrpp.args$residuals <- res
    betas <- lapply(1:perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      Map(function(h,y) h %*% y, Hf, Yi)
    })
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
      rrpp.args$fitted <- fitted
      rrpp.args$residuals <- res
    }
    Qf <- lapply(fit$wXfs, qr)
    Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
    
    betas <- lapply(1:perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      Map(function(h,y) h%*%y, Hf, Yi)
    })
  }
  beta.mats <- lapply(1:k, function(j){
    result <- lapply(betas, function(x) x[[j]])
    result
  })
  beta.mat.d <- lapply(1:k, function(j){
    b <- beta.mats[[j]]
    kk <- length(b)
    result <- sapply(1:kk, function(jj){
      bb <- b[[jj]]
      if(p == 1) res <- abs(bb) else res <- sqrt(diag(tcrossprod(bb)))
    })
    rownames(result) <- rownames(b[[1]])
    colnames(result) <- c("obs", paste("iter", seq(1,(perms-1),1), sep = "."))
    result
  })
  
  names(beta.mats) <- names(beta.mat.d) <- trms
  
  if(k == 1) d.stitched <- beta.mat.d[[1]]
  
  if(k > 1) {
    d.stitched <- beta.mat.d[[k]]
    
    beta.match <- lapply(1:k, function(j){
      cr <- colnames(fit$wQRs.reduced[[j]]$qr)
      cf <- colnames(fit$wQRs.full[[j]]$qr)
      cf[!(cf %in% cr)]
    })
    
    for(i in 1:(k-1)){
      beta.check <- beta.match[[i]]
      d.check <- beta.mat.d[[i]]
      target <- which(rownames(d.stitched) %in% beta.check)
      d.stitched[target, ] <- d.check[rownames(d.check) %in% beta.check,]
    }
  }
  
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  out <- list(random.coef = beta.mats,
              random.coef.distances = d.stitched)
  out
}

beta.iterPP <- function(fit, ind, P = NULL, RRPP = TRUE, print.progress = TRUE) {
  cl <- detectCores()-1
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  if(print.progress){
    cat(paste("\nCoefficients estimation:", perms, "permutations.\n"))
    cat(paste("Progress bar not possible with parallel processing, but this shouldn't take long...\n"))
  }
  perms <- length(ind)
  fitted <- fit$wFitted.reduced
  res <- fit$wResiduals.reduced
  Y <- fit$wY
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  trms <- fit$term.labels
  k <- length(trms)
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  rrpp.args <- list(fitted = fitted, residuals = res,
                    ind.i = NULL, w = NULL, o = NULL)
  if(offset) rrpp.args$o <- o
  if(gls){
    Y <- crossprod(P, Y)
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Qr <- lapply(Xr, qr)
    Qf <- lapply(Xf, qr)
    Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
    Ur <- lapply(Qr, function(x) qr.Q(x))
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      fitted <- Map(function(u) crossprod(tcrossprod(u), Y), Ur)
      res <- lapply(fitted, function(f) Y - f)
    }
    rrpp.args$fitted <- fitted
    rrpp.args$residuals <- res
    betas <- mclapply(1:perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      Map(function(h,y) h %*% y, Hf, Yi)
    }, mc.cores = cl)
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
      rrpp.args$fitted <- fitted
      rrpp.args$residuals <- res
    }
    Qf <- lapply(fit$wXfs, qr)
    Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
    betas <- mclapply(1:perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      Map(function(h,y) h%*%y, Hf, Yi)
    }, mc.cores = cl)
  }
  beta.mats <- lapply(1:k, function(j){
    result <- lapply(betas, function(x) x[[j]])
    result
  })
  beta.mat.d <- lapply(1:k, function(j){
    b <- beta.mats[[j]]
    kk <- length(b)
    result <- sapply(1:kk, function(jj){
      bb <- b[[jj]]
      if(p == 1) res <- abs(bb) else res <- sqrt(diag(tcrossprod(bb)))
    })
    rownames(result) <- rownames(b[[1]])
    colnames(result) <- c("obs", paste("iter", seq(1,(perms-1),1), sep = "."))
    result
  })
  
  names(beta.mats) <- names(beta.mat.d) <- trms
  
  if(k == 1) d.stitched <- beta.mat.d[[1]]
  
  if(k > 1) {
    d.stitched <- beta.mat.d[[k]]
    
    beta.match <- lapply(1:k, function(j){
      cr <- colnames(fit$wQRs.reduced[[j]]$qr)
      cf <- colnames(fit$wQRs.full[[j]]$qr)
      cf[!(cf %in% cr)]
    })
    
    for(i in 1:(k-1)){
      beta.check <- beta.match[[i]]
      d.check <- beta.mat.d[[i]]
      target <- which(rownames(d.stitched) %in% beta.check)
      d.stitched[target, ] <- d.check[rownames(d.check) %in% beta.check,]
    }
  }
  out <- list(random.coef = beta.mats,
              random.coef.distances = d.stitched)
  out
}


beta.iter.null <- function(fit, ind, P = NULL, RRPP = TRUE, print.progress = TRUE) {
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  fitted <- fit$wFitted.full
  res <- fit$wResiduals.full
  Y <- fit$wY
  X <- fit$wX
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  k <- 1
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  if(print.progress){
    cat(paste("\nCoefficients estimation:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  if(gls){
    Y <- crossprod(P, Y)
    X <- crossprod(P, X)
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      U <- qr.Q(qr(crossprod(P, matrix(1, n))))
      fitted <- crossprod(tcrossprod(U), Y)
      res <- lapply(fitted, function(f) Y - f)
    }
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    }
  }    
  Q <- qr(X)
  H <- tcrossprod(solve(qr.R(Q)), qr.Q(Q))
  betas <- lapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    x <-ind[[j]]
    y <- fitted[[1]] + res[[1]][x,]
    H %*% y
  })
  coef.d <- function(b) sqrt(sum(b^2))
  beta.d <- sapply(betas, coef.d)
  betas <- lapply(1:perms, function(j){
    names(betas[[j]]) <- "Intercept"
  })
  iter.names <- c("obs", paste("iter", seq(1,(perms-1),1), sep = "."))
  names(betas) <- iter.names
  names(beta.d) <- iter.names
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  out <- list(random.coef = betas,
              random.coef.distances = beta.d)
  out
}

# beta.boot
# for calculating confidence intervals for predicted values
# used in predcit.lm.rrpp

beta.boot <- function(f, r, h, ind.i){
  y <- f + r[ind.i,]
  h %*% y
}

# ellipse.points
# A helper function for plotting elipses from non-parametric CI data
# Used in predict.lm.rrpp

ellipse.points <- function(m, pr, conf) {
  m <- as.matrix(m)
  p <- NCOL(m)
  z <- qnorm((1 - conf)/2, lower.tail = FALSE)
  angles <- seq(0, 2*pi, length.out=200)
  ell    <- z* cbind(cos(angles), sin(angles))
  del <- lapply(pr, function(x) x[,1:p] - m)
  vcv <- lapply(1:NROW(m), function(j){
    x <- sapply(1:length(del), function(jj){
      del[[jj]][j,]
    })
    tcrossprod(x)/ncol(x)
  })
  R <- lapply(vcv, chol)
  ellP <- lapply(R, function(x) ell %*% x)
  np <- NROW(ellP[[1]])
  ellP <- lapply(1:length(ellP), function(j){
    x <- m[j,]
    mm <- matrix(rep(x, each = np), np, length(x))
    mm + ellP[[j]]
  })
  ellP <- simplify2array(ellP)
  pc1lim <- c(min(ellP[,1,]), max(ellP[,1,]))
  pc2lim <- c(min(ellP[,2,]), max(ellP[,2,]))
  list(ellP = ellP, pc1lim = pc1lim, pc2lim = pc2lim,
       means = m)
}

# refit
# finds rrpp.fit from lm.rrpp object
# used in Deviance
refit <- function(fit){
  dat <- fit$LM$data
  dims <- dim(dat$Y)
  pdf.args <- list(data=dat,
                   x = model.matrix(fit$LM$Terms, data = dat),
                   w = fit$LM$weights,
                   offset = fit$LM$offset,
                   terms = fit$LM$Terms,
                   formula = NULL,
                   SS.type = fit$ANOVA$SS.type)
  form <- fit$LM$form
  form <- update.formula(form, Y ~.)
  pdf.args$formula <- form
  if(sum(attr(pdf.args$x, "assign")) == 0)
    out <- rrpp.fit.int(pdf.args) else
      out <- rrpp.fit.lm(pdf.args)
  out
}

# aov.single.model
# performs ANOVA on a single model
# used in anova.lm.rrpp
aov.single.model <- function(object, ...,
                             effect.type = c("F", "cohenf", "SS", "MS", "Rsq"),
                             error = NULL) {
  x <- object$ANOVA
  df <- x$df
  k <- length(df)-2
  SS <- x$SS
  MS <- x$MS 
  RSS <-x$RSS
  TSS <- x$TSS
  perms <- object$PermInfo$perms
  pm <- object$PermInfo$perm.method
  trms <- object$LM$term.labels
  
  if(!is.null(error)) {
    if(!inherits(error, "character")) stop("The error description is illogical.  It should be a string of character values matching ANOVA terms.",
                                           call. = FALSE)
    kk <- length(error)
    if(kk != k) stop("The error description should match in length the number of ANOVA terms (not including Residuals)",
                     call. = FALSE)
    MSEmatch <- match(error, c(trms, "Residuals"))
    if(any(is.na(MSEmatch))) stop("At least one of the error terms is not an ANOVA term",
                                  call. = FALSE)
  } else MSEmatch <- NULL
  if(k >= 1) {
    Fs <- x$Fs
    
    if(!is.null(MSEmatch)){
      Fmatch <- which(MSEmatch <= k)
      Fs[Fmatch,] <- MS[Fmatch,]/MS[MSEmatch[Fmatch],]
      F.effect.adj <- apply(Fs, 1, effect.size)
    }
    
    effect.type <- match.arg(effect.type)
    
    if(object$LM$gls) {
      est <- "GLS"
      
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
    } else est <- "OLS"
    
    if(effect.type == "F") Z <- Fs
    if(effect.type == "SS") Z <- x$SS
    if(effect.type == "MS") Z <- x$MS
    if(effect.type == "Rsq") Z <- x$Rsq
    if(effect.type == "cohenf") Z <- x$cohenf
    if(effect.type == "Rsq") effect.type = "R-squared"
    if(effect.type == "cohenf") effect.type = "Cohen's f-squared"
    Fs <- Fs[,1]
    SS <- SS[,1]
    MS <- MS[,1]
    Rsq <- x$Rsq[,1]
    cohenf <- x$cohenf[,1]
    if(!is.null(Z)) {
      if(!is.matrix(Z)) Z <- matrix(Z, 1, length(Z))
      P.val <- apply(Z, 1, pval) 
      Z <- apply(log(Z), 1, effect.size)
      } else P.val <- NULL
    
    Residuals <- c(df[k+1], RSS[[1]], RSS[[1]]/df[k+1], RSS[[1]]/TSS[[1]], rep(NA, 3))
    Total <- c(df[k+2], TSS[[1]], rep(NA, 5))
    tab <- data.frame(Df=df[1:k], SS=SS, MS = MS, Rsq = Rsq, F = Fs, Z = Z, P.val = P.val)
    tab <- rbind(tab, Residuals = Residuals, Total = Total)
    colnames(tab)[NCOL(tab)] <- paste("Pr(>", effect.type, ")", sep="")
    class(tab) = c("anova", class(tab))
    SS.type <- x$SS.type
    
    out <- list(table = tab, perm.method = pm, perm.number = perms,
                est.method = est, SS.type = SS.type, effect.type = effect.type,
                call = object$call)
    
      } else {
        RSS.model <- x$RSS.model
        tab <- data.frame(df = df, SS = RSS.model[[1]], MS = RSS.model[[1]]/df)
        rownames(tab) <- c("Residuals")
        class(tab) = c("anova", class(tab))
        if(object$LM$gls) est <-"GLS" else est <- "OLS"
        out <- list(table = tab, perm.method = pm, perm.number = perms,
                    est.method = est, SS.type = NULL, effect.type = NULL,
                    call = object$call)
        
      }
  class(out) <- "anova.lm.rrpp"
  out
  }

# aov.multi.model
# performs ANOVA on multiple models
# used in anova.lm.rrpp
aov.multi.model <- function(object, lm.list,
                            effect.type = c("F", "cohenf", "SS", "MS", "Rsq"),
                            print.progress = TRUE) {
  
  effect.type <- match.arg(effect.type)
  
  if(inherits(object, "lm.rrpp")) refModel <- object else 
    stop("The reference model is not a class lm.rrpp object")
  ind <- refModel$PermInfo$perm.schedule
  perms <- length(ind)
  
  X <- as.matrix(refModel$LM$X * sqrt(refModel$LM$weights))
  B <- refModel$LM$coefficients
  Y <- refModel$LM$Y
  U <- as.matrix(qr.Q(qr(X)))
  n <- refModel$LM$n
  p <- refModel$LM$p
  if(refModel$LM$gls) {
    P <- refModel$LM$Pcov
    B <- refModel$LM$gls.coefficients
    X <- as.matrix(crossprod(P, X))
    Y <- crossprod(P, Y)
    U <- as.matrix(qr.Q(qr(X)))
  }
  Yh <- as.matrix(fastFit(U, Y, n, p))
  R <- as.matrix(Y) - Yh
  
  K <- length(lm.list)
  Ulist <- lapply(1:K, function(j){
    m <- lm.list[[j]]
    X <- m$LM$X
    w <- sqrt(m$LM$weights)
    X <- X * w
    if(m$LM$gls){
      P <- m$LM$Pcov
      X <- crossprod(P, X)
    }
    qr.Q(qr(X))
  })
  
  if(print.progress){
    if(K > 1)
    cat(paste("\nSums of Squares calculations for", K, "models:", perms, "permutations.\n")) else
      cat(paste("\nSums of Squares calculations for", K, "model:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+5, initial = 0, style=3)
  }

  int <- attr(refModel$LM$Terms, "intercept")
  if(refModel$LM$gls) {
    int <- crossprod(refModel$LM$Pcov, rep(int, n))
  } else int <- rep(int, n)
  
  U0 <- as.matrix(qr.Q(qr(int * sqrt(refModel$LM$weights))))
  yh0 <- as.matrix(fastFit(U0, Y, n, p))
  r0 <- as.matrix(Y) - yh0
  
  rY <- function(ind.i) Yh + R[ind.i,]
  rY0 <- function(ind.i) yh0 + r0[ind.i,]
  
  RSS <- function(ind.i, U, Ul, K, n, p, Y, yh0, r0) {
    y <- as.matrix(rY(ind.i))
    rss0  <- sum(y^2) - sum(crossprod(U, y)^2)
    rss <- lapply(1:K, function(j){
      u <- Ul[[j]]
      sum(y^2) - sum(crossprod(u, y)^2)
    })
    
    RSSp <- c(rss0, unlist(rss))
    
    rss <- lapply(1:K, function(j){
      y <- as.matrix(rY0(ind.i))
      u <- Ul[[j]]
      sum(y^2) - sum(crossprod(u, y)^2)
    })
    RSSy <- c(rss0, unlist(rss))
    
    c(RSSp, RSSy)
  }
  
  rss.list <- list(ind.i = NULL, U = U, 
                   Ul = Ulist, K = K, n = n , p = p, Y = Y,
                   yh0 = yh0, r0 = r0)
  
  RSSp <- sapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    rss.list$ind.i <- ind[[j]]
    do.call(RSS, rss.list)
  })
  
  RSSy <- RSSp[-(1:(K+1)),]
  RSSp <- RSSp[1:(K+1),]
  
  fit.names <- c(refModel$call[[2]], lapply(1:K, function(j) lm.list[[j]]$call[[2]]))
  rownames(RSSp) <- rownames(RSSy) <- fit.names
  
  SS <- rep(RSSp[1,], each = K + 1) - RSSp

  SSY <- sapply(1:perms, function(j){
    y <- yh0 + r0[ind[[j]],]
    sum(y^2) - sum(crossprod(U0, y)^2)
  })
  
  Rsq <-  1 - (RSSp / rep(SSY, each = K + 1))
  dfe <- n - c(object$LM$wQR$rank, unlist(lapply(1:K, 
                                                 function(j) lm.list[[j]]$LM$wQR$rank)))
  df <- dfe[1] - dfe
  df[1] <- 1
  
  MS <- SS/rep(df, perms)
  MSE <- RSSy/matrix(rep(dfe, perms), length(dfe), perms)
  Fs <- MS/MSE
  
  SS[which(zapsmall(SS) == 0)] <- 1e-32
  MS[which(zapsmall(MS) == 0)] <- 1e-32
  Rsq[which(zapsmall(Rsq) == 0)] <- 1e-32
  Fs[which(zapsmall(Fs) == 0)] <- 1e-32
  
  if(effect.type == "SS") {
    Pvals <- apply(SS, 1, pval)
    Z <- apply(log(SS), 1, effect.size)
  } else   if(effect.type == "MS") {
    Pvals <- apply(MS, 1, pval)
    Z <- apply(log(MS), 1, effect.size)
  } else   if(effect.type == "Rsq") {
    Pvals <- apply(Rsq, 1, pval)
    Z <- apply(log(Rsq), 1, effect.size)
  } else{
    Pvals <- apply(Fs, 1, pval)
    Z <- apply(log(Fs), 1, effect.size)
  }
  SS[1,] <- NA
  MS[1,] <- NA
  Fs[1,] <- NA
  Pvals[1] <- NA
  Z[1] <- NA
  
  RSS.obs <- RSSp[,1]
  SS.obs <- SS[,1]
  MS.obs <- MS[,1]
  Rsq.obs <- Rsq[,1]
  F.obs <- Fs[,1]
  
  tab <- data.frame(ResDf = dfe, Df = df, RSS = RSS.obs, SS = SS.obs, MS = MS.obs,
                    Rsq = Rsq.obs, F = F.obs, Z = Z, P = Pvals)
  tab$DF[1] <- NA
  tab$Rsq <- zapsmall(tab$Rsq, digits = 8)
  tab <-  rbind(tab, c(n-1, NA, SSY[1], NA, NA, NA, NA, NA, NA, NA))
  rownames(tab)[NROW(tab)] <- "Total"
  
  if(effect.type == "SS") p.type <- "Pr(>SS)" else
    if(effect.type == "MS") p.type <- "Pr(>MS)" else
      if(effect.type == "Rsq") p.type <- "Pr(>Rsq)" else
        if(effect.type == "cohenf") p.type <- "Pr(>cohenf)" else p.type <- "Pr(>F)" 
  names(tab)[length(names(tab))] <- p.type
  rownames(tab)[1] <- paste(rownames(tab)[1], "(Null)")
  class(tab) <- c("anova", class(tab))
  
  step <- perms + 5
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  pm <- "RRPP"
  if(refModel$LM$gls) est <- "GLS" else est <- "OLS"
  
  out <- list(table = tab, perm.method = pm, perm.number = perms,
              est.method = est, SS.type = NULL, effect.type = effect.type,
              SS = SS[-1,], MS = MS[-1,], Rsq = Rsq[-1,], F = Fs[-1,],
              call = object$call)
  
class(out) <- "anova.lm.rrpp"
out

}

# getSlopes
# gets the slopes for groups from a lm.rrpp fit
# used in pairwise
getSlopes <- function(fit, x, g){
  k <- length(fit$LM$term.labels)
  p <- fit$LM$p
  beta <- fit$LM$random.coef[[k]]
  X <- fit$LM$X
  X <- X[, colnames(X)  %in% rownames(beta[[1]])]
  getFitted <- function(b) X %*% b
  fitted <- lapply(beta, getFitted)
  Xn <- model.matrix(~ g * x + 0)
  Q <- qr(Xn)
  H <- tcrossprod(solve(qr.R(Q)), qr.Q(Q))
  getCoef <- function(f) H %*% f
  Coef <- lapply(fitted, getCoef)
  group.slopes <- function(B){ # B of form ~ group * x + 0
    gp <- qr(model.matrix(~g))$rank
    gnames <- rownames(B)[1:gp]
    B <- as.matrix(B[-(1:gp),])
    B[2:gp, ] <- B[2:gp, ] + matrix(rep(B[1,], each = gp -1), gp -1, p)
    rownames(B) <- gnames
    B
  }
  slopes <- lapply(Coef, group.slopes)
  rename <- function(x) {
    dimnames(x)[[1]] <- levels(g)
    x
  }
  slopes <- lapply(slopes, rename)
  slopes
}

# getLSmeans
# gets the LS means for groups from a lm.rrpp fit, after constaining covariates to mean values
# used in pairwise
getLSmeans <- function(fit, g){
  k <- length(fit$LM$term.labels)
  n <- fit$LM$n
  beta <- fit$LM$random.coef[[k]]
  dat <- fit$LM$data
  covCheck <- sapply(dat, class)
  for(i in 1:length(covCheck)) if(covCheck[i] == "numeric") dat[[i]] <- mean(dat[[i]])
  L <- model.matrix(fit$LM$Terms, data = dat)
  L <- L[, colnames(L)  %in% rownames(beta[[1]])]
  getFitted <- function(b) L %*% b
  fitted <- lapply(beta, getFitted)
  Xn <- model.matrix(~ g + 0)
  Q <- qr(Xn)
  H <- tcrossprod(solve(qr.R(Q)), qr.Q(Q))
  getCoef <- function(f) H %*% f
  means <- lapply(fitted, getCoef)
  rename <- function(x) {
    dimnames(x)[[1]] <- levels(g)
    x
  }
  means <- lapply(means, rename)
  means
}


#' Support function for RRPP
#'
#' Calculate vector correlations for a matrix (by rows).  Used for pairwise comparisons.
#'
#' @param M Matrix for vector correlations.
#' @keywords utilities
#' @export
#' @author Michael Collyer
vec.cor.matrix <- function(M) {
  options(warn = -1)
  M = as.matrix(M)
  w = 1/sqrt(rowSums(M^2))
  vc = tcrossprod(M*w)
  diag(vc) <- 1
  options(warn = 0)
  vc
}

# Pval.list
# P-values across a list
# used in pairwise
Pval.list <- function(M){
  pvals <- M[[1]]
  n <- length(M)
  for(i in 1:length(pvals)) {
    y <- sapply(1:n, function(j) M[[j]][i])
    pvals[i] <- pval(y)
  }
  diag(pvals) <- 1
  pvals
}

# effect.list
# effect size across a list
# used in pairwise
effect.list <- function(M){
  Z <- M[[1]]
  n <- length(M)
  for(i in 1:length(Z)) {
    y <- sapply(1:n, function(j) M[[j]][i])
    Z[i] <- effect.size(y)
  }
  diag(Z) <- 0
  Z
}

# percentile.list
# find percentiles across a list
# used in pairwise
percentile.list <- function(M, confidence = 0.95){
  P <- M[[1]]
  n <- length(M)
  for(i in 1:length(P)) {
    y <- sapply(1:n, function(j) M[[j]][i])
    P[i] <- quantile(y, confidence, na.rm = TRUE)
  }
  P
}

# d.summary.from.list
# find distance statistics from a list
# used in pairwise
d.summary.from.list <- function(M, confidence = 0.95){
  P <- Pval.list(M)
  Z <- effect.list(M)
  CL <- percentile.list(M, confidence)
  list(D=M[[1]], P=P, Z=Z, CL=CL, confidence = confidence)
}

# d.summary.from.list
# find vec correlation statistics from a list
# used in pairwise
r.summary.from.list <- function(M, confidence = 0.95){
  options(warn = -1)
  acos.mat <- function(x){
    a <- acos(x)
    diag(a) <- 1
    a
  }
  A <- lapply(M, acos.mat)
  options(warn = 0)
  P <- Pval.list(A)
  Z <- effect.list(A)
  aCL <- percentile.list(A, confidence)
  angle = A[[1]]
  diag(angle) <- 0
  list(r=M[[1]], angle = angle,
       P=P, 
       Z=Z, aCL=aCL, confidence = confidence)
}

# makePWDTable
# arrange distance statistics into a table
# used in pairwise
makePWDTable <- function(L) { # List from summary.from.list
  nms <- rownames(L$D)
  DD <- as.dist(L$D)
  DP <- as.dist(L$P)
  DZ <- as.dist(L$Z)
  DC <- as.dist(L$CL)
  nam.com <- combn(length(nms), 2)
  name.list <- list()
  for(i in 1:NCOL(nam.com)) 
    name.list[[i]]  <- paste(nms[nam.com[1,i]], nms[nam.com[2,i]], sep =":")
  name.list <- unlist(name.list)
  tab <- data.frame(d = as.vector(DD),
                    UCI = as.vector(DC),
                    Z = as.vector(DZ), 
                    P = as.vector(DP))
  rownames(tab) <- name.list
  colnames(tab)[2] <- paste("UCL (", L$confidence*100,"%)", sep = "")
  colnames(tab)[4] <- "Pr > d"
  tab
}

# makePWCorTable
# arrange vec cor statistics into a table
# used in pairwise
makePWCorTable <- function(L){
  nms <- rownames(L$r)
  DR <- as.dist(L$r)
  DA <- as.dist(L$angle)
  DP <- as.dist(L$P)
  DaZ <- as.dist(L$Z)
  DaC <- as.dist(L$aCL)
  nam.com <- combn(length(nms), 2)
  name.list <- list()
  for(i in 1:NCOL(nam.com)) 
    name.list[[i]]  <- paste(nms[nam.com[1,i]], nms[nam.com[2,i]], sep =":")
  name.list <- unlist(name.list)
  tab <- data.frame(r = as.vector(DR),
                    angle = as.vector(DA),
                    UCL = as.vector(DaC),
                    Z = as.vector(DaZ),
                    P = as.vector(DP))
  rownames(tab) <- name.list
  colnames(tab)[3] <- paste("UCL (", L$confidence*100,"%)", sep = "")
  colnames(tab)[5] <- "Pr > angle"
  tab

}

# leaveOneOut
# set-up for jackknife classification
# used in classify
leaveOneOut <- function(X, Y, n.ind) {
  x <- X[-n.ind,]
  QR <- qr(x)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  H <- tcrossprod(solve(R), Q)
  y <- Y[-n.ind,]
  H %*% y
}

# RiReg
# Ridge Regularization of a covariance matrix, if needed
# used in classify
RiReg <- function(Cov, residuals){
  leads <- seq(0,1,0.005)[-1]
  leads <- leads[-length(leads)]
  I <- diag(1, NROW(Cov))
  N <- NROW(residuals)
  p <- NCOL(residuals)
  
  Covs <- lapply(1:length(leads), function(j){
    lambda <- leads[[j]]
    (lambda * Cov + (1 - lambda) * I)
  })
  logL <- sapply(1:length(leads), function(j){
    C <- Covs[[j]]
    N*p*log(2*pi) + N * determinant(C, logarithm = TRUE)$modulus[1] + sum(diag(
      residuals %*% fast.solve(C) %*% t(residuals) 
    ))
  })
  
  Covs[[which.min(logL)]]
  
}


logL <- function(fit){
  n <- fit$LM$n
  p <- fit$LM$p.prime
  X <- as.matrix(fit$LM$X * sqrt(fit$LM$weights))
  Y <- as.matrix(fit$LM$Y)
  rdf <- fit$LM$data
  if(fit$LM$gls){
    Sig <- (crossprod(fit$LM$gls.residuals, 
                      fast.solve(fit$LM$Cov)) %*%
              fit$LM$gls.residuals) /n
    s <- svd(Sig)
    pr <- which(cumsum(s$d)/sum(s$d) < 0.999)
    pp <- length(pr)
    if(p > 1) P <- as.matrix(Y %*% s$v[,pr]) else P <- center(Y)
    fit$LM$data$Y <- P
    pfit <- lm.rrpp(formula(fit$LM$Terms), print.progress = FALSE, 
                    Cov = fit$LM$Cov, data = fit$LM$data, 
                    weights = fit$LM$weights, iter = 0)
    Sig <- as.matrix((crossprod(pfit$LM$gls.residuals, fast.solve(pfit$LM$Cov)) %*%
              pfit$LM$gls.residuals)) / n
    if(kappa(Sig) > 1e10) Sig <- RiReg(Sig, pfit$LM$gls.residuals)
    
    ll <- -0.5*(n*pp + n*determinant(Sig, logarithm = TRUE)$modulus[1] + 
      pp*determinant(pfit$LM$Cov, logarithm = TRUE)$modulus[1] + n*pp*log(2*pi))
  }  else {
    
    Sig <- crossprod(fit$LM$wResiduals) /n
    s <- svd(Sig) 
    pr <- which(cumsum(s$d)/sum(s$d) < 0.999)
    pp <- length(pr)
    if(p > 1) P <- as.matrix(Y %*% s$v[,pr]) else P <- center(Y)
    fit$LM$data$Y <- P
    pfit <- lm.rrpp(formula(fit$LM$Terms), print.progress = FALSE, 
                    data = fit$LM$data, 
                    weights = fit$LM$weights, iter = 0)
    Sig <- as.matrix(crossprod(pfit$LM$residuals)) / n
    if(kappa(Sig) > 1e10) Sig <- RiReg(Sig, pfit$LM$residuals)
    ll <- -0.5*(n * pp + n * determinant(Sig, logarithm = TRUE)$modulus[1] + n * pp * log(2*pi))
  }
  
  ll
  
}

cov.trace <- function(fit) {
  n <- fit$LM$n
  p <- fit$LM$p.prime
  if(fit$LM$gls){
    Sig <- (crossprod(fit$LM$gls.residuals, 
                      fast.solve(fit$LM$Cov)) %*%
              fit$LM$gls.residuals) /n
    
  }  else {
    
    Sig <- crossprod(fit$LM$wResiduals) /n
    
  }
  
  sum(Sig^2)
  
}

z.test <- function(aov.mm){
  effect.type = aov.mm$effect.type
  if(effect.type == "F") stat <- aov.mm$F
  if(effect.type == "cohenf") stat <- aov.mm$F
  if(effect.type == "SS") stat <- aov.mm$SS
  if(effect.type == "MS") stat <- aov.mm$MS
  if(effect.type == "Rsq") stat <- aov.mm$Rsq
  
  perms <- ncol(stat)
  m <- nrow(stat)
  stat.c <- sapply(1:m, function(j){
    s <- stat[j,]
    center(s)
  })
  
  index <- combn(m, 2)
  Dz <- Pz <- dist(matrix(0, m))

  zdj <- function(x, y, j) {
    obs <- x[j] - y[j]
    sigd <- sqrt(var(x) + var(y))
    obs/sigd
  }
  
  for(i in 1:ncol(index)) {
    x <- stat.c[,index[1,i]]
    y <- stat.c[,index[2,i]]
    res <- array(NA, perms)
    for(j in 1: perms) res[j] <- zdj(x, y, j)
    Dz[i] <- abs(res[1])
    Pz[i] <- pval(abs(res))
  }

Z = as.matrix(Dz)
P = as.matrix(Pz)

options(warn = -1)
mds <- cmdscale(Z, m-1, eig = TRUE)
options(warn = 0)
    
list(Z = as.matrix(Dz), P = as.matrix(Pz), 
     mds = mds,
     form.names = rownames(aov.mm$table)[1:m],
     model.names = paste("m", 0:(m-1), sep = ""))
  
}

wilks <- function(e) {
  e <- zapsmall(e)
  e <- e[e > 0]
  prod(1/(1 + e))
}

pillai <- function(e) {
  e <- zapsmall(e)
  e <- e[e > 0]
  sum(e/(1 + e))
}

hot.law <- function(e) {
  e <- zapsmall(e)
  e <- e[e > 0]
  sum(e)
}

# trajsize
# find path distance of trajectory
# used in: trajectory.analysis
trajsize <- function(y) {
  k <- NROW(y[[1]])
  tpairs <- cbind(1:(k-1),2:k)
  sapply(1:length(y), function(j) {
    d <- as.matrix(dist(y[[j]]))
    sum(d[tpairs])
  })
}

# trajorient
# find trajectory correlations from first PCs
# used in: trajectory.analysis
trajorient <- function(y, tn) {
  m <- t(sapply(1:tn, function(j){
    x <- y[[j]]
    La.svd(center.scale(x)$coords, 0, 1)$vt
  }))
  vec.cor.matrix(m)
}

# trajshape
# find shape differences among trajectories
# used in: trajectory.analysis
trajshape <- function(y){
  y <- Map(function(x) center.scale(x)$coords, y)
  M <- Reduce("+",y)/length(y)
  z <- apply.pPsup(M, y)
  z <- t(sapply(z, function(x) as.vector(t(x))))
  as.matrix(dist(z))
}
