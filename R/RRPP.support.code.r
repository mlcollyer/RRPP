#' @name RRPP-package
#' @docType package
#' @aliases RRPP
#' @title Linear model evaluation with Randomized Residual Permutation Procedures
#' @author Michael Collyer and Dean Adams
#'
#' @description Functions in this package allow one to evaluate linear models with residual randomization.
#' The name, "RRPP", is an acronym for, "Randomized Residual Permutation Procedure."  Through
#' the various functions in this package, one can use randomization of residuals to generate empirical probability
#' distributions for linear model effects, for high-dimensional data or distance matrices.
#'
#' @import parallel
#' @importFrom stats as.dist 
#' @importFrom stats cmdscale 
#' @importFrom stats formula
#' @importFrom stats lm.fit
#' @importFrom stats lm.wfit 
#' @importFrom stats model.matrix
#' @importFrom stats prcomp
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @section RRPP TOC:
#' RRPP-package
NULL


#' Landmarks on pupfish
#'
#' @name pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprindon pecosensis body shapes, with indication of Sex and
#' Population from which fish were sampled (Marsh or Sinkhole).
#' @details These data were previously aligned
#' with GPA.  Centroid size (CS) is also provided.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic
#' change for phenotypes described by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
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

rrpp.data.frame<- function(...){
  dots <- list(...)
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
    if(is.data.frame(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[2]]
    if(is.vector(dots[[i]])) dots.ns[i] <- length(dots[[i]])
    if(is.factor(dots[[i]])) dots.ns[i] <- length(dots[[i]])
    if(is.logical(dots[[i]])) dots.ns[i] <- length(dots[[i]])
  }
  if(any(is.na(dots.ns))) stop("Some input is either dimensionless or inappropriate for data frames")
  if(length(unique(dots.ns)) > 1) stop("Inputs have different numbers of observations")
  class(dots) <- c("rrpp.data.frame")
  dots
}
#####----------------------------------------------------------------------------------------------------

# SUPPORT FUNCTIONS

# center
# centers a matrix faster than scale()
# used in various functions where mean-centering is required

center <- function(x){
  if(is.vector(x)) x- mean(x) else {
    x <- as.matrix(x)
    dims <- dim(x)
    fast.center(x, dims[1], dims[2])
  }
}

fast.center <- function(x, n, p){
  m <- colMeans(x)
  x - rep.int(m, rep_len(n, p))
}

fast.scale <- function(x, n, p){
  if(p > 1) {
    x <- fast.center(x, n, p)
    scale <- apply(x, 2, sd)
    x / rep.int(scale, rep_len(n, p))
  } else {
    x <- x - mean(x)
    x/sd(x)
  }
}

# fast.ginv
# same as ginv, but without traps (faster)
# used in any function requiring a generalized inverse
fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
  k <- ncol(X)
  Xsvd <- La.svd(X, k, k)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
  v <-t(Xsvd$vt)[, Positive, drop = FALSE]
  v%*%rtu
}

# fast.solve
# chooses between fast.ginv or qr.solve, when det might or might not be 0
# used in any function requiring a matrix inverse where the certainty of
# singular matrices is in doubt
fast.solve <- function(x) if(det(x) > 1e-8) chol2inv(chol(x)) else fast.ginv(x)

# pcoa
# acquires principal coordinates from distance matrices
# used in all linear model functions with data input
pcoa <- function(D){
  options(warn=-1)
  if(class(D) != "dist") stop("function only works with distance matrices")
  cmd <- cmdscale(D, k=attr(D, "Size") -1, eig=TRUE)
  options(warn=0)
  p <- length(cmd$eig[zapsmall(cmd$eig) > 0])
  Yp <- cmd$points[,1:p]
  Yp
}

# perm.index
# creates a permutation index for resampling
# used in all functions with a resampling procedure

perm.index <-function(n, iter, seed=NULL){
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      ind <- c(list(1:n),(Map(function(x) sample.int(n,n), 1:iter)))
      rm(.Random.seed, envir=globalenv())
      attr(ind, "seed") <- seed
      ind
}

# boot.index
# creates a bootstrap index for resampling
# used in lm.rrpp for intercept models
boot.index <-function(n, iter, seed=NULL){
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      ind <- c(list(1:n),(Map(function(x) sample.int(n, n, replace = TRUE), 1:iter)))
      rm(.Random.seed, envir=globalenv())
      attr(ind, "seed") <- seed
      ind
}

# fastFit
# calculates fitted values for a linear model, after decomoposition of X to get U
# used in SS.iter
fastFit <- function(U,y,n,p){
  if(!is.matrix(y)) y <- as.matrix(y)
  if(p > n) tcrossprod(U)%*%y else
    U%*%crossprod(U,y)
}

# fastLM
# calculates fitted values and residuals, after fastFit
# placeholder in case needed later
fastLM<- function(U,y){
  p <- dim(y)[2]; n <- dim(y)[1]
  yh <- fastFit(U,y,n,p)
  list(fitted = yh, residuals = y-yh)
}

# pval
# P-values form random outcomes
# any analytical function
pval = function(s){# s = sampling distribution
  p = length(s)
  r = rank(s)[1]-1
  pv = 1-r/p
  pv
}

# effect.size
# Effect sizes (standard deviates) form random outcomes
# any analytical function
effect.size <- function(x, center = TRUE) {
  z = scale(x, center=center)
  n <- length(z)
  z[1]*sqrt((n-1)/(n))
}

# rrpp.fit + subfunctions
# lm-like fit modified for all submodels
# general workhorse for all 'rrpp.lm' functions
# used in all 'rrpp.lm' functions

# rrpp.fit.lm
# base for rrpp.fit
rrpp.fit.lm <- function(a){
  # set-up
  w <- a$w
  if(any(w <= 0)) stop("Weights must be positive")
  o <-a$offset
  dat <- a$data
  Y <- dat$Y
  X <- a$x
  n <- NROW(Y)
  SS.type <- a$SS.type
  dat <- a$data
  wY <- Y*sqrt(w); wX <- X*sqrt(w)
  # data and design matrix
  Terms <- a$terms
  X.k <- attr(X, "assign")
  QRx <- qr(X)
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0,X.k))
  k <- length(attr(Terms, "term.labels"))
  # SS types: reduced and full X matrices
  if(SS.type == "III"){
    Xrs <- lapply(2:length(uk), function(j)  X[, X.k %in% uk[-j]])
    Xfs <- lapply(2:length(uk), function(j)  X)
  }
  if(SS.type == "II") {
    fac <- crossprod(attr(Terms, "factor"))
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
    Xs <- lapply(1:length(uk), function(j)  Xj <- X[, X.k %in% uk[1:j]])
    Xrs <- Xs[1:k]
    Xfs <- Xs[2:(k+1)]
  }
  # unweighted output
  QRs.reduced <- lapply(Xrs, function(x) qr(x))
  fits.reduced <- lapply(Xrs, function(x) lm.fit(as.matrix(x),Y, offset = o))
  fitted.reduced <- lapply(fits.reduced, function(x) as.matrix(x$fitted.values))
  residuals.reduced <- lapply(fits.reduced, function(x) as.matrix(x$residuals))
  coefficients.reduced <- lapply(fits.reduced, function(x) as.matrix(x$coefficients))
  QRs.full <- lapply(Xfs, function(x) qr(x))
  fits.full <- lapply(Xfs, function(x) lm.fit(as.matrix(x),Y, offset = o))
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
  term.labels <- attr(Terms, "term.labels")
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
              SS.type = SS.type,
              Terms = Terms, term.labels = term.labels)
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
  fits.full <- lm.fit(as.matrix(X),Y, offset = o)
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
              SS.type = NULL,
              Terms = Terms, term.labels = term.labels)
  class(out) <- "procD.fit"
  invisible(out)
}

# rrpp.fit
# calls one of previous functions, depending on conditions
rrpp.fit <- function(f1, keep.order=FALSE, pca=TRUE,
                     SS.type = NULL, data = NULL, ...){
  dots <- list(...)
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","II", "III")))) SS.type <- "I"
  if(any(inherits(f1, "lm"))) {
    d <- f1$model
    form <- formula(terms(f1), keep.order = keep.order)
    form.adj <- update(form, Y ~.)
    form[[2]] <- form.adj[[2]]
    Terms <- terms(form, keep.order = keep.order)
    tl <- attr(Terms, "term.labels")
    if(length(tl) == 0){
      dat <- data.frame(Y = 1:NROW(d))
      dat$Y <- as.matrix(d)
    } else {
      dat <- lapply(1:length(tl),
                    function(j) try(get(as.character(tl[j]),
                                        d), silent = TRUE))
      names(dat) <- tl
      dat <- as.data.frame(dat)
    }
    x <- model.matrix(Terms, data = dat)
    w <- f1$weights
    o <- f1$offset
    t <- f1$terms
    f <- form
    pdf.args <- list(data=dat, x=x, w=w, offset=o, terms=t, formula=f,
                     SS.type = SS.type)
  } else {
    form.in <- formula(f1)
    d <- list()
    d$Y <- eval(form.in[[2]], data, parent.frame())
    if(inherits(d$Y, "dist")) d$Y <- pcoa(d$Y) else
      if ((is.matrix(d$Y) || is.data.frame(d$Y))
          && isSymmetric(unname(as.matrix(d$Y)))) d$Y <- pcoa(as.dist(d$Y)) else
            d$Y <- as.matrix(d$Y)
    n <- NROW(d$Y)
    form <- formula(terms(f1), keep.order = keep.order)
    form.adj <- update(form, Y ~.)
    form[[2]] <- form.adj[[2]]
    Terms <- terms(form, keep.order=keep.order)
    tl <- unique(unlist(strsplit(attr(Terms, "term.labels"), ":")))
    log.check <- grep("log", tl)
    scale.check <- grep("scale", tl)
    exp.check <- grep("exp", tl)
    poly.check <- grep("poly", tl)
    if(length(log.check) > 0) for(i in 1:length(log.check)){
      tlf <- tl[log.check[i]]
      tlf <- gsub("log\\(", "", tlf)
      tlf <- gsub("\\)", "", tlf)
      tl[log.check[i]]<- tlf
    }
    if(length(scale.check) > 0) for(i in 1:length(scale.check)){
      tlf <- tl[scale.check[i]]
      tlf <- gsub("scale\\(", "", tlf)
      tlf <- gsub("\\)", "", tlf)
      tl[scale.check[i]]<- tlf
    }
    if(length(exp.check) > 0) for(i in 1:length(exp.check)){
      tlf <- tl[exp.check[i]]
      tlf <- gsub("exp\\(", "", tlf)
      tlf <- gsub("\\)", "", tlf)
      tl[exp.check[i]]<- tlf
    }
    if(length(poly.check) > 0) for(i in 1:length(poly.check)){
      tlf <- tl[poly.check[i]]
      tlf <- gsub("poly\\(", "", tlf)
      tlf <- gsub("\\)", "", tlf)
      tlf <- strsplit(tlf, "")[[1]][1]
      tl[poly.check[i]]<- tlf
    }
    if(length(tl) == 0){
      dat <- data.frame(Y = 1:n)
      dat$Y <- d$Y
    } else {
      if(is.null(data))
        dat <- lapply(1:length(tl), function(j) {
          try(get(as.character(tl[j]), parent.frame()),
              silent = TRUE)
        }) else
          dat <- lapply(1:length(tl), function(j) {
            try(get(as.character(tl[j]), data),
                silent = TRUE)
          })
        check <- (sapply(dat, NROW) == n)
        dat <- dat[check]
        tl <- tl[check]
        names(dat) <- tl
        dat <- as.data.frame(dat)
    }
    if(pca) d$Y <- prcomp(d$Y)$x
    dat$Y <- d$Y
    pdf.args <- list(data=dat,
                     x = model.matrix(Terms, data = dat),
                     w = dots$weights,
                     offset = dots$offset,
                     terms = Terms,
                     formula = form,
                     SS.type = SS.type)
  }
  if(is.null(pdf.args$w))
    pdf.args$w <- rep(1, NROW(pdf.args$data))
  if(is.null(pdf.args$offset))
    pdf.args$offset <- rep(0, NROW(pdf.args$data))
  if(sum(attr(pdf.args$x, "assign")) == 0)
    out <- rrpp.fit.int(pdf.args) else
      out <- rrpp.fit.lm(pdf.args)
  out
}

rrpp.fit.from.lm.rrpp <- function(fit){
  # get the rrpp.args needed
  # then decide to use rrpp.fit.lm or rrpp.fit.int
}



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
  if(!is.null(w) && !is.null(o)) rrpp.w.o(fitted, residuals, ind.i, w, o)
  if(!is.null(w) && is.null(o)) rrpp.w(fitted, residuals, ind.i, w)
  if(is.null(w) && !is.null(o)) rrpp.o(fitted, residuals, ind.i, o)
  if(is.null(w) && is.null(o)) rrpp.basic(fitted, residuals, ind.i)
}

# Cov.proj
# generates projection matrix from covariance matrix
# used in lm.rrpp

Cov.proj <- function(Cov, id){
  if(is.null(id)) id <- 1:NCOL(Cov)
  Cov <- Cov[id, id]
  invC <- fast.solve(Cov)
  eigC <- eigen(Cov)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[,1:(length(lambda))]
  P <- fast.solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect))
  dimnames(P) <- dimnames(Cov)
  P
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
  fitted <- fit$fitted.reduced
  res <- fit$residuals.reduced
  Y <- fit$Y
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  trms <- fit$term.labels
  k <- length(trms)
  fF <- fit$fitted.full[[k]]
  rF <- fit$residuals.full[[k]]
  if(!RRPP) {
    fitted <- lapply(fitted, function(.) matrix(0, n, p))
    res <- lapply(res, function(.) Y)
  }
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(w) != n) weighted = TRUE else weighted = FALSE
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  rrpp.args <- list(fitted = fitted, residuals = res,
                    ind.i = NULL, w = NULL, o = NULL)
  if(weighted) rrpp.args$w <- w
  if(offset) rrpp.args$o <- o
  if(print.progress){
    cat(paste("\n\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  if(gls){
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Ur <- lapply(Xr, function(x) crossprod(P,qr.Q(qr(x))))
    Uf <- lapply(Xf, function(x) crossprod(P,qr.Q(qr(x))))
    Ufull <- Uf[[k]]
    SS <- lapply(1: perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      if(weighted) y <- (fF + rF[x,])*w else y <- fF + rF[x,]
      py <- crossprod(P,y); pyy <- sum(py^2)
      c(Map(function(y, ur, uf) sum(crossprod(uf,y)^2) - sum(crossprod(ur,y)^2),
            Yi, Ur, Uf),
        pyy - sum(crossprod(Ufull, y)^2), pyy - SS.mean(py, n))
    })
  } else {
    Ur <- lapply(fit$wQRs.reduced, qr.Q)
    Uf <- lapply(fit$wQRs.full, qr.Q)
    Ufull <- Uf[[k]]
    SS <- lapply(1: perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      if(weighted) y <- (fF + rF[x,])*w else y <- fF + rF[x,]
      yy <- sum(y^2)
      c(Map(function(y, ur, uf) sum(crossprod(uf,y)^2) - sum(crossprod(ur,y)^2),
            Yi, Ur, Uf),
        yy - sum(crossprod(Ufull, y)^2), yy - SS.mean(y, n))
    })
  }
  SS <- matrix(unlist(SS), k+2, perms)
  rownames(SS) <- c(trms, "Residuals", "Total")
  colnames(SS) <- c("obs", paste("iter", 1:(perms-1), sep=":"))
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  SS
}

# SS.iterPP
# workhorse for lm.rrpp, same as SS.iter, but with parallel processing
# used in lm.rrpp

SS.iterPP <- function(fit, ind, P = NULL, RRPP = TRUE, print.progress = TRUE) {
  cl <- detectCores()-1
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  if(print.progress){
    cat(paste("\n\nSums of Squares calculations:", perms, "permutations.\n"))
    cat(paste("Progress bar not possible with parallel processing, but this shouldn't take long...\n"))
  }
  fitted <- fit$fitted.reduced
  res <- fit$residuals.reduced
  Y <- fit$Y
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  trms <- fit$term.labels
  k <- length(trms)
  fF <- fit$fitted.full[[k]]
  rF <- fit$residuals.full[[k]]
  if(!RRPP) {
    fitted <- lapply(fitted, function(.) matrix(0, n, p))
    res <- lapply(res, function(.) Y)
  }
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(sum(w) != n) weighted = TRUE else weighted = FALSE
  if(sum(o) != 0) offset = TRUE else offset = FALSE
  rrpp.args <- list(fitted = fitted, residuals = res,
                    ind.i = NULL, w = NULL, o = NULL)
  if(weighted) rrpp.args$w <- w
  if(offset) rrpp.args$o <- o
  step <- 2
  if(gls){
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Ur <- lapply(Xr, function(x) crossprod(P,qr.Q(qr(x))))
    Uf <- lapply(Xf, function(x) crossprod(P,qr.Q(qr(x))))
    Ufull <- Uf[[k]]
    Unull <- matrix(rowSums(crossprod(P,qr.Q(qr(matrix(rowSums(P))))))*w)
    SS <- mclapply(1: perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      if(weighted) y <- (fF + rF[x,])*w else y <- fF + rF[x,]
      py <- crossprod(P,y); pyy <- sum(py^2)
      c(Map(function(y, ur, uf) sum(crossprod(uf,y)^2) - sum(crossprod(ur,y)^2),
            Yi, Ur, Uf),
        pyy - sum(crossprod(Ufull, y)^2), pyy - SS.mean(py, n))
    }, mc.cores = cl)
  } else {
    Ur <- lapply(fit$wQRs.reduced, qr.Q)
    Uf <- lapply(fit$wQRs.full, qr.Q)
    SS <- mclapply(1: perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      if(weighted) y <- (fF + rF[x,])*w else y <- fF + rF[x,]; yy <- sum(y^2)
      c(Map(function(y, ur, uf) sum(crossprod(uf,y)^2) - sum(crossprod(ur,y)^2),
            Yi, Ur, Uf),
        yy - sum(crossprod(Uf[[k]], y)^2), yy - SS.mean(y, n))
    }, mc.cores = cl)
  }
  SS <- matrix(unlist(SS), k+2, perms)
  rownames(SS) <- c(trms, "Residuals", "Total")
  colnames(SS) <- c("obs", paste("iter", 1:(perms-1), sep=":"))
  SS
}

SS.iter.null <- function(fit, ind, P = NULL, print.progress = TRUE) {
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  if(print.progress){
    cat(paste("\n\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  fitted <- fit$fitted.reduced
  res <- fit$residuals.reduced
  Y <- fit$Y
  dims <- dim(as.matrix(Y))
  n <- dims[1]; p <- dims[2]
  k <- 1
  w <- sqrt(fit$weights)
  o <- fit$offset
  if(gls){
    SS <- lapply(1:perms, function(j){
      x <-ind[[j]]
      y <- Y[x,]*w; py <- crossprod(P,y); pyy <- sum(py^2)
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      pyy - SS.mean(py, n)
    })
  } else {
    SS <- lapply(1:perms, function(j){
      x <-ind[[j]]
      y <- Y[x,]*w; yy <- sum(y^2)
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      yy - SS.mean(y, n)
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
  dft <- sum(df, dfe)
  df <- c(df, dfe, dft)
  MS <- SS/df
  Fs <- (apply(MS[1:(k+1),], 2, function(x) x/x[length(x)]))[1:k,]
  Rsq <- apply(SS, 2, function(x) x/x[length(x)])
  if(SS.type == "III") {
    etas <- apply(SS, 2, function(x) x/(x + x[k+1]))
    cohenf <- (etas/(1-etas))[1:k,]
  } else {
    etas <- Rsq[1:k,]
    if(k == 1) unexp <- 1 - etas else unexp <- 1 - apply(etas, 2, cumsum)
    cohenf <- etas/unexp
  }
  if(k == 1) {
    P.val <- pval(Fs)
    SS.effect <- effect.size(log(SS[1,]))
    MS.effect <- effect.size(log(MS[1,]))
    Rsq.effect <- effect.size(log(Rsq[1,]))
    F.effect <- effect.size(log(Fs))
    cohenf.effect <- effect.size(log(cohenf))
  } else{
    P.val <- apply(Fs, 1, pval)
    SS.effect <- apply(log(SS[1:k,]), 1, effect.size)
    MS.effect <- apply(log(MS[1:k,]), 1, effect.size)
    Rsq.effect <- apply(log(Rsq[1:k,]), 1, effect.size)
    cohenf.effect <- apply(log(cohenf), 1, effect.size)
    F.effect <- apply(log(Fs), 1, effect.size)
  }
  Fs <- rbind(Fs, NA, NA)
  cohenf <- rbind(cohenf, NA, NA)
  rownames(Fs) <- rownames(cohenf) <- rownames(SS)
  out <- list(SS.type = SS.type, SS = SS, MS = MS, Rsq = Rsq,
              Fs = Fs, cohenf = cohenf, P.val = P.val, SS.effect = SS.effect,
              MS.effect = MS.effect, Rsq.effect = Rsq.effect, F.effect = F.effect,
              cohenf.effect = cohenf.effect,
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

# fastPC
# get PC scores in p dimensions faster than going through prcomp.  Good for HD data
# used in compare.models

fastPC <- function(y, p){
  y <- center(y)
  s <-La.svd(y, p, p)
  tcrossprod(y, s$vt)
}

# detPC
# finds determinant of singular covariance matrix.  Finds the relevant eigenvalues for computation.
# 99.9% is used as a cutoff to prevent singular dimensions
# used in compare.models

detPC <- function(Sig){
  s <- svd(Sig)
  d <- s$d
  dp <- d[which(cumsum(d/sum(d)) <= 0.999)]
  prod(dp)
}



# DevianceOLS
# calculates Deviances from RSS output in lm.rrpp, based on OLS
# used in compare.models

# refit
# finds rrpp.fit from lm.rrpp object
# used in Deviance
refit <- function(fit){
  dat <- fit$LM$data
  dims <- dim(dat$Y)
  dat$Y <- prcomp(dat$Y)$x[,1:min(dims)]
  pdf.args <- list(data=dat,
                   x = model.matrix(fit$LM$Terms, data = dat),
                   w = fit$LM$weights,
                   offset = fit$LM$offset,
                   terms = fit$LM$Terms,
                   formula = NULL,
                   SS.type = fit$ANOVA$SS.type)
  form <- formula(deparse(fit$call[[2]]))
  form[[2]] <- "Y"
  pdf.args$formula <- form
  if(sum(attr(pdf.args$x, "assign")) == 0)
    out <- rrpp.fit.int(pdf.args) else
      out <- rrpp.fit.lm(pdf.args)
  out
}

# remodel.set
# set up for remodel.SS, from bootstrap iterations
# used in Deviance
remodel.set <- function(fit){
    k <- length(fit$term.labels)
    if(k == 0) k = 1
    f <- fit$fitted.full[[k]]
    r <- fit$residuals.full[[k]]
    w <- sqrt(fit$weights)
    o <- fit$offset
    U <- qr.Q(fit$wQRs.full[[k]])
  list(f = f, r =r, w = w, o = o, 
       n = NROW(r), p = NCOL(r), U = U, ind.i = NULL)
}


remodel.res<- function(f, r, w, o, n, p, U, ind.i){
  y <- f + r[ind.i,] - o
  y - fastFit(U, (y * w), n, p)
}

remodel <- function(f, r, w, o, n, p, U, ind.i){
  y <- f + r[ind.i,] - o
  f <- fastFit(U, (y * w), n, p)
  r <- y - f
  list(y = y + o, f = f + o, r = r)
}


remodel.Y<- function(f, r, w, o, n, p, U, ind.i){
  (f + r[ind.i,])*w - o
}

remodel.logL <- function(res){ # n and p part of Deviance function
  n <- NROW(res)
  p <- qr(res)$rank
  Sig <- crossprod(res)/n
  -0.5*n*(log(detPC(Sig)) + p*log(2*pi))
}

remodel.logL.GLS <- function(res, C){ # n and p part of Deviance function
  n <- NROW(res)
  p <- qr(res)$rank
  Sig <- crossprod(res)/n
  -0.5*n*(log(detPC(Sig)) + p*log(2*pi)) -0.5*log(detPC(C))
}


DevianceOLS<- function(fit){
  ind <- fit$PermInfo$perm.schedule
  perms <- length(ind)
  n <- length(ind[[1]])
  fit.args <- remodel.set(refit(fit))
  U <- fit.args$U
  fit.args$r <- res <- as.matrix(fit.args$r)
  wRes <- res*fit.args$w
  k <- fit$LM$QR$rank
  d <- sapply(1:perms, function(j){
    fit.args$ind.i <- ind[[j]]
    nM <- do.call(remodel, fit.args)
    dev <- sum(nM$r^2)
    ss <- sum(nM$f^2)
    c(deviance = dev, SS = ss)
  })
  SS <- d[2,]
  Dev <- d[1,]
  eta.sq <- SS/(SS + Dev)
  eta.sq.a <- eta.sq - (1-eta.sq)*(k-1)/(n-k)
  p <- qr(center(fit$LM$Y))$rank
  pA <- p*k + p*(p+1)/2
  logL.obs <- remodel.logL(wRes)
  AIC <-  -2*logL.obs + 2*pA
  list(Deviance = Dev, Eta.sq = eta.sq, Eta.sq.a = eta.sq.a, logL=logL.obs, AIC=AIC,
       p = p, k = k, n = n)
}

# DevianceGLS
# calculates Deviances from RSS output in lm.rrpp, based on GLS
# used in compare.models

DevianceGLS <- function(fit){
  ind <- fit$PermInfo$perm.schedule
  perms <- length(ind)
  n <- length(ind[[1]])
  fit.args <- remodel.set(refit(fit))
  P <- fit$LM$Pcov
  C <- fit$LM$Cov
  X <- fit$LM$X
  X <- crossprod(P, X*fit.args$w)
  U <- crossprod(P, qr.Q(qr(X)))
  fit.args$U <- U
  fit.args$r <- as.matrix(fit.args$r)
  Y <- (fit.args$f + fit.args$r)*fit.args$w
  wRes <- Y - fastFit(U, Y, n, fit.args$p)
  k <- fit$LM$QR$rank
  d <- sapply(1:perms, function(j){
    fit.args$ind.i <- ind[[j]]
    nM <- do.call(remodel, fit.args)
    dev <- sum(nM$r^2)
    ss <- sum(nM$f^2)
    c(deviance = dev, SS = ss)
  })
  SS <- d[2,]
  Dev <- d[1,]
  eta.sq <- SS/(SS + Dev)
  eta.sq.a <- eta.sq - (1-eta.sq)*(k-1)/(n-k)
  p <- qr(center(Y))$rank
  r <- fit$LM$Cov.par
  pA <- p*k + p*(p+1)/2 + r
  logL.obs <- remodel.logL.GLS(wRes, C)
  AIC <-  -2*logL.obs + 2*pA
  list(Deviance = Dev, Eta.sq = eta.sq, Eta.sq.a = eta.sq.a, logL=logL.obs, AIC=AIC,
       p = p, k = k, n = n)
}
