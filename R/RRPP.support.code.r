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
#' @importFrom stats lm
#' @importFrom stats loess
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats qnorm
#' @importFrom graphics arrows
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar

#' @section RRPP TOC:
#' RRPP-package
NULL


#' Landmarks on pupfish
#'
#' @name Pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprindon pecosensis body shapes, with indication of Sex and
#' Population from which fish were sampled (Marsh or Sinkhole).
#' @details These data were previously aligned with GPA.  Centroid size (CS) is also provided.  
#' See the \link[geomorph]{geomorph} package for details.
#' 
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic
#' change for phenotypes described by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
NULL

#' Landmarks on pupfish heads
#'
#' @name PupfishHeads
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprindon pecosensis head shapes, with variables for 
#' sex, month and year sampled, locality, head size, and coordinates of landmarks for head shape,
#' per specimen.  These data are a subset of a larger data set.
#' @details The variable, "coords", are data that were previously aligned
#' with GPA.  The variable, "headSize", is the Centroid size of each vector of coordinates.
#' See the \link[geomorph]{geomorph} package for details.
#' @references Gilbert, M.C. 2016. Impacts of habitat fragmentation on the cranial morphology of a 
#' threatened desert fish (Cyprinodon pecosensis). Masters Thesis, Western Kentucky University.

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
      if(is.data.frame(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[2]]
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
  d <- cmd$eig
  min.d <- min(d)
  if(min.d < 0) {
    dcor <- abs(min.d)
    Dc <- D^2 + 2*dcor
    Dc <- sqrt(Dc)
    cmd.c <- cmdscale(Dc, k=attr(Dc, "Size") -1, eig=TRUE)
    options(warn=0)
    d <- cmd.c$eig
  } else cmd.c <- cmd
  p <- length(cmd.c$eig[zapsmall(d) > 0])
  Yp <- cmd.c$points[,1:p]
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
  if(!is.null(w) && !is.null(o)) r <- rrpp.w.o(fitted, residuals, ind.i, w, o)
  if(!is.null(w) && is.null(o)) r <- rrpp.w(fitted, residuals, ind.i, w)
  if(is.null(w) && !is.null(o)) r <- rrpp.o(fitted, residuals, ind.i, o)
  if(is.null(w) && is.null(o)) r <- rrpp.basic(fitted, residuals, ind.i)
  r
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
    cat("\nWarning: singular covariance matrix. Proceed with caution\n")
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
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  if(gls){
    Y <- crossprod(P, Y)
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Ur <- lapply(Xr, function(x) qr.Q(qr(x)))
    Uf <- lapply(Xf, function(x) qr.Q(qr(x)))
    Ufull <- Uf[[k]]
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      fitted <- Map(function(u) crossprod(tcrossprod(u), Y), Ur)
      res <- lapply(fitted, function(f) Y - f)
    }
    rrpp.args$fitted <- fitted
    rrpp.args$residuals <- res
    SS <- lapply(1: perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- Yi[[1]]
      pyy <- sum(y^2)
      c(Map(function(y, ur, uf) sum(crossprod(uf,y)^2) - sum(crossprod(ur,y)^2),
            Yi, Ur, Uf),
        pyy - sum(crossprod(Ufull, y)^2), pyy - SS.mean(y, n))
    })
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
      rrpp.args$fitted <- fitted
      rrpp.args$residuals <- res
    }
    Ur <- lapply(fit$wQRs.reduced, qr.Q)
    Uf <- lapply(fit$wQRs.full, qr.Q)
    Ufull <- Uf[[k]]
    SS <- lapply(1: perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- Yi[[1]]
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
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    cat(paste("Progress bar not possible with parallel processing, but this shouldn't take long...\n"))
  }
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
  step <- 2
  if(gls){
    Y <- crossprod(P, Y)
    Xr <- lapply(fit$wXrs, function(x) crossprod(P, as.matrix(x)))
    Xf <- lapply(fit$wXfs, function(x) crossprod(P, as.matrix(x)))
    Ur <- lapply(Xr, function(x) qr.Q(qr(x)))
    Uf <- lapply(Xf, function(x) qr.Q(qr(x)))
    Ufull <- Uf[[k]]
    Unull <- Ufull[,1]
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    } else {
      fitted <- Map(function(u) crossprod(tcrossprod(u), Y), Ur)
      res <- lapply(fitted, function(f) Y - f)
    }
    rrpp.args$fitted <- fitted
    rrpp.args$residuals <- res
    SS <- mclapply(1: perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- Yi[[1]]
      pyy <- sum(y^2)
      c(Map(function(y, ur, uf) sum(crossprod(uf,y)^2) - sum(crossprod(ur,y)^2),
            Yi, Ur, Uf),
        pyy - sum(crossprod(Ufull, y)^2), pyy - SS.mean(y, n))
    }, mc.cores = cl)
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
      rrpp.args$fitted <- fitted
      rrpp.args$residuals <- res
    }
    Ur <- lapply(fit$wQRs.reduced, qr.Q)
    Uf <- lapply(fit$wQRs.full, qr.Q)
    SS <- mclapply(1: perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- Yi[[1]]; yy <- sum(y^2)
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

SS.iter.null <- function(fit, ind, P = NULL, RRPP=TRUE, print.progress = TRUE) {
  if(!is.null(P)) gls = TRUE else gls = FALSE
  perms <- length(ind)
  if(print.progress){
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  fitted <- fit$wFitted.full
  res <- fit$wResiduals.full
  Y <- fit$wY
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
      U <- qr.Q(qr(crossprod(P, matrix(1, n))))
      fitted <- crossprod(tcrossprod(U), Y)
      res <- Y - fitted
    }
    SS <- lapply(1:perms, function(j){
      x <-ind[[j]]
      y <- fitted + res[x,]; pyy <- sum(y^2)
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      pyy - SS.mean(y, n)
    })
  } else {
    if(!RRPP) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) Y)
    }
    
    SS <- lapply(1:perms, function(j){
      x <-ind[[j]]
      y <- Y[x,]; yy <- sum(y^2)
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
    SS.effect <- effect.size(log(SS[1,]))
    MS.effect <- effect.size(log(MS[1,]))
    Rsq.effect <- effect.size(log(Rsq[1,]))
    F.effect <- effect.size(log(Fs))
    cohenf.effect <- effect.size(log(cohenf))
  } else{
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
              Fs = Fs, cohenf = cohenf, SS.effect = SS.effect,
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

# beta.lm.iter
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
    Ur <- lapply(Qr, function(x) qr.Q(x))
    Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
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
  beta.names <- lapply(1:length(beta.mats), function(j) rownames(beta.mats[[j]][[1]]))
  names(beta.mats) <- names(beta.mat.d) <- trms
  beta.hyp <- list()
  beta.hyp[[1]] <- beta.names[[1]]
  if(k > 1) for(i in 2:k) 
    beta.hyp[[i]] <- beta.names[[i]][which(!beta.names[[i]]%in%beta.names[[i-1]])]
  d.stitched <- lapply(1:k, function(j){
    b <- beta.mat.d[[j]]
    b <- b[match(beta.hyp[[j]], row.names(b)),]
    if(is.matrix(b)) t(b) else b
  })
  d.stitched <- t(matrix(unlist(d.stitched), perms, length(beta.names[[k]])))
  rownames(d.stitched) <- beta.names[[k]]
  colnames(d.stitched) <- colnames(beta.mat.d[[1]])
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
    Ur <- lapply(Qr, function(x) qr.Q(x))
    Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
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
  beta.names <- lapply(1:length(beta.mats), function(j) rownames(beta.mats[[j]][[1]]))
  names(beta.mats) <- names(beta.mat.d) <- trms
  beta.hyp <- list()
  beta.hyp[[1]] <- beta.names[[1]]
  if(k > 1) for(i in 2:k) 
    beta.hyp[[i]] <- beta.names[[i]][which(!beta.names[[i]]%in%beta.names[[i-1]])]
  d.stitched <- lapply(1:k, function(j){
    b <- beta.mat.d[[j]]
    b <- b[match(beta.hyp[[j]], row.names(b)),]
    if(is.matrix(b)) t(b) else b
  })
  d.stitched <- t(matrix(unlist(d.stitched), perms, length(beta.names[[k]])))
  rownames(d.stitched) <- beta.names[[k]]
  colnames(d.stitched) <- colnames(beta.mat.d[[1]])
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
