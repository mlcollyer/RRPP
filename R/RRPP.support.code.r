#' @name RRPP-package
#' @docType package
#' @aliases RRPP
#' @title Linear Model Evaluation with Randomized Residual Permutation Procedures
#' @author Michael Collyer and Dean Adams
#' @return Key functions for this package:
#' \item{\code{\link{lm.rrpp}}}{Fits linear models, using RRPP.}
#' \item{\code{\link{anova.lm.rrpp}}}{ANOVA on linear models, using RRPP, 
#' plus model comparisons.}
#' \item{\code{\link{coef.lm.rrpp}}}{Extract coefficients or perform test 
#' on coefficients, using RRPP.}
#' \item{\code{\link{predict.lm.rrpp}}}{Predict values from lm.rrpp fits 
#' and generate bootstrapped confidence intervals.}
#' \item{\code{\link{pairwise}}}{Perform pairwise tests, based on lm.rrpp 
#' model fits.}
#' 
#' @description Functions in this package allow one to evaluate linear models 
#' with residual randomization.
#' The name, "RRPP", is an acronym for, "Randomization of Residuals in a 
#' Permutation Procedure."  Through
#' the various functions in this package, one can use randomization of 
#' residuals to generate empirical probability
#' distributions for linear model effects, for high-dimensional data or 
#' distance matrices.
#' 
#' An especially useful option of this package is to fit models with 
#' either ordinary or generalized
#' least squares estimation (OLS or GLS, respectively), using theoretic 
#' covariance matrices.  Mixed linear
#' effects can also be evaluated.
#' 
#' @import parallel
#' @import stats
#' @import graphics
#' @import utils
#' @import ggplot2 
#' @import Matrix
#' @importFrom ape multi2di.phylo
#' @importFrom ape root.phylo
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
NULL
#' @section RRPP TOC:
#' RRPP-package
NULL

#' Landmarks on pupfish
#'
#' @name Pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprinodon pecosensis body shapes, with 
#' indication of Sex and
#' Population from which fish were sampled (Marsh or Sinkhole).
#' @details These data were previously aligned with GPA.  Centroid size (CS) 
#' is also provided.  
#' See the \pkg{geomorph} package for details.
#' 
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for 
#' analysis of phenotypic
#' change for phenotypes described by high-dimensional data. Heredity. 113: 
#' doi:10.1038/hdy.2014.75.
NULL

#' Landmarks on pupfish heads
#'
#' @name PupfishHeads
#' @docType data
#' @author Michael Collyer
#' @description Landmark data from Cyprinodon pecosensis head shapes, with 
#' variables for 
#' sex, month and year sampled, locality, head size, and coordinates of 
#' landmarks for head shape,
#' per specimen.  These data are a subset of a larger data set.
#' @details The variable, "coords", are data that were previously aligned
#' with GPA.  The variable, "headSize", is the Centroid size of each vector 
#' of coordinates.
#' See the \pkg{geomorph} package for details.
#' @references Gilbert, M.C. 2016. Impacts of habitat fragmentation on the 
#' cranial morphology of a 
#' threatened desert fish (Cyprinodon pecosensis). Masters Thesis, 
#' Western Kentucky University.
NULL

#' Plethodon comparative morphological data 
#'
#' @name PlethMorph
#' @docType data
#' @author Michael Collyer and Dean Adams
#' @keywords datasets
#' @description Data for 37 species of plethodontid salamanders.  
#' Variables include snout to vent length
#' (SVL) as species size, tail length, head length, snout to eye length, 
#' body width, forelimb length,
#' and hind limb length, all measured in mm.  A grouping variable is also 
#' included for functional guild size.  
#' The data set also includes a phylogenetic covariance matrix based on a 
#' Brownian model of evolution, to assist in 
#' generalized least squares (GLS) estimation.
#' @details The covariance matrix was estimated with the vcv.phylo function 
#' of the R package, ape, based on the tree
#' described in Adams and Collyer (2018).
#' @references Adams, D.C and Collyer, M.L. 2018. Multivariate phylogenetic 
#' anova: group-clade aggregation, biological 
#' challenges, and a refined permutation procedure. Evolution, 72: 1204-1215.
NULL

#' Simulated motion paths
#'
#' @name motionpaths
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for 
#' the analysis of phenotypic
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @keywords datasets
NULL

#####----------------------------------------------------------------------------------------------------

# HELP FUNCTIONs

#' Create a data frame for lm.rrpp analysis
#'
#' Create a data frame for lm.rrpp analysis, when covariance or distance 
#' matrices are used
#'
#' This function is not much different than \code{\link{data.frame}} but is 
#' more flexible to allow
#' distance matrices and covariance matrices to be included.  Essentially, 
#' this function creates a list,
#' much like an object of class \code{data.frame} is also a list.  However, 
#' \code{rrpp.data.frame} is
#' less concerned with coercing the list into a matrix and more concerned 
#' with matching the number of observations (n).
#' It is wise to use this function with any \code{lm.rrpp} analysis so that 
#' \code{\link{lm.rrpp}} does not have to search
#' the global environment for needed data.
#'
#' It is assumed that multiple data sets for the same subjects are in the 
#' same order.
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
    cat("\nWarning: Some geomorph.data.frame objects might not be 
        compatible with RRPP functions.")
    cat("\nIf any part of the geomorph.data.frame conatins a 3D array,")
    cat("\nconsider converting it to a matrix before attempting to 
        make an rrpp.data.frame.")
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
      if(is.matrix(dots[[i]])) {
        dots.ns[i] <- dim(dots[[i]])[[1]]
        dt.nms <- rownames(dots[[i]])
        if(dots.ns[i] == length(dots[[i]])) {
          dots[[i]] <- as.vector(dots[[i]])
          names(dots[[i]]) <- dt.nms
        }
      }
      
      if(inherits(dots[[i]], "dist")) dots.ns[i] <- attr(dots[[i]], "Size")
      if(is.data.frame(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[1]]
      if(is.vector(dots[[i]])) dots.ns[i] <- length(dots[[i]])
      if(is.factor(dots[[i]])) dots.ns[i] <- length(dots[[i]])
      if(is.logical(dots[[i]])) dots.ns[i] <- length(dots[[i]])
    }
    if(any(is.na(dots.ns))) 
      stop("Some input is either dimensionless or inappropriate for data frames")
    if(length(unique(dots.ns)) > 1) 
      stop("Inputs have different numbers of observations")
    class(dots) <- c("rrpp.data.frame")
  }
  dots
}

#####--------------------------------------------------------------\

# SUPPORT FUNCTIONS


# lm.rrpp subfunctions
# lm-like fit modified for all submodels
# general workhorse for all 'rrpp.lm' functions
# used in all 'rrpp.lm' functions

makeDF <- function(form, data, n) {
  
  if(!is.list(data)) 
    stop("\nThe data frame provide is not class rrpp.data.frame, 
         data.frame, or list, and therefore, unusable.\n", call. = FALSE)
  
  dat <- data
  class(dat) <- "list"
  
  form <- try(as.formula(form), silent = TRUE)
  if(inherits(form, "try-error"))
    stop("Formula is not coercible into a formula object.  
         Please fix the formula.\n",
         call. = FALSE)
              
  var.names <- if(length(form) == 3) all.vars(form)[-1] else if(length(form) == 2) 
      all.vars(form) else
        stop("\nFormula is not appropriately formatted.\n", call. = FALSE)
  
  dat <- dat[names(dat) %in% var.names]
  
  if(length(dat) > 0) {
    
    var.check <- sapply(seq_len(length(dat)), function(j) {
      x <- dat[[j]]
      length(x) == n
    })
    
    if(any(!var.check)) 
      stop("One or more independent variables does not match the number 
           of observations in the data.\n",
           call. = FALSE)
  }
   
  dat <- if(length(dat) == 0)  NULL else as.data.frame(dat)
  
  dat
}
  

lm.args.from.formula <- function(cl){
  
  lm.args <- list(formula = NULL, data = NULL, subset = NULL, weights = NULL,
                  na.action = na.omit, method = "qr", model = TRUE, 
                  qr = TRUE,
                  singular.ok = TRUE, contrasts = NULL, offset = NULL, tol = 1e-7)
  
  lm.nms <- names(lm.args)
  
  m1 <- match(names(cl), lm.nms)
  m2 <- match(lm.nms, names(cl))
  lm.args[na.omit(m1)] <- cl[na.omit(m2)]
  lm.args$x <- lm.args$y <- TRUE
  
  form <- lm.args$formula
  if(is.null(form))
    stop("The formula is either missing or not formatted correctly.\n", 
         call. = FALSE)
  
  Dy <- NULL
  Y <- try(eval(lm.args$formula[[2]], lm.args$data, parent.frame()),
           silent = TRUE)
  
  nms <- if(is.vector(Y)) names(Y) else if(inherits(Y, "dist")) attr(Y, "Labels") else
    if(inherits(Y, "matrix")) rownames(Y) else dimnames(Y)[[3]]
  
  if(inherits(Y, "try-error"))
    stop("Data are missing from either the data frame or global environment.\n", 
         call. = FALSE)
  
  if(is.vector(Y)) {
    Y <- as.matrix(Y)
    if(!is.null(nms)) rownames(Y) <- nms
    Dy <- NULL
  }
  
  if(inherits(Y, "matrix") || is.data.frame(Y)) {
    if(isSymmetric(Y)) {
      Dy <- Y <- as.dist(Y)
      if(any(Dy < 0)) stop("Distances in distance matrix cannot be less than 0\n",
                           call. = FALSE)
      lm.args$formula <- update(lm.args$formula, Y ~ .)
    } else Dy <- NULL
  }
  
  if(inherits(Y, "dist")) {
    if(any(Y < 0)) stop("Distances in distance matrix cannot be less than 0")
    Dy <- Y
    Y <- pcoa(Y)
  }
  
  if(is.array(Y) && length(dim(Y)) > 2) 
    stop("Data are arranged in an array rather than a matrix.  
         Please convert data to a matrix. \n", 
         call. = FALSE)
  
  if(form[[3]] == ".") {
    xs <- paste(names(lm.args$data), collapse = "+")
    form <- as.formula(noquote(c("~", xs)))
  }
  form <- update(form, Y ~.,)
  lm.args$formula <- form
  
  n <- NROW(Y)
  
  if(!is.null(lm.args$data)) {
    lm.args$data <- makeDF(form, lm.args$data, n)
    lm.args$data$Y <- as.matrix(Y)
  }
  
  if(is.null(lm.args$data)) {
    lm.args$data <- data.frame(Int = rep(1, n))
    lm.args$data$Y <- as.matrix(Y)
    lm.args$data <- lm.args$data[-1]
  }
  
  rownames(lm.args$data) <- nms
  
  f <- try(do.call(lm, lm.args), silent = TRUE)
  
  if(inherits(f, "try-error")) 
    stop("Variables or data are missing from either the data frame or 
           global environment\n", 
         call. = FALSE)
  
  Y = as.matrix(f$y)
  rownames(Y) <- nms
  out <- list(Terms = f$terms, model = f$model, Y = Y)
  if(!is.null(Dy)) {
    d <- as.matrix(Dy)
    if(nrow(d) != NROW(out$Y)) d <- d[rownames(Y), rownames(Y)]
    out$D <- as.dist(d)
  }
  
  out
}

lm.fits <- function(Terms, Y, offset = NULL, tol = 1e-7, SS.type = "I", model) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  Y <- as.matrix(Y)
  o <- if(!is.null(offset)) offset else NULL
  fit.args <- list(x = X, y = Y, offset = o, tol = tol, method = "qr", 
                   singular.ok = TRUE)
  fit <- do.call(lm.fit, fit.args)
  Pcov <- NULL
  
  X.n.k.obs <- length(X.k.obs)
  QRx <- fit$qr
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
  
  if(k > 0) {
    
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
    get.lm <- function(x) {
      fit.args$x <- as.matrix(x)
      do.call(lm.fit, fit.args)
    }
    
    reduced <- lapply(Xrs, get.lm)
    full <-lapply(Xfs, get.lm)
    
  } else {
    reduced <- full <- fit
    terms.r <- terms.f <- Terms
    Xrs <- Xfs <- X
  }
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, 
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
}

lm.fits.exchange <- function(Terms, Y, offset = NULL, tol = 1e-7, 
                             SS.type = "I", model) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  Y <- as.matrix(Y)
  o <- if(!is.null(offset)) offset else NULL
  fit.args <- list(x = X, y = Y, offset = o, tol = tol, method = "qr", 
                   singular.ok = TRUE)
  fit <- do.call(lm.fit, fit.args)
  Pcov <- NULL
  
  X.n.k.obs <- length(X.k.obs)
  QRx <- fit$qr
  X.n.k <- QRx$rank
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0, X.k))
  
  if(X.n.k < X.n.k.obs) fix <- TRUE else fix <- FALSE
  
  if(fix) Terms <- Terms[uk]
  
  term.labels <- attr(Terms, "term.labels")
  k <- length(term.labels)
  
  if(k > 0) {
    
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
    get.lm <- function(x) {
      fit.args$x <- as.matrix(x)
      do.call(lm.fit, fit.args)
    }
    
    reduced <- lapply(Xrs, get.lm)
    full <-lapply(Xfs, get.lm)
    
  } else {
    reduced <- full <- fit
    terms.r <- terms.f <- Terms
    Xrs <- Xfs <- X
  }
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, 
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
}

lm.wfits <- function(Terms,Y, w, offset = NULL, tol = 1e-7, 
                     SS.type = "I", model) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  Y <- as.matrix(Y)
  o <- if(!is.null(offset)) offset else NULL
  fit.args <- list(x = X, y = Y, w = w, offset = o, tol = tol, 
                   method = "qr", singular.ok = TRUE)
  fit <- do.call(lm.wfit, fit.args)
  Pcov <- NULL
  
  X.n.k.obs <- length(X.k.obs)
  QRx <- fit$qr
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
  
  if(k > 0) {
    
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
    get.wlm <- function(x) {
      fit.args$x <- as.matrix(x)
      do.call(lm.wfit, fit.args)
    }
    
    reduced <- lapply(Xrs, get.wlm)
    full <-lapply(Xfs, get.wlm)
    
  } else {
    reduced <- full <- fit
    terms.r <- terms.f <- Terms
    Xrs <- Xfs <- X
  }
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset,weights = w^2,
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
}


lm.wfits.exchange <- function(Terms,Y, w, offset = NULL, tol = 1e-7, 
                              SS.type = "I", model) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  Y <- as.matrix(Y)
  o <- offset
  weighted <- if(!is.null(w)) TRUE else FALSE
  offst <- if(!is.null(o)) TRUE else FALSE
  if(weighted && offst) o <- o * w
  
  fit.args <- list(x = X, y = Y, offset = o, tol = tol, 
                   method = "qr", singular.ok = TRUE)
  fitf.args <- fit.args
  fitf.args$w <- w
  fit <- do.call(lm.wfit, fitf.args)
  Pcov <- NULL
  
  X.n.k.obs <- length(X.k.obs)
  QRx <- fit$qr
  X.n.k <- QRx$rank
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0, X.k))
  
  if(X.n.k < X.n.k.obs) fix <- TRUE else fix <- FALSE
  
  if(fix) Terms <- Terms[uk]
  
  term.labels <- attr(Terms, "term.labels")
  k <- length(term.labels)
  
  if(k > 0) {
    
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
    
    
    fit.args$y <- fit.args$y * sqrt(w)
    if(!is.null(o)) fit.args$offset <- o * sqrt(w)
    get.lm <- function(x) {
      fit.args$x <- as.matrix(x * sqrt(w))
      do.call(lm.fit, fit.args)
    }
    
    reduced <- lapply(Xrs, get.lm)
    full <-lapply(Xfs, get.lm)
    
  } else {
    reduced <- full <- fit
    terms.r <- terms.f <- Terms
    Xrs <- Xfs <- X
  }
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, weights = w^2,
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
}


lm.glsfits.exchange <- function(Terms,Y, Cov, offset = NULL, tol = 1e-7, 
                                SS.type = "I", model) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  Y <- as.matrix(Y)
  Pcov <- Cov.proj(Cov, rownames(Y))
  PX <- Pcov %*% X
  PY <- Pcov %*% Y
  o <- if(!is.null(offset)) offset else NULL
  fit.args <- list(x = PX, y = PY, offset = o,  tol = tol, method = "qr", 
                   singular.ok = TRUE)
  fit <- do.call(lm.fit, fit.args)
  
  X.n.k.obs <- length(X.k.obs)
  QRx <- fit$qr
  X.n.k <- QRx$rank
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0, X.k))
  
  if(X.n.k < X.n.k.obs) fix <- TRUE else fix <- FALSE
  
  if(fix) Terms <- Terms[uk]
  
  term.labels <- attr(Terms, "term.labels")
  k <- length(term.labels)
  
  if(k > 0) {
    
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
    get.lm <- function(x) {
      fit.args$x <- Pcov %*% as.matrix(x)
      do.call(lm.fit, fit.args)
    }
    
    reduced <- lapply(Xrs, get.lm)
    full <-lapply(Xfs, get.lm)
    
  } else {
    reduced <- full <- fit
    terms.r <- terms.f <- Terms
    Xrs <- Xfs <- X
  }
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, 
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
}



lm.glsfits <- function(Terms,Y, Cov, offset = NULL, tol = 1e-7, 
                       SS.type = "I", model) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  Y <- as.matrix(Y)
  X.n.k.obs <- length(X.k.obs)
  QRx <- qr(X)
  X.n.k <- QRx$rank
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0, X.k))
  
  Pcov <- Cov.proj(Cov, rownames(Y))
  PX <- Pcov %*% X
  PY <- Pcov %*% Y
  o <- if(!is.null(offset)) offset else NULL
  fit.args <- list(x = PX, y = PY, offset = o,  tol = tol, method = "qr", 
                   singular.ok = TRUE)
  fit <- do.call(lm.fit, fit.args)
  
  
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
  
  if(k > 0) {
    
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
    
    exchanges <- lm.glsfits.exchange(Terms, Y, Cov, offset = o,
                                     tol = tol, SS.type = SS.type, model)
    fits <- lm.fits(Terms, Y, offset = o,
                    tol = tol, SS.type = SS.type, model = model)
    refit <- function(L, E, Xs) {
      lapply(1:length(L), function(j){
        f <- L[[j]]
        e <- E[[j]]
        x <- as.matrix(Xs[[j]])
        coef.j <- e$coefficients
        f$coefficients <- coef.j
        f$fitted.values <- if(length(coef.j > 0)) as.matrix(x %*% coef.j) else 
          f$fitted.values
        if(!is.null(e$effects)) f$effects <- qr.qty(e$qr, PY) else {
          Q <- qr(rep(0, NROW(X)))
          f$effects <- qr.qty(Q, PY)
        }
        f$residuals <- if(!is.null(f$fitted.values)) Y - f$fitted.values else Y
        f$qr <- if(!is.null(e$effects)) e$qr else Q
        f$qraux <- e$qraux
        f
      })
    }
    
    reduced <- refit(fits$reduced, exchanges$reduced, Xrs)
    full <- refit(fits$full, exchanges$full, Xfs)
    
  } else {
    reduced <- full <- fit
    terms.r <- terms.f <- Terms
    Xrs <- Xfs <- X
  }
  
  model.sets <- list(terms.r = terms.r, terms.f = terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, 
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
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


SS.iter <- function(exchange, ind, RRPP = TRUE, print.progress = TRUE) {
  reduced <- exchange$reduced
  full <- exchange$full
  Terms <- exchange$Terms
  trms <- attr(Terms, "term.labels")
  k <- length(trms)
  Y <- as.matrix(exchange$model[[1]])
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  w <- if(!is.null(exchange$weights)) exchange$weights else NULL
  perms <- length(ind)
  
  if(k > 0) {
    fitted <- lapply(reduced, function(x) as.matrix(x$fitted.values))
    res <- lapply(reduced, function(x) as.matrix(x$residuals))
  } else {
    fitted <- list(fitted = as.matrix(reduced$fitted.values))
    res <- list(residuals = as.matrix(reduced$residuals))
  }
  
  if(!RRPP) {
    if(k > 0) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) as.matrix(Y))
    } else {
      fitted <- list(fitted = matrix(0, n, p))
      res <- list(res = as.matrix(Y))
    }
  } 
  
  Pcov <- exchange$Pcov
  
  
  if(k > 0) {
    Qr <- lapply(reduced, function(x) if(!is.null(x$qr)) x$qr else qr(rep(0, n)))
    Qf <- lapply(full, function(x) x$qr)
  } else {
    Qr <- list(Qr = if(!is.null(reduced$qr)) reduced$qr else qr(rep(0, n)))
    Qf <- list(Qf = full$qr)
  }
  
  
  Ur <- lapply(Qr, function(x) qr.Q(x))
  Uf <- lapply(Qf, function(x) qr.Q(x))
  Urs <- lapply(Ur, function(x) Matrix(round(x, 15), sparse = TRUE))
  Ufs <- lapply(Uf, function(x) Matrix(round(x, 15), sparse = TRUE))
  Hr <- lapply(Ur, function(x) 
    forceSymmetric(Matrix(round(tcrossprod(x), 15), sparse = TRUE)))
  Hf <- lapply(Uf, function(x) 
    forceSymmetric(Matrix(round(tcrossprod(x), 15), sparse = TRUE)))
  
  # Linear model checkers
  for(i in 1:max(1, k)) {
    o.ur <- object.size(Ur[[i]])
    o.urs <- object.size(Urs[[i]])
    if(o.urs < o.ur) Ur[[i]] <- Urs[[i]]
    o.uf <- object.size(Uf[[i]])
    o.ufs <- object.size(Ufs[[i]])
    if(o.ufs < o.uf) Uf[[i]] <- Ufs[[i]]
  }
  
  for(i in 1:max(1, k)) {
    h <- Hr[[i]]
    sparse <- length(h@x)/h@Dim[1]/h@Dim[2]
    if(sparse < 0.2) Ur[[i]] <- Hr[[i]]
    h <- Hf[[i]]
    sparse <- length(h@x)/h@Dim[1]/h@Dim[2]
    if(sparse < 0.2) Uf[[i]] <- Hf[[i]]
  }
  
  Ufull <-Uf[[max(1, k)]]
  
  int <- attr(Terms, "intercept")
  Unull <- if(!is.null(Pcov)) 
    qr.Q(qr(crossprod(Pcov, rep(int, n)))) else if(!is.null(w)) 
      qr.Q(qr(rep(int, n) * sqrt(w))) else
        qr.Q(qr(rep(int, n)))
  if(!is.null(Pcov)) Y <- Pcov %*% Y
  if(!is.null(w)) Y <- Y * sqrt(w)
  
  yh0 <- fastFit(Unull, Y, n, p)
  r0 <- Y - yh0
  FR <- lapply(1:max(1, k), function(j) list(fitted = fitted[[j]], 
                                             residuals = res[[j]]))
  rrpp.args <- list(FR = FR, ind.i = NULL)
  
  rrpp <- function(FR, ind.i) {
    lapply(FR, function(x) x$fitted + x$residuals[ind.i, ])
  }
  
  if(print.progress) {
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  
  ss <- function(ur, uf, y) c(sum(crossprod(ur, y)^2), sum(crossprod(uf, y)^2), 
                              sum(y^2) - sum(crossprod(Ufull, y)^2))
  
  result <- lapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    x <-ind[[j]]
    rrpp.args$ind.i <- x
    Yi <- do.call(rrpp, rrpp.args)
    y <- yh0 + r0[x,]
    yy <- sum(y^2)
    
    if(k > 0) {
      res <- vapply(1:k, function(j){
        ss(Ur[[j]], Uf[[j]], Yi[[j]])
      }, numeric(3))
      
      SSr <- res[1, ]
      SSf <- res[2, ]
      RSS <- res[3, ]
      
      TSS <- yy - sum(crossprod(Unull, y)^2)
      TSS <- rep(TSS, k)
      SS = SSf - SSr
    } else SSr <- SSf <- SS <- RSS <- TSS <- NA
    RSS.model <- yy - sum(crossprod(Ufull, y)^2)
    if(k == 0) TSS <- RSS.model
    list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
  })
  
  SS <- matrix(sapply(result, "[[", "SS"), max(1, k), perms)
  RSS <- matrix(sapply(result, "[[", "RSS"), max(1, k), perms)
  TSS <- matrix(sapply(result, "[[", "TSS"), max(1, k), perms)
  RSS.model <- matrix(sapply(result, "[[", "RSS.model"), max(1, k), perms)
  
  if(print.progress)  close(pb)
    
  res.names <- list(if(k > 0) trms else "Intercept", 
                    c("obs", paste("iter", 1:(perms-1), sep=".")))
  dimnames(SS) <- dimnames(RSS) <- dimnames(TSS) <- 
    dimnames(RSS.model) <- res.names
  
  if(all(is.na(SS))) RSS <- SS <- NULL
  
  list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
  
}


# SS.iterPP
# workhorse for lm.rrpp, same as SS.iter, but with parallel processing
# used in lm.rrpp

SS.iterPP <- function(exchange, ind, RRPP = TRUE, print.progress = TRUE, 
                      ParCores = TRUE) {
  
  if(is.logical(ParCores)) no_cores <- detectCores() - 1 else
    no_cores <- min(detectCores() - 1, ParCores)
  
  Unix <- .Platform$OS.type == "unix"
  reduced <- exchange$reduced
  full <- exchange$full
  Terms <- exchange$Terms
  trms <- attr(Terms, "term.labels")
  k <- length(trms)
  Y <- as.matrix(exchange$model[[1]])
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  w <- if(!is.null(exchange$weights)) exchange$weights else NULL
  
  perms <- length(ind)
  
  if(k > 0) {
    fitted <- lapply(reduced, function(x) as.matrix(x$fitted.values))
    res <- lapply(reduced, function(x) as.matrix(x$residuals))
  } else {
    fitted <- list(fitted = as.matrix(reduced$fitted.values))
    res <- list(residuals = as.matrix(reduced$residuals))
  }
  
  if(!RRPP) {
    if(k > 0) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) as.matrix(Y))
    } else {
      fitted <- list(fitted = matrix(0, n, p))
      res <- list(res = as.matrix(Y))
    }
  } 
  
  Pcov <- exchange$Pcov
  
  if(k > 0) {
    Qr <- lapply(reduced, function(x) if(!is.null(x$qr)) x$qr else qr(rep(0, n)))
    Qf <- lapply(full, function(x) x$qr)
  } else {
    Qr <- list(Qr = if(!is.null(reduced$qr)) reduced$qr else qr(rep(0, n)))
    Qf <- list(Qf = full$qr)
  }
  
  Ur <- lapply(Qr, function(x) qr.Q(x))
  Uf <- lapply(Qf, function(x) qr.Q(x))
  Urs <- lapply(Ur, function(x) Matrix(round(x, 15), sparse = TRUE))
  Ufs <- lapply(Uf, function(x) Matrix(round(x, 15), sparse = TRUE))
  Hr <- lapply(Ur, function(x) 
    forceSymmetric(Matrix(round(tcrossprod(x), 15), sparse = TRUE)))
  Hf <- lapply(Uf, function(x) 
    forceSymmetric(Matrix(round(tcrossprod(x), 15), sparse = TRUE)))
  
  # Linear model checkers
  for(i in 1:max(1, k)) {
    o.ur <- object.size(Ur[[i]])
    o.urs <- object.size(Urs[[i]])
    if(o.urs < o.ur) Ur[[i]] <- Urs[[i]]
    o.uf <- object.size(Uf[[i]])
    o.ufs <- object.size(Ufs[[i]])
    if(o.ufs < o.uf) Uf[[i]] <- Ufs[[i]]
  }
  
  for(i in 1:max(1, k)) {
    h <- Hr[[i]]
    sparse <- length(h@x)/h@Dim[1]/h@Dim[2]
    if(sparse < 0.2) Ur[[i]] <- Hr[[i]]
    h <- Hf[[i]]
    sparse <- length(h@x)/h@Dim[1]/h@Dim[2]
    if(sparse < 0.2) Uf[[i]] <- Hf[[i]]
  }
  
  Ufull <-Uf[[max(1, k)]]
  
  int <- attr(Terms, "intercept")
  Unull <- if(!is.null(Pcov)) 
    qr.Q(qr(crossprod(Pcov, rep(int, n)))) else if(!is.null(w)) 
      qr.Q(qr(rep(int, n) * sqrt(w))) else
        qr.Q(qr(rep(int, n)))
  if(!is.null(Pcov)) Y <- Pcov %*% Y
  if(!is.null(w)) Y <- Y * sqrt(w)
  
  yh0 <- fastFit(Unull, Y, n, p)
  r0 <- Y - yh0
  
  if(print.progress)
    cat("\nProgress bar not available for sums of squares calculations...\n")
  
  FR <- lapply(1:max(1, k), function(j) list(fitted = fitted[[j]], 
                                             residuals = res[[j]]))
  rrpp.args <- list(FR = FR, ind.i = NULL)
  
  rrpp <- function(FR, ind.i) {
    lapply(FR, function(x) x$fitted + x$residuals[ind.i, ])
  }
  
  ss <- function(ur, uf, y) c(sum(crossprod(ur, y)^2), sum(crossprod(uf, y)^2), 
                              sum(y^2) - sum(crossprod(Ufull, y)^2))
  
  if(Unix) {
    result <- mclapply(1:perms, mc.cores = no_cores, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- yh0 + r0[x,]
      yy <- sum(y^2)
      if(k > 0) {
        res <- vapply(1:k, function(j){
          ss(Ur[[j]], Uf[[j]], Yi[[j]])
        }, numeric(3))
        
        SSr <- res[1, ]
        SSf <- res[2, ]
        RSS <- res[3, ]
        
        TSS <- yy - sum(crossprod(Unull, y)^2)
        TSS <- rep(TSS, k)
        SS = SSf - SSr
      } else SSr <- SSf <- SS <- RSS <- TSS <- NA
      RSS.model <- yy - sum(crossprod(Ufull, y)^2)
      if(k == 0) TSS <- RSS.model
      list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
    })
    
  } else {
    
    cl <- makeCluster(no_cores)
    clusterExport(cl, "ind",
                  envir=environment())
    
    result <- parLapply(cl, 1:perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      y <- yh0 + r0[x,]
      yy <- sum(y^2)
      if(k > 0) {
        res <- vapply(1:k, function(j){
          ss(Ur[[j]], Uf[[j]], Yi[[j]])
        }, numeric(3))
        
        SSr <- res[1, ]
        SSf <- res[2, ]
        RSS <- res[3, ]
        
        TSS <- yy - sum(crossprod(Unull, y)^2)
        TSS <- rep(TSS, k)
        SS = SSf - SSr
      } else SSr <- SSf <- SS <- RSS <- TSS <- NA
      RSS.model <- yy - sum(crossprod(Ufull, y)^2)
      if(k == 0) TSS <- RSS.model
      list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
    })
    stopCluster(cl)
  }
  
  
  SS <- matrix(sapply(result, "[[", "SS"), max(1, k), perms, byrow = TRUE)
  RSS <- matrix(sapply(result, "[[", "RSS"), max(1, k), perms, byrow = TRUE)
  TSS <- matrix(sapply(result, "[[", "TSS"), max(1, k), perms, byrow = TRUE)
  RSS.model <- matrix(sapply(result, "[[", "RSS.model"), max(1, k), 
                      perms, byrow = TRUE)
  
  res.names <- list(if(k > 0) trms else "Intercept", 
                    c("obs", paste("iter", 1:(perms-1), sep=".")))
  dimnames(SS) <- dimnames(RSS) <- dimnames(TSS) <- dimnames(RSS.model) <- 
    res.names
  
  if(all(is.na(SS))) RSS <- SS <- NULL
  
  list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
  
}


# anova.parts
# construct an ANOVA tablefrom random SS output
# used in lm.rrpp

anova.parts <- function(exchange, SS){
  SS.type <- exchange$SS.type
  perms <- NCOL(SS)
  reduced <- exchange$reduced
  full <- exchange$full
  Terms <- exchange$Terms
  trms <- attr(Terms, "term.labels")
  k <- length(trms)
  n <- NROW(exchange$model)
  if(k > 0) {
    QRr <- lapply(reduced, function(x) if(!is.null(x$qr)) x$qr else qr(rep(0, n)))
    QRf <- lapply(full, function(x) x$qr)
  } else {
    QRr <- if(!is.null(reduced$qr)) reduced$qr else qr(rep(0, n))
    QRf <- full$qr
  }
  
  dims <- dim(as.matrix(exchange$model[[1]]))
  n <- dims[1]; p <- dims[2]
  
  if(k > 0) {
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
    dft <- n - 1
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
    
  } else {
    RSS <- NULL
    TSS <- SS$TSS
    RSS.model <- SS$RSS.model
    SS <- NULL
    Rsq <- NULL
    MS <- NULL
    RMS <- NULL
    Fs <- NULL
    cohenf <- NULL
    df <- n - 1
    SS.type <- NULL
    
  }
  
  
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
# not used but was used in compare.models, during development.  Retained
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

# beta.boot.iter
# gets appropriate beta vectors for coefficients via bootstrap
# used in predict.lm.rrpp

beta.boot.iter <- function(fit, ind) {
  Terms <- fit$LM$Terms
  trms <- fit$LM$term.labels
  k <- length(trms)
  gls <- fit$LM$gls
  
  fitted <- if(gls) fit$LM$gls.fitted else fit$LM$fitted
  res <- if(gls) fit$LM$gls.residuals else fit$LM$residuals
  
  
  w <- fit$LM$weights
  if(!is.null(w)) weighted = TRUE else weighted = FALSE
  
  if(weighted) {
    fitted <- fitted * sqrt(w)
    res <- res * sqrt(w)
  }
  
  Y <- fitted + res
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  perms <- length(ind)
  
  Pcov <- fit$LM$Pcov
  
  rrpp.args <- list(fitted = as.matrix(fitted), 
                    residuals = as.matrix(res),
                    ind.i = NULL)
  
  rrpp <- function(fitted, residuals, ind.i) as.matrix(fitted + residuals[ind.i,])
  
  Qf <- fit$LM$QR
  Hf <- tcrossprod(solve(qr.R(Qf)), qr.Q(Qf))
  Hfs <- Matrix(round(Hf, 15), sparse = TRUE)
  if(object.size(Hfs) < object.size(Hf)) Hf <- Hfs
  
  betas <- lapply(1:perms, function(j){
    x <-ind[[j]]
    rrpp.args$ind.i <- x
    Yi <- do.call(rrpp, rrpp.args)
    if(!is.null(Pcov)) Yi <- crossprod(Pcov, Yi)
    as.matrix(Hf %*% Yi)
  })
  betas
  
}

# beta.iter
# gets appropriate beta vectors for random permutations in lm.rrpp
# generates distances as statistics for summary


beta.iter <- function(exchange, ind, RRPP = TRUE, print.progress = TRUE) {
  reduced <- exchange$reduced
  full <- exchange$full
  Terms <- exchange$Terms
  trms <- attr(Terms, "term.labels")
  k <- length(trms)
  Y <- as.matrix(exchange$model[[1]])
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  perms <- length(ind)
  if(k > 0) {
    fitted <- lapply(reduced, function(x) as.matrix(x$fitted.values))
    res <- lapply(reduced, function(x) as.matrix(x$residuals))
  } else {
    fitted <- list(fitted = as.matrix(reduced$fitted.values))
    res <- list(residuals = as.matrix(reduced$residuals))
  }
  
  if(!RRPP) {
    if(k > 0) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) as.matrix(Y))
    } else {
      fitted <- list(fitted = matrix(0, n, p))
      res <- list(res = as.matrix(Y))
    }
  } 
  
  o <- exchange$offset
  if(!is.null(o)) offset = TRUE else offset = FALSE
  
  Pcov <- exchange$Pcov
  
  FR <- lapply(1:max(1, k), function(j) list(fitted = fitted[[j]], 
                                             residuals = res[[j]]))
  rrpp.args <- list(FR = FR, ind.i = NULL, offset = offset, o = o)
  
  rrpp <- function(FR, ind.i, offset, o) {
    if(offset) lapply(FR, function(x) x$fitted + x$residuals[ind.i, ] - o) else 
      lapply(FR, function(x) x$fitted + x$residuals[ind.i, ])
  }
  
  if(print.progress){
    cat(paste("\nCoefficients estimation:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  
  Qf <- if(k > 0) lapply(full, function(x) x$qr) else list(Qf = full$qr)
  Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
  Hfs <- lapply(Hf, function(x) Matrix(round(x, 15), sparse = TRUE))
  for(i in 1:max(1,k)) {
    if(object.size(Hfs[[i]]) < object.size(Hf[[i]]))
      Hf[[i]] <- Hfs[[i]]
  }
  

  betas <- lapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    x <-ind[[j]]
    rrpp.args$ind.i <- x
    Yi <- do.call(rrpp, rrpp.args)
    Map(function(h, y)  as.matrix(h %*% y), Hf, Yi) 
  
  })

beta.mats <- betas

if(k > 0) {
  
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
    if(is.vector(result)) result <- matrix(result, 1, kk)
    rownames(result) <- rownames(b[[1]])
    colnames(result) <- c("obs", paste("iter", seq(1,(perms-1),1), sep = "."))
    result
  })
  
  names(beta.mats) <- names(beta.mat.d) <- trms
  
  if(k == 1) d.stitched <- beta.mat.d[[1]]
  
  if(k > 1) {
    d.stitched <- beta.mat.d[[k]]
    
    beta.match <- lapply(1:k, function(j){
      cr <- rownames(reduced[[j]]$coefficients)
      cf <- rownames(full[[j]]$coefficients)
      cf[!(cf %in% cr)]
    })
    
    for(i in 1:(k-1)){
      beta.check <- beta.match[[i]]
      d.check <- beta.mat.d[[i]]
      target <- which(rownames(d.stitched) %in% beta.check)
      d.stitched[target, ] <- d.check[rownames(d.check) %in% beta.check,]
    }
  }
} else {
  beta.mat.d <- d.stitched <- NULL
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


beta.iterPP <- function(exchange, ind, RRPP = TRUE, print.progress = TRUE, 
                        ParCores =  TRUE) {
  
  if(is.logical(ParCores)) no_cores <- detectCores() - 1 else
    no_cores <- min(detectCores() - 1, ParCores)
  
  Unix <- .Platform$OS.type == "unix"
  reduced <- exchange$reduced
  full <- exchange$full
  Terms <- exchange$Terms
  trms <- attr(Terms, "term.labels")
  k <- length(trms)
  Y <- as.matrix(exchange$model[[1]])
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  perms <- length(ind)
  if(k > 0) {
    fitted <- lapply(reduced, function(x) as.matrix(x$fitted.values))
    res <- lapply(reduced, function(x) as.matrix(x$residuals))
  } else {
    fitted <- list(fitted = as.matrix(reduced$fitted.values))
    res <- list(res = as.matrix(reduced$residuals))
  }
  
  if(!RRPP) {
    if(k > 0) {
      fitted <- lapply(fitted, function(.) matrix(0, n, p))
      res <- lapply(res, function(.) as.matrix(Y))
    } else {
      fitted <- list(fitted = matrix(0, n, p))
      res <- list(residuals = as.matrix(Y))
    }
  } 
  
  o <- exchange$offset
  if(!is.null(o)) offset = TRUE else offset = FALSE
  
  Pcov <- exchange$Pcov
  
  FR <- lapply(1:max(1, k), function(j) list(fitted = fitted[[j]], 
                                             residuals = res[[j]]))
  rrpp.args <- list(FR = FR, ind.i = NULL, offset = offset, o = o)
  
  rrpp <- function(FR, ind.i, offset, o) {
    if(offset) lapply(FR, function(x) x$fitted + x$residuals[ind.i, ] - o) else 
      lapply(FR, function(x) x$fitted + x$residuals[ind.i, ])
  }
  
  if(print.progress){
    cat("\nProgress bar not available for coefficients estimation...\n")
  }
  
  Qf <- if(k > 0) lapply(full, function(x) x$qr) else list(Qf = full$qr)
  Hf <- lapply(Qf, function(x) tcrossprod(solve(qr.R(x)), qr.Q(x)))
  Hfs <- lapply(Hf, function(x) Matrix(round(x, 15), sparse = TRUE))
  for(i in 1:max(1,k)) {
    if(object.size(Hfs[[i]]) < object.size(Hf[[i]]))
      Hf[[i]] <- Hfs[[i]]
  }
  
  if(Unix) {
    betas <- mclapply(1:perms, mc.cores = no_cores, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      Map(function(h, y)  as.matrix(h %*% y), Hf, Yi) 
      
    })
    
  } else {
    
    cl <- makeCluster(no_cores)
    clusterExport(cl, "rrpp.args")
    
    betas <- parLapply(cl, 1:perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      Yi <- do.call(rrpp, rrpp.args)
      Map(function(h, y)  h %*% y, Hf, Yi) 
      
    })
    stopCluster(cl)
  }

  
  beta.mats <- betas
  
  if(k > 0) {
    
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
      if(is.vector(result)) result <- matrix(result, 1, kk)
      rownames(result) <- rownames(b[[1]])
      colnames(result) <- c("obs", paste("iter", seq(1,(perms-1),1), sep = "."))
      result
    })
    
    names(beta.mats) <- names(beta.mat.d) <- trms
    
    if(k == 1) d.stitched <- beta.mat.d[[1]]
    
    if(k > 1) {
      d.stitched <- beta.mat.d[[k]]
      
      beta.match <- lapply(1:k, function(j){
        cr <- rownames(reduced[[j]]$coefficients)
        cf <- rownames(full[[j]]$coefficients)
        cf[!(cf %in% cr)]
      })
      
      for(i in 1:(k-1)){
        beta.check <- beta.match[[i]]
        d.check <- beta.mat.d[[i]]
        target <- which(rownames(d.stitched) %in% beta.check)
        d.stitched[target, ] <- d.check[rownames(d.check) %in% beta.check,]
      }
    }
  } else {
    beta.mat.d <- d.stitched <- NULL
  }
  
  out <- list(random.coef = beta.mats,
              random.coef.distances = d.stitched)
  out
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
    if(!inherits(error, "character")) 
      stop("The error description is illogical.  It should be a string of character 
           values matching ANOVA terms.",
                                           call. = FALSE)
    kk <- length(error)
    if(kk != k) 
      stop("The error description should match in length the number of ANOVA terms 
           (not including Residuals)",
                     call. = FALSE)
    MSEmatch <- match(error, c(trms, "Residuals"))
    if(any(is.na(MSEmatch))) 
      stop("At least one of the error terms is not an ANOVA term",
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
      Z <- apply(Z, 1, effect.size)
      } else P.val <- NULL
    
    Residuals <- c(df[k+1], RSS[[1]], RSS[[1]]/df[k+1], 
                   RSS[[1]]/TSS[[1]], rep(NA, 3))
    Total <- c(df[k+2], TSS[[1]], rep(NA, 5))
    tab <- data.frame(Df=df[1:k], SS=SS, MS = MS, Rsq = Rsq, 
                      F = Fs, Z = Z, P.val = P.val)
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
  
  if(refModel$LM$gls) {
    X <- if(!is.null(refModel$LM$Pcov)) crossprod(refModel$LM$Pcov, 
                                                  refModel$LM$X) else
      refModel$LM$X * sqrt(refModel$LM$weights)
  } else X <- refModel$LM$X
  
  if(refModel$LM$gls) {
    Y <- if(!is.null(refModel$LM$Pcov)) crossprod(refModel$LM$Pcov, 
                                                  refModel$LM$Y) else
      refModel$LM$Y * sqrt(refModel$LM$weights)
  } else Y <- refModel$LM$Y
  
  B <- if(refModel$LM$gls) refModel$LM$gls.coefficients else 
    refModel$LM$coefficients
  U <- as.matrix(qr.Q(refModel$LM$QR))
  n <- refModel$LM$n
  p <- refModel$LM$p
  Yh <- as.matrix(fastFit(U, Y, n, p))
  R <- as.matrix(Y) - Yh
  
  K <- length(lm.list)
  Ulist <- lapply(1:K, function(j){
    m <- lm.list[[j]]
    as.matrix(qr.Q(m$LM$QR))
  })
  
  if(print.progress){
    if(K > 1)
    cat(paste("\nSums of Squares calculations for", K, "models:", 
              perms, "permutations.\n")) else
      cat(paste("\nSums of Squares calculations for", K, "model:", 
                perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+5, initial = 0, style=3)
  }

  int <- attr(refModel$LM$Terms, "intercept")
  if(refModel$LM$gls) {
    int <- if(!is.null(refModel$LM$Pcov))  crossprod(refModel$LM$Pcov, 
                                                     rep(int, n)) else
      sqrt(refModel$LM$weights)
  } else int <- rep(int, n)
  
  U0 <- as.matrix(qr.Q(qr(int)))
  yh0 <- as.matrix(fastFit(U0, Y, n, p))
  r0 <- as.matrix(Y) - yh0
  
  rY <- function(ind.i) Yh + R[ind.i,]
  rY0 <- function(ind.i) yh0 + r0[ind.i,]
  
  RSS <- function(ind.i, U0, U, Ul, K, n, p, Y, yh0, r0) {
    y <- as.matrix(rY(ind.i))
    y0 <- as.matrix(rY0(ind.i))
    rss0  <- sum(y^2) - sum(crossprod(U, y)^2)
    rss00 <- sum(y0^2) - sum(crossprod(U0, y)^2)
    
    rss <- lapply(1:K, function(j){
      u <- Ul[[j]]
      sum(y^2) - sum(crossprod(u, y)^2)
    })
    tss <- sum(y0^2) - sum(crossprod(U0, y)^2)
    
    RSSp <- c(rss0, unlist(rss), tss)
    RSSp
  }
  
  rss.list <- list(ind.i = NULL, U0 = U0, U = U, 
                   Ul = Ulist, K = K, n = n , p = p, Y = Y,
                   yh0 = yh0, r0 = r0)
  
  RSSp <- sapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    rss.list$ind.i <- ind[[j]]
    do.call(RSS, rss.list)
  })
  
  RSSy <- as.matrix(RSSp[nrow(RSSp),])
  RSSp <- as.matrix(RSSp[-(nrow(RSSp)),])
  RSSy <- as.matrix(matrix(RSSy, nrow(RSSp), perms, byrow = TRUE))
  
  fit.names <- c(refModel$call[[2]], lapply(1:K, function(j) 
    lm.list[[j]]$call[[2]]))
  rownames(RSSp) <- rownames(RSSy) <- fit.names
  
  SS <- rep(RSSp[1,], each = K + 1) - RSSp
  
  Rsq <-  SS / RSSy
  
  dfe <- n - c(object$LM$QR$rank, unlist(lapply(1:K, 
                                    function(j) lm.list[[j]]$LM$QR$rank)))
  df <- dfe[1] - dfe
  df[1] <- 1
  
  MS <- SS/rep(df, perms)
  MSE <- RSSp/matrix(rep(dfe, perms), length(dfe), perms)
  Fs <- MS/MSE
  
  SS[which(zapsmall(SS) == 0)] <- 1e-32
  MS[which(zapsmall(MS) == 0)] <- 1e-32
  Rsq[which(zapsmall(Rsq) == 0)] <- 1e-32
  Fs[which(zapsmall(Fs) == 0)] <- 1e-32
  
  if(effect.type == "SS") {
    Pvals <- apply(SS, 1, pval)
    Z <- apply(SS, 1, effect.size)
  } else   if(effect.type == "MS") {
    Pvals <- apply(MS, 1, pval)
    Z <- apply(MS, 1, effect.size)
  } else   if(effect.type == "Rsq") {
    Pvals <- apply(Rsq, 1, pval)
    Z <- apply(Rsq, 1, effect.size)
  } else{
    Pvals <- apply(Fs, 1, pval)
    Z <- apply(Fs, 1, effect.size)
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
  tab <-  rbind(tab, c(n-1, NA, RSSy[1], NA, NA, NA, NA, NA, NA, NA))
  rownames(tab)[NROW(tab)] <- "Total"
  
  if(effect.type == "SS") p.type <- "Pr(>SS)" else
    if(effect.type == "MS") p.type <- "Pr(>MS)" else
      if(effect.type == "Rsq") p.type <- "Pr(>Rsq)" else
        if(effect.type == "cohenf") p.type <- "Pr(>cohenf)" else 
          p.type <- "Pr(>F)" 
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
# gets the LS means for groups from a lm.rrpp fit, 
# after constaining covariates to mean values
# used in pairwise
getLSmeans <- function(fit, g){
  k <- length(fit$LM$term.labels)
  n <- fit$LM$n
  beta <- fit$LM$random.coef[[k]]
  dat <- fit$LM$data
  covCheck <- sapply(dat, class)
  for(i in 1:length(covCheck)) if(covCheck[i] == "numeric") 
    dat[[i]] <- mean(dat[[i]])
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
  names(means) <- c("obs", paste("iter", 1:(length(means) - 1), sep = "."))
  means
}


#' Support function for RRPP
#'
#' Calculate vector correlations for a matrix (by rows).  
#' Used for pairwise comparisons.
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


logL <- function(fit, tol = NULL, pc.no = NULL){
  if(is.null(tol)) tol = 0
  n <- fit$LM$n
  p <- fit$LM$p.prime
  k <- if (!is.null(pc.no)) {
    stopifnot(length(pc.no) == 1, is.finite(pc.no), as.integer(pc.no) > 0)
    min(as.integer(pc.no), n, p)
  } else min(n, p)
  X <- as.matrix(fit$LM$X)
  Y <- as.matrix(fit$LM$Y)
  rdf <- fit$LM$data
  gls <- fit$LM$gls
  w <- fit$LM$weights
  Pcov <- fit$LM$Pcov
  
  if(gls){
    if(!is.null(Pcov)) Sig <- crossprod(Pcov %*% fit$LM$gls.residuals)/n else
      Sig <- crossprod(fit$LM$gls.residuals * sqrt(w))/n
    
  }  else {
    
    Sig <- crossprod(fit$LM$residuals) /n
    
  }
  
  s <- svd(Sig, nu = 0, nv = k)
  sdev <- s$d/sqrt(max(1, n - 1))
  rank <- min(sum(sdev > (sdev[1L] * tol)), k)
  if(rank < k){
    pr <- seq_len(rank)
    s$v <- s$v[, pr, drop = FALSE]
  }
  
  Yc <- if(fit$LM$gls) scale(Y, scale = F, center = fit$LM$gls.mean) else
    center(Y)
  
  if(p > 1) P <- as.matrix(Yc %*% s$v) else P <- Yc
  
  fit$LM$data$Y <- P
  form <- formula(fit$LM$Terms)
  form <- update(form, Y ~ .)
  
  if(!is.null(w)) {
    pfit <- lm.rrpp(form, print.progress = FALSE, 
                    Cov = fit$LM$Cov, data = fit$LM$data, 
                    iter = 0, weights = w )
  } else {
    pfit <- lm.rrpp(form, print.progress = FALSE, 
                    Cov = fit$LM$Cov, data = fit$LM$data, 
                    iter = 0)
  }
  
  if(gls){
    if(!is.null(Pcov)) Sig <- crossprod(Pcov %*% pfit$LM$gls.residuals)/n else
      Sig <- crossprod(pfit$LM$gls.residuals * sqrt(w))/n
    if(kappa(Sig) > 1e10) Sig <- RiReg(Sig, pfit$LM$gls.residuals)
    
  }  else {
    
    Sig <- crossprod(pfit$LM$residuals) /n
    if(kappa(Sig) > 1e10) Sig <- RiReg(Sig, pfit$LM$residuals)
    
  }
  
  if(gls) {
    if(!is.null(Pcov)) {
      ll <- -0.5*(n * rank + n * determinant(Sig, logarithm = TRUE)$modulus[1] + 
                    rank * determinant(pfit$LM$Cov, logarithm = TRUE)$modulus[1] + 
                    n * rank * log(2*pi))
    } else {
      ll <- -0.5*(n * rank + n * determinant(Sig, logarithm = TRUE)$modulus[1] - 
                    rank * sum(log(w)) + n * rank * log(2*pi))
      
    }
  } else 
    ll <- -0.5*(n * rank + n * determinant(Sig, logarithm = TRUE)$modulus[1] + 
                  n * rank * log(2*pi))
  
  list(logL = ll, rank = rank)
  
}

cov.trace <- function(fit) {
  n <- fit$LM$n
  p <- fit$LM$p.prime
  if(fit$LM$gls){
    if(!is.null(fit$LM$Pcov)) Sig <- crossprod(fit$LM$Pcov %*% 
                                                 fit$LM$gls.residuals)/n else
      Sig <- crossprod(fit$LM$gls.residuals * sqrt(fit$LM$weights))/n
    
  }  else {
    
    Sig <- crossprod(fit$LM$residuals) /n
    
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

# lda.prep
# sets up lm.rrpp object to use in lda
# used in: lda.rrpp
lda.prep <- function(fit, tol = 1e-7, PC.no = NULL, newdata = NULL){
  dat <- fit$LM$data
  gls <- fit$LM$gls
  w <- fit$LM$weights
  Pcov <- fit$LM$Pcov
  
  dat.class <- sapply(dat, class)
  fac.list <- which(dat.class == "factor")
  if(length(fac.list) == 0) group <- factor(rep(1, n)) else {
    datf <- dat[fac.list]
    group <- factor(apply(datf, 1, paste, collapse = "."))
  }
  
  
  newY <- function(fit, group, w = NULL, Pcov=NULL) {
    Xg <- model.matrix(~group)
    Yg <- fit$LM$Y
    res <- if(gls) fit$LM$gls.residuals else fit$LM$residuals
    
    if(!is.null(w)) {
      nfit <- lm.fit(Xg * sqrt(w), Yg * sqrt(w))
      fitted <- Yg * sqrt(w) - res * sqrt(w)
      nY <- (fitted + nfit$residuals)/sqrt(w)
    } else if(!is.null(Pcov)) {
      PY <- crossprod(Pcov, Yg)
      nfit <- lm.fit(crossprod(Pcov, Xg), PY)
      fitted <- PY - crossprod(Pcov, res)
      nY <- crossprod(fast.solve(Pcov), fitted + nfit$residuals)
    } else {
      nfit <- lm.fit(Xg, Yg)
      nY <- fit$LM$fitted + nfit$res
    }
    nY
  }
  
  Yn <- newY(fit, group, w, Pcov)
  dims <- dim(Yn)
  n <- dims[1]
  p <- dims[2]
  
  Yc <- if(gls) Yn - matrix(fit$LM$gls.mean, n, p, byrow = TRUE) else
    center(Yn)
  
  V <- crossprod(Yc)/n
  if(!is.null(Pcov)) V <- crossprod(crossprod(Pcov, Yc))/n
  if(!is.null(w)) V <-  crossprod(Yc * sqrt(w))/n
  
  if(gls) s <- svd(V)
  
  PCA <- prcomp(Yn, tol = tol)
  if(gls) {
    PCA$rotation <- s$v
    PCA$sdev <- sqrt(s$d)
    PCA$x <- Yc %*% s$v
  }
  
  d <- zapsmall(PCA$sdev^2)
  d <- 1:length(d[d > 0])
  
  if(!is.null(PC.no)) {
    if(PC.no < length(d)) d <- 1:PC.no
  }
  
  PCA$rotation <- PCA$rotation[, d]
  PCA$sdev <- PCA$sdev[d]
  PCA$x <- as.matrix(PCA$x[, d])
  colnames(PCA$x) <- colnames(PCA$rotation) <- paste("PC", d, sep = "")
  
  Yn <- PCA$x
  
  if(!is.null(newdata)) {
    Yt <- newdata
    if(NCOL(newdata) != p) 
      stop("\nNumber of variables in newdata does not match the number 
           for the data in lm.rrpp fit.\n",
           call. = FALSE)
    Ytc <- if(gls) Yt - matrix(fit$LM$gls.mean, NROW(Yt), p, byrow = TRUE) else
      center(Yt)
    Yt <- Ytc %*% PCA$rotation
    
  } else Yt <- NULL
  
  list(Yn = Yn, Yt = Yt, group = group, gls = gls)
  
}

