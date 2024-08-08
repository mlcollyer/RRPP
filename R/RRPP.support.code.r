#' @name RRPP-package
#' @aliases RRPP
#' @title Linear Model Evaluation with Randomized Residual Permutation Procedures
#' @author Michael Collyer and Dean Adams
#' @return Key functions for this package:
#' \item{\code{\link{lm.rrpp}}}{Fits linear models, using RRPP.
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
#' @import ggplot2 
#' @import Matrix
#' @importFrom ape multi2di.phylo
#' @importFrom ape root.phylo
#' @importFrom ape collapse.singles
#' @importFrom stats na.omit anova as.dist as.formula cmdscale coef delete.response dist formula lm
#' lm.fit lm.wfit loess logLik model.frame.default model.matrix optimise prcomp qnorm quantile 
#' resid sd spline var model.frame
#' @importFrom graphics abline arrows axis legend lines par plot.default points text title
#' @importFrom utils combn object.size setTxtProgressBar txtProgressBar
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
#' @docType package
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
#' included for functional guild size.  A variable for species names is also 
#' included.
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

#' Simulated fish data for measurement error analysis
#'
#' @name fishy
#' @docType data
#' @author Michael Collyer
#' @references Collyer, M.L, and D.C. Adams. 2024. 
#' Interrogating random and systematic measurement 
#' error in morphometric data. Evolutionary Biology. 
#' In press.
#' @details Data as simulated in Collyer and Adams (2024),
#' resembling a fish shape, comprising Procrustes coordinates
#' for 11 anatomical landmarks.  Data represent 120 
#' configurations for 60 subjects, each with two replicates of
#' measurement.  The 60 subjects represent 20 subjects each 
#' from three groups.
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
#' fit <- lm.rrpp(d ~ x, data = rdf, iter = 99)
#' summary(fit)

rrpp.data.frame<- function(...) {
  dots <- list(...)
  if(length(dots) == 1 && is.data.frame(dots[[1]])) {
    dots <- dots[[1]]
    class(dots) <- "rrpp.data.frame"
  } else if(length(dots) == 1 && inherits(dots[[1]], "geomorph.data.frame")) {
    dots <- dots[[1]]
    
    warning(
      paste(
        "\nThis is not an error!  It is a friendly warning.\n",
        "\nSome geomorph.data.frame objects might not be compatible with RRPP functions.",
        "\nIf any part of the geomorph.data.frame conatins a 3D array,",
        "consider converting it to a matrix before attempting to make an rrpp.data.frame\n",
        "Use suppressWarnings() to turn off these warnings.\n\n", sep = " "),
      noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 

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
      if(inherits(dots[[i]], "matrix")) {
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

# Parallel.setup
# Function to help parallel library adjust to settings
# for parallel processing

Parallel.setup <- function(Parallel){
  
  ParLog <- is.logical(Parallel)
  usecluster <- inherits(Parallel, "cluster")
  if(usecluster){
    cluster <- Parallel
    ParLog <- Parallel <- TRUE
  } else cluster <- FALSE
  ParCores <- NULL
  
  if(is.numeric(Parallel)) {
    ParCores <- Parallel
    ParLog <- TRUE
    Parallel <- TRUE
  }
  
  if(ParLog && is.null(ParCores)) {
    ParCores <- detectCores() - 1
    if(usecluster) ParCores <- length(cluster)
    if(!Parallel) ParCores <- 1
  }
  
  if(is.numeric(ParCores)) {
    if(ParCores > detectCores() - 1) ParCores <- detectCores() - 1
  }
  
  Unix <- .Platform$OS.type == "unix"
  forking <- Unix && !usecluster
  
  if(is.null(ParCores)) ParCores <- 1
  
  if(ParCores == 1) {
    ParLog <- FALSE
    forking <- FALSE
    usecluster <- FALSE
    cluster <- NULL
  }
    
  if(ParCores > 1) {
    if(!Unix && !usecluster)
      cluster <- makeCluster(ParCores)
    if(Unix && usecluster)
      Unix <- FALSE
  }
  
  list(Parallel = ParLog,
       Unix = Unix, 
       forking = forking, 
       ParCores = ParCores, 
       usecluster = usecluster,
       cluster = cluster)
  
}

# lm.rrpp subfunctions
# lm-like fit modified for all submodels
# general workhorse for all 'lm.rrpp' functions

get.names.from.list <- function(L) {
  temp <- lapply(L, get.names)
  check <- which(!sapply(temp, is.null))
  temp <- if(length(check) > 0) temp[check] else NULL
  nms <- if(!is.null(temp)) temp[[1]] else NULL
  nms
}

add.names <- function(Y, nms) {
  if(is.vector(Y)) names(Y) <- nms
  if(inherits(Y, "matrix")) rownames(Y) <- nms
  if(inherits(Y, "dist")) attr(Y, "Labels") <- nms
  Y
}

makeDF <- function(form, data, n, nms) {
  
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
      nn <- if(is.matrix(x)) nrow(x) else length(x)
      nn == n
    })
    
    if(any(!var.check)) 
      stop("One or more independent variables does not match the number 
           of observations in the data.\n",
           call. = FALSE)
  }
   
  dat <- if(length(dat) == 0)  NULL else as.data.frame(dat)
  
  if(!is.null(dat)) rownames(dat) <- nms
  
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
  nmsY <- get.names(Y)
  
  if(inherits(Y, "try-error"))
    stop("Data are missing from either the data frame or global environment.\n", 
         call. = FALSE)
  
  if(is.vector(Y)) {
    Y <- as.matrix(Y)
    Dy <- NULL
  }
  
  if(inherits(Y, "matrix") || is.data.frame(Y)) {
    if(isSymmetric(as.matrix(Y))) {
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
    nmsDF <- if(inherits(lm.args$data, "data.frame"))  
      attr(lm.args$data, "row.names") else 
        get.names.from.list(lm.args$data)

    lm.args$data <- makeDF(form, lm.args$data, n, nms = NULL)
    
  }
  
  if(is.null(lm.args$data)) {
    lm.args$data <- data.frame(Int = rep(1, n))
    lm.args$data$Y <- as.matrix(Y)
    lm.args$data <- lm.args$data[-1]
  }
  
  lm.args$data$Y <- Y
  
  
  model <- try(model.frame(form, data = lm.args$data),
               silent = TRUE)
  Terms <- try(attr(model, "terms"),
               silent = TRUE)
  
  if(inherits(model, "try-error") || inherits(Terms, "try-error"))
  stop("Variables or data might be missing from either the data frame or 
           global environment, or a linear model fit just does not work...\n", 
       call. = FALSE)
  
  Y <- add.names(Y, nmsY)
  out <- list(Terms = Terms, model = model, Y = Y)
  if(!is.null(Dy)) {
    d <- as.matrix(Dy)
    if(nrow(d) != NROW(out$Y)) d <- d[rownames(Y), rownames(Y)]
    out$D <- as.dist(d)
  }
  
  out
}

.getTerms <- function(fit = NULL, Terms = NULL, SS.type = NULL) {
  if(is.null(Terms)) Terms <-  fit$LM$Terms
  if(is.null(SS.type)) SS.type <- fit$ANOVA$SS.type
  if(length(attr(Terms, "factors")) == 0)
    SS.type <- "I"
  trms <- attr(Terms, "term.labels")
  k <- length(trms)
  mod.k <- if(k > 0) c(0, seq(1, k, 1)) else 0
  
  if(!is.null(fit) && SS.type == "Within-subject II") {
    SS.type <- "IIws"
    WS <- TRUE
  } else WS <- FALSE
  
  if(k > 0) {
    if(SS.type == "III"){
      k3 <- mod.k[-1]
      modf <- lapply(as.list(k3), function(.) Terms)
      modr <- lapply(as.list(k3), function(j) Terms[-j])
    } else if(SS.type == "II" || SS.type == "IIws"){
      k2 <- mod.k[-1]
      fac <- crossprod(attr(Terms, "factor"))
      modr <- lapply(as.list(k2), function(j){
        ind <- as.logical(ifelse(fac[j,] < fac[j,j], 1, 0))
        Terms[ind]
      })
      modf <- lapply(as.list(k2), function(j){
        ind <- ifelse(fac[j,] < fac[j,j], 1, 0)
        ind[j] <- 1
        ind <- as.logical(ind)
        Terms[ind]
      })
      
      
    } else {
      kf <- mod.k[-1]
      kr <- mod.k[-(max(mod.k) + 1)]
      modf <- lapply(as.list(kf), function(j) Terms[1:j])
      modr <- lapply(as.list(kr), function(j) Terms[0:j])
    }
    
    names(modf) <- names(modr) <- trms
    
    if(WS) {
      if(fit$subTest){
        trms.change <- which(trms == fit$subjects.var)
        modf[[trms.change]] <- Terms
        modr[[trms.change]] <- Terms[!trms %in% fit$subjects.var]
      }
    }
    
  } else {
    modr <- modf <- list("Intercept" = Terms)
  }
  list(terms.r = modr, terms.f = modf)
}

LM.fit <- function(x, y, offset = NULL, tol = 1e-07) {
  if(inherits(x, "matrix")) 
    x.s <- Matrix(x, sparse = TRUE) else {
      x.s <- x
      x <- as.matrix(x.s)
    }
  osx <- length(x)
  osxs <- length(x.s@x)
  X <- if(osx < osxs) x else x.s
  x <- x.s <- NULL
  Q <- QRforX(X, tol = tol)
  if(!is.null(offset)) y <- y - offset
  dims <- dim(y)
  n <- dims[1]
  p <- dims[2]
  fitted.values <- fastFit (Q$Q, y, n, p)
  residuals <- y - fitted.values
  if(!is.null(offset)) fitted.values <- fitted.values + offset
  list(qr = Q, fitted.values = as.matrix(fitted.values),
       residuals = as.matrix(residuals))
}

removeRedundant <- function(X){
  if(NCOL(X) > 1){
    QR <- QRforX(X, returnQ = FALSE,
                 reduce = TRUE, reQR = FALSE)
    X <- QR$X
  }
  as.matrix(X)
}

getPivot <- function(Xf, Xr){
  nf <- colnames(Xf)
  nr <- colnames(Xr)
  k <- seq(1, NCOL(Xf))
  p <- match(nr, nf)[1:NCOL(Xr)]
  k <- which(!k %in% p)
  c(p, k)
}

getXs <- function(Terms, Y, SS.type, tol = 1e-7,
                  model, subjects.term = NULL) {
  
  X <- model.matrix(Terms, data = model)
  X.k <- X.k.obs <- attr(X, "assign")
  X.n.k.obs <- length(X.k.obs)
  Xred <- removeRedundant(X)
  X.rank <- NCOL(Xred)
  pivot <- getPivot(X, Xred)
  uk <- unique(c(0, X.k))
  term.labels <- attr(Terms, "term.labels")
  k <- length(term.labels)
  fix <- (X.rank < X.n.k.obs) 
  delete <- NULL
  if(fix) {
    Terms <- Terms[uk]
    warning(
      paste(
        "\nThis is not an error!  It is a friendly warning.\n",
        "\nBecause variables in the linear model are redundant,",
        "\nthe linear model design has been truncated (via QR decomposition).",
        "\nOriginal X columns:", X.n.k.obs,
        "\nFinal X columns (rank):", X.rank,
        "\nCheck coefficients or degrees of freedom in ANOVA to see changes.\n",
        "\nUse suppressWarnings() to turn off these warnings. \n\n", sep = " "),
      noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
  } 
  
  if(fix) delete <- pivot[-(1:X.rank)]
  
  if(k > 0){
    if(SS.type == "III"){
      fix <- !(length(delete) == 0)
      Xfs <- lapply(2:length(uk), function(j)  if(fix) X[, -delete] else X)
      Xrs <- lapply(2:length(uk), function(j){
        rmove <- c(which(X.k == (j - 1)), delete)
        as.matrix(X[, - rmove])
      })
      
    } else if(SS.type == "II" || SS.type == "IIws"){
      fac <- as.matrix(crossprod(attr(Terms, "factor")))
      
      Xrs <- lapply(1:NROW(fac), function(j){
        index <- fac[, j]
        m <- max(index)
        ind <-  as.logical(ifelse(c(0, index) < m, 1, 0))
        rmove <- which(!X.k %in% uk[ind])
        rmove <- unique(c(rmove, delete))
        as.matrix(X[, -rmove])
      })
      
      Xfs <- lapply(1:NROW(fac), function(j){
        index <- fac[, j]
        m <- max(index)
        ind <-  as.logical(ifelse(c(0, index) < m, 1, 0))
        ind[j + 1] <- TRUE
        rmove <- which(!X.k %in% uk[ind])
        rmove <- unique(c(rmove, delete))
        fix <- !(length(rmove) == 0)
        if(fix) as.matrix(X[, -rmove]) else X
      })
      
    } else {
      Xrs <- lapply(2:length(uk), function(j){
        rmove <- which(!X.k %in% uk[1:(j - 1)])
        rmove <- unique(c(rmove, delete))
        as.matrix(X[, -rmove])
      })
      
      Xfs <- lapply(2:length(uk), function(j){
        rmove <- which(!X.k %in% uk[1:j])
        rmove <- unique(c(rmove, delete))
        fix <- !(length(rmove) == 0)
        if(fix) as.matrix(X[, -rmove]) else X
      })
    }
    names(Xrs) <- names(Xfs) <- term.labels
  } else {
    Xrs <- Xfs <- list(Intercept = X)
  }
  
  if(!is.null(subjects.term)) {
    Xfs[[subjects.term]] <- if(fix) X[, -delete] else X
    rmove <- c(which(X.k == subjects.term), delete)
    Xrs[[subjects.term]] <- as.matrix(X[, - rmove])
  }
  
  list(Xrs = Xrs, Xfs = Xfs)
}


lm.rrpp.fit <- function(x, y, Pcov = NULL, w = NULL, offset = NULL, tol = 1e-07){
  
  getGLSlm<- function(x, y, Pcov, offset = offset, method = "qr", tol = tol){
    PX <- as.matrix(Pcov %*% x)
    PY <- as.matrix(Pcov %*% y)
    fit <- LM.fit(x = PX, y = PY, offset = offset, tol = tol)
    fit
  }
  
  z <- if(!is.null(Pcov)) getGLSlm(x, y, Pcov, offset = offset, tol = tol) else 
    if(!is.null(w)) lm.wfit(x = as.matrix(x), y = y, w = w, offset = offset, tol = tol) else
      LM.fit(x = x, y = y, offset = offset, tol = tol)
  
  if(!is.null(Pcov)) {
    z$residuals <- fast.solve(Pcov) %*% z$residuals
    z$fitted.values <- y - z$residuals
  }
  z
}

# NO LONGER USED but retained for potential future use
lm.rrpp.exchange <- function(x, y, Pcov = NULL, w = NULL, offset = NULL, tol = 1e-07){
  
  getGLSlm<- function(x, y, Pcov, offset = offset, method = "qr", tol = tol){
    PX <- Pcov %*% x
    PY <- Pcov %*% y
    fit <- LM.fit(x = PX, y = PY, offset = offset, tol = tol)
    fit
  }
  
  getWlm<- function(x, y, w, offset = offset, method = "qr", tol = tol){
    wts <- sqrt(w)
    fit <- lm.fit(x = as.matrix(x * wts), 
                  y = as.matrix(y * wts), offset = offset, tol = tol)
    fit
  }
  
  z <- if(!is.null(Pcov)) getGLSlm(x = x, y = y, Pcov, offset = offset, tol = tol) else 
    if(!is.null(w)) getWlm(x = x, y = y, w = w, offset = offset, tol = tol) else
      LM.fit(x = x, y = y, offset = offset, tol = tol)
  z
}

# NO LONGER USED but retained for potential future use
package.exchanges <- function(Y, mods, Xs, Terms, model, 
                             offset = NULL, w = NULL,
                             Pcov = NULL, tol = 1e-7, SS.type) {
  Xrs <- Xs$Xrs
  Xfs <- Xs$Xfs
  exchange.args <- list(x = Xrs[[1]], y = Y, Pcov = Pcov, w = w,
                        offset = offset, tol = tol)
  
  reduced = lapply(Xrs, function(x) {
    exchange.args$x <- x
    do.call(lm.rrpp.exchange, exchange.args)
  })
  
  full = lapply(Xfs, function(x) {
    exchange.args$x <- x
    do.call(lm.rrpp.exchange, exchange.args)
  })
  
  model.sets <- list(terms.r = mods$terms.r, terms.f = mods$terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, weights = w,
       Terms = Terms, model = model, Pcov = Pcov, SS.type = SS.type, 
       model.sets = model.sets)
  
}

# NO LONGER USED but retained for potential future use
package.fits <- function(Y, mods, Xs, Terms, model, 
                         offset = NULL, w = NULL,
                         Pcov = NULL, tol = 1e-7, SS.type) {
  Xrs <- Xs$Xrs
  Xfs <- Xs$Xfs
  fit.args <- list(x = Xrs[[1]], y = Y, Pcov = Pcov, w = w,
                   offset = offset, tol = tol)
  
  reduced = lapply(Xrs, function(x) {
    fit.args$x <- x
    do.call(lm.rrpp.fit, fit.args)
  })
  
  full = lapply(Xfs, function(x) {
    fit.args$x <- x
    do.call(lm.rrpp.fit, fit.args)
  })
  
  model.sets <- list(terms.r = mods$terms.r, terms.f = mods$terms.f,
                     Xrs = Xrs, Xfs = Xfs)
  list(reduced = reduced, full = full, offset = offset, weights = w,
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

# getHb
# function to find "hat" matrix for coefficients
# used in lm.rrpp/SS.iter/beta.iter
getHb <- function(Q) {
  Qnames <- Q$dimnames[[2]]
  k <- Q$rank
  R <- Q$R
  U <- Q$Q
  Rs <- try(fast.solve(R), silent = TRUE)
  if(inherits(Rs, "try-error")){
    Rs <- 1
  } 

  res <- as.matrix(tcrossprod(Rs, U))
  
  if(length(Rs) > 1 && !is.null(rownames(res)) &&
     !identical(rownames(Rs), Qnames))
    res <- res[Qnames, ]
  
  if(is.null(rownames(res))) rownames(res) <- Qnames
 
  res
  
}


# checkers
# algorithms to facilitate RRPP iteration stats calculations
# used in lm.rrpp/SS.iter/beta.iter
checkers <- function(Y, Qs, Xs, turbo = FALSE, 
                     Terms, Pcov = NULL, w = NULL) {
  k <- length(attr(Terms, "term.labels"))
  Qr <- Qs$reduced
  Qf <- Qs$full
  n <- NROW(Y)
  Ur <- lapply(Qr, function(q) q$Q)
  Uf <- lapply(Qf, function(q) q$Q)
  kk <- length(Uf)
  
  if(k > 0 && k != kk) k <- kk

  if(!turbo) {
    Hbf <- lapply(Qf, getHb)
    Hbfs <- lapply(Hbf, function(x) Matrix(round(x, 12), sparse = TRUE))
    Hbr <- lapply(Qr, getHb)
    Hbrs <- lapply(Hbr, function(x) Matrix(round(x, 12), sparse = TRUE))  
    
    for(i in 1:max(1,k)) {
      if(object.size(Hbfs[[i]]) < object.size(Hbf[[i]]))
        Hbf[[i]] <- Hbfs[[i]]
      if(object.size(Hbrs[[i]]) < object.size(Hbr[[i]]))
        Hbr[[i]] <- Hbrs[[i]]
    }
    
  } else Hbr <- Hbf <- NULL
  
  Ufull <- Uf[[max(1, k)]]
  
  int <- attr(Terms, "intercept")
  intercept <- rep(int, n)
  Qint <- if(!is.null(Pcov))
    QRforX(Pcov %*% intercept, reduce = FALSE) else if(!is.null(w))
      QRforX(intercept * sqrt(w), reduce = FALSE) else
        QRforX(intercept, reduce = FALSE)

  Hbnull <- tcrossprod(fast.solve(Qint$R), Qint$Q) 
  
  out <- list(Y = Y, Ur = Ur, Uf = Uf, Unull = Qint$Q, Ufull = Ufull,
              Hbr = Hbr, Hbf = Hbf, Hbnull = Hbnull, QR = Qs, k = k,
              realized.trms = names(Xs$Xfs))
  
  out
  
}

# SS.iter
# three functions: main, and two for whether PP is used
# gets appropriate SS vectors for random permutations in lm.rrpp
# generates ANOVA stats

SS.iter <- function(checkrs, ind,  ind_s = NULL,
                    subTest = FALSE, STerm = NULL,
                    print.progress = TRUE, 
                    Parallel.args) {
  
  forking <- Parallel.args$forking
  usecluster <- Parallel.args$usecluster
  cluster <- Parallel.args$cluster
  no_cores <- Parallel.args$ParCores
  
    SS.iter.main(checkrs = checkrs, ind = ind,
                 ind_s = ind_s, subTest = subTest, 
                 STerm = STerm,
                 print.progress = print.progress,
                 no_cores = no_cores, usecluster = usecluster,
                 cluster = cluster, forking = forking)
    
  
}

SS.iter.main <- function(checkrs, ind, ind_s, subTest, STerm,
                         print.progress = TRUE, 
                         no_cores, cluster = NULL, forking = FALSE,
                         usecluster = FALSE) {
  
  Ur <- checkrs$Ur
  Uf <- checkrs$Uf
  Unull <- checkrs$Unull
  Ufull <- checkrs$Ufull
  FR <- checkrs$FR
  Y <- checkrs$Y
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  yh0 <- as.matrix(fastFit(Unull, Y, n, p))
  r0 <- as.matrix(Y - yh0)
  
  k <- checkrs$k
  trms <- checkrs$realized.trms
  rm(checkrs)
  
  perms <- length(ind)
  
  rrpp.args <- list(FR = FR, ind.i = NULL,
                    ind_s.i = NULL, subTest = subTest,
                    STerm = STerm)
  
  rrpp <- function(FR, ind.i, ind_s.i, subTest, STerm) {
    result <- lapply(FR, function(x) x$fitted + x$residuals[ind.i, ])
    if(subTest)
      result[[STerm]] <- FR[[STerm]]$fitted + 
        FR[[STerm]]$residuals[ind_s.i,]
    result
  }
  
  
  ss <- function(ur, uf, y) c(sum(crossprod(ur, y)^2), sum(crossprod(uf, y)^2), 
                              sum(y^2) - sum(crossprod(Ufull, y)^2))
  
  pbbar <- FALSE
  if(print.progress && no_cores > 1){
    cat("\nProgress bar not available for Sums of Squares calculations...\n")
  } else   if(print.progress && no_cores == 1){
    cat(paste("\nSums of Squares calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms, initial = 0, style=3)
    pbbar <- TRUE
  }
  
  if(forking && no_cores > 1) {
    result <- mclapply(1:perms, mc.cores = no_cores, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      if(subTest) rrpp.args$ind_s.i <- ind_s[[j]]
      Yi <- do.call(rrpp, rrpp.args)
      y <- as.matrix(yh0 + r0[x,])
      yy <- sum(y^2)
      if(k > 0) {
        res <- vapply(1:max(1, k), function(j){
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
    
  } else if(usecluster && no_cores > 1) {
  
    result <- parLapply(cluster, 1:perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      if(subTest) rrpp.args$ind_s.i <- ind_s[[j]]
      Yi <- do.call(rrpp, rrpp.args)
      y <- as.matrix(yh0 + r0[x,])
      yy <- sum(y^2)
      if(k > 0) {
        res <- vapply(1:max(1, k), function(j){
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
    stopCluster(cluster)
    
  } else {
    result <- lapply(1:perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      if(subTest) rrpp.args$ind_s.i <- ind_s[[j]]
      Yi <- do.call(rrpp, rrpp.args)
      y <- as.matrix(yh0 + r0[x,])
      yy <- sum(y^2)
      
      if(k > 0) {
        res <- vapply(1:max(1, k), function(j){
          ss(Ur[[j]], Uf[[j]], Yi[[j]])
        }, numeric(3))
        
        SSr <- res[1, ]
        SSf <- res[2, ]
        RSS <- res[3, ]
        
        TSS <- yy - sum(crossprod(Unull, y)^2)
        RSS.model <- yy - sum(crossprod(Ufull, y)^2)
        SS = SSf - SSr
      } else {
        SSr <- SSf <- SS <- RSS <- TSS <- NA
        RSS.model <- TSS <- yy - sum(crossprod(Ufull, y)^2)
      }
      list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
    })
  }
  
  SS <- matrix(sapply(result, "[[", "SS"), max(1, k), perms)
  RSS <- matrix(sapply(result, "[[", "RSS"), max(1, k), perms)
  TSS <- matrix(sapply(result, "[[", "TSS"), max(1, k), perms, 
                byrow = TRUE)
  RSS.model <- matrix(sapply(result, "[[", "RSS.model"), 
                      max(1, k), perms, byrow = TRUE)
  
  res.names <- list(if(k > 0) trms else "Intercept", 
                    c("obs", paste("iter", 1:(perms-1), sep=".")))
  dimnames(SS) <- dimnames(RSS) <- dimnames(TSS)  <- dimnames(RSS.model) <-
    res.names
  
  if(all(is.na(SS))) RSS <- SS <- NULL
  
  if(pbbar) close(pb)
  
  list(SS = SS, RSS = RSS, TSS = TSS, RSS.model = RSS.model)
}


# anova_parts
# construct an ANOVA table from random SS output
# used in lm.rrpp

getRank <- function(Q) {
  if(inherits(Q, "QR")){
    r <- Q$rank
  } else {
    r <- NCOL(removeRedundant(Q))
  }
  return(r)
}

anova_parts <- function(checkrs, SS, full.resid = FALSE){
  SS.type <- checkrs$SS.type
  perms <- NCOL(SS)
  trms <- checkrs$terms
  k <- checkrs$k
  dims <- dim(checkrs$Y)
  n <- dims[1]
  p <- dims[2]
  QRf <- checkrs$QR$full
  QRr <- checkrs$QR$reduced
  
  if(k > 0) {
    df <- unlist(Map(function(qf, qr) 
      getRank(qf) - getRank(qr), 
      QRf, QRr))
    dfe <- n - getRank(QRf[[k]])
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
    
    if(full.resid) {
      SS <- abs(SS - cbind(0, matrix(SS[,1], nrow(SS), ncol(SS) - 1)))
      MS <- abs(MS - cbind(0, matrix(MS[,1], nrow(MS), ncol(MS) - 1)))
      Rsq <- abs(Rsq - cbind(0, matrix(Rsq[,1], nrow(Rsq), ncol(Rsq) - 1)))
      Fs <- MS/RMS
      rownames(Fs) <- rownames(SS)
      cohenf <- NULL
    }
    
    
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
              Fs = Fs, cohenf = cohenf, full.resid = full.resid,
              n = n, p = p, df=df
  )
  out
}

# SS.mean
# support function for calculating SS quickly in SS.iter
# used in lm.rrpp
# NO LONGER USED but retained for potential future use

SS.mean <- function(x, n) if(is.vector(x)) sum(x)^2/n else sum(colSums(x)^2)/n


# getXfromNewData
# make a new X matrix from new data, for a lm.rrpp fit
# used in predict.lm.rrpp
getXfromNewData <- function(fit, newdata){
  Terms <- fit$LM$Terms
  o_vars <- all.vars(Terms)
  Terms <- delete.response(Terms)
  mf_vars <- all.vars(Terms)
  mf <- fit$LM$data
  mf <- mf[o_vars %in% mf_vars]
  nd_vars <- names(newdata)
  factors <- attr(Terms, "factors")
  X <- fit$LM$X
  if(is.null(attr(X, "assign")))
     X <- model.matrix(Terms, mf)
  asn <- attr(X, "assign")
  m <- colMeans(X)
  
  mfn <- as.data.frame(
    matrix(NA, NROW(newdata), NCOL(mf)))
  dimnames(mfn) <- list(rownames(newdata), 
                        mf_vars)
  mfn[, mf_vars %in% nd_vars] <- newdata[, nd_vars]
  
  NA_id <- which(is.na(mfn[1,]))
  if(length(NA_id) > 0) {
    for(i in 1:length(NA_id)) mfn[, NA_id[i]] <- mf[1, NA_id[i]]
  }
  
  N <- NROW(mf)
  mfnn <- rbind(mf, model.frame(Terms, mfn))
  mfnn[is.na(mfnn)] <- 0
  attr(mfnn, "terms") <- Terms
  
  fix <- which(!mf_vars %in% nd_vars)
  term_fix <- which(colSums(as.matrix(factors[fix,])) > 0)
  asn_fix <- which(asn %in% term_fix)
  
  b_ind <- rep(1, length(asn))
  b_ind[asn_fix] <- 0
  m_ind <- 1 - b_ind
  
  nXb <- model.matrix(Terms, mfnn)[-(1:N),]
  mfnn <- mfnn[-(1:N),]
  b_mat <- matrix(b_ind, NROW(mfnn), length(asn), byrow = TRUE)
  m_mat <- matrix(m_ind, NROW(mfnn), length(asn), byrow = TRUE)
  nXm <- matrix(m, NROW(mfnn), length(asn), byrow = TRUE)
  nX <- nXb * b_mat + nXm * m_mat
  rownames(nX) <- rownames(newdata)
  if(is.null(colnames(nX))) colnames(nX) <- colnames(X)
  nX
}


# beta.boot.iter
# gets appropriate beta vectors for coefficients via bootstrap
# used in predict.lm.rrpp

beta.boot.iter <- function(fit, ind) {

  gls <- fit$LM$gls
  QR <- getModels(fit, "qr")
  Qf <- QR$full[[length(QR$full)]]
  rm(QR)
  id <- Qf$dimnames[[2]]
  
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
  if(is.null(Pcov) && !is.null(fit$LM$Cov))
    Pcov <- Cov.proj(fit$LM$Cov)
  
  rrpp.args <- list(fitted = as.matrix(fitted), 
                    residuals = as.matrix(res),
                    ind.i = NULL)
  
  rrpp <- function(fitted, residuals, ind.i) as.matrix(fitted + residuals[ind.i,])
  
  Hf <- tcrossprod(fast.solve(Qf$R), Qf$Q)
  Hfs <- Matrix(round(Hf, 15), sparse = TRUE)
  if(object.size(Hfs) < object.size(Hf)) Hf <- Hfs
  b <- if(gls) fit$LM$gls.coefficients else
    fit$LM$coefficients
  Hf <- Hf[rownames(b),]
  
  betas <- lapply(1:perms, function(j){
    x <-ind[[j]]
    rrpp.args$ind.i <- x
    Yi <- do.call(rrpp, rrpp.args)
    if(!is.null(Pcov)) Yi <- Pcov %*% Yi
    res <- as.matrix(Hf %*% Yi)
    rownames(res) <- id
    res
  })
  betas
  
}

# beta.iter
# three functions: main, and two for whether PP is used
# gets appropriate beta vectors for random permutations in lm.rrpp
# generates distances as statistics for summary
beta.iter <- function(checkrs, ind, ind_s = NULL,
                      subTest = FALSE, STerm = NULL,
                      print.progress = TRUE, 
                      Parallel.args) {
  
  forking <- Parallel.args$forking
  cluster <- Parallel.args$cluster
  usecluster <- Parallel.args$usecluster
  no_cores <- Parallel.args$ParCores
    
    beta.iter.main(checkrs = checkrs, ind = ind, ind_s = ind_s, 
                   subTest = subTest, STerm = STerm,
                   print.progress = print.progress,
                   no_cores = no_cores, usecluster = usecluster,
                   cluster = cluster, forking = forking)
}


beta.iter.main <- function(checkrs, ind, ind_s, subTest,
                           STerm, print.progress = TRUE, 
                           no_cores, cluster = NULL, forking = FALSE,
                           usecluster = FALSE) {
  
  k <- checkrs$k
  trms <- checkrs$realized.trms
    
  Hr <- checkrs$Hbr 
  Hf <- checkrs$Hbf
  pert.rows <- lapply(1:max(1, k), function(j){
    br.nms <- try(rownames(Hr[[j]]), silent = TRUE)
    if(inherits(br.nms, "try-error"))
      br.nms <- "(Intercept)"
    bf.nms <- try(rownames(Hf[[j]]), silent = TRUE)
    if(inherits(bf.nms, "try-error"))
      bf.nms <- "(Intercept)"
    if(length(bf.nms) > length(br.nms)) {
      b <- bf.nms
      a <- br.nms
    } else {
      a <- bf.nms
      b <- br.nms
    }
    b.keep <- setdiff(bf.nms, br.nms)
    res <- which(b %in% b.keep)
    list(res = res, b.keep = b.keep)
  })
  
  b.names <- lapply(pert.rows, function(x) x$b.keep)
  pert.rows <- lapply(pert.rows, function(x) x$res)
  names(b.names) <- names(pert.rows) <- trms
  
  FR <- checkrs$FR
  o <- checkrs$offset
  offst <- !is.null(o)
  perms <- length(ind)
  Y <- checkrs$Y
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  
  rrpp.args <- list(FR = FR, ind.i = NULL, 
                    ind_s.i = NULL, subTest = subTest,
                    STerm = STerm, offst = offst, o = o)
  
  rrpp <- function(FR, ind.i, ind_s.i, subTest, STerm, offst, o) {
    result <- if(offst)
      lapply(FR, function(x) x$fitted + x$residuals[ind.i, ] - o) else
        lapply(FR, function(x) x$fitted + x$residuals[ind.i, ]) 
    if(subTest) {
      result[[STerm]] <- FR[[STerm]]$fitted + 
        FR[[STerm]]$residuals[ind_s.i,]
      if(offst) result[[STerm]] <- result[[STerm]] - o
    }
    result
  }
  
  pbbar <- FALSE
  if(print.progress && no_cores > 1){
    cat("\nProgress bar not available for coefficients estimation...\n")
  } else   if(print.progress && no_cores == 1){
    cat(paste("\nCoefficients estimation:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms, initial = 0, style=3)
    pbbar <- TRUE
  }
  
  getBetas <- function(Hf, Yi, pert.rows){
    Bf <- Map(function(hf, y, p) {
      B <- as.matrix(hf %*% y)
      d <- rowSums(B^2)[p]
      res = list(B = B, d = sqrt(d))
      res
    }, Hf, Yi, pert.rows) 
  } 
  
  if(forking && no_cores > 1) {
    betas <- mclapply(1:perms, mc.cores = no_cores, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      if(subTest) rrpp.args$ind_s.i <- ind_s[[j]]
      Yi <- do.call(rrpp, rrpp.args)
      getBetas(Hf, Yi = Yi, pert.rows)
      
    })
    
  } else if(usecluster && no_cores > 1) {
    
    betas <- parLapply(cluster, 1:perms, function(j){
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      if(subTest) rrpp.args$ind_s.i <- ind_s[[j]]
      Yi <- do.call(rrpp, rrpp.args)
      getBetas(Hf, Yi = Yi, pert.rows)
      
    })
    
  } else {
    betas <- lapply(1:perms, function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      x <-ind[[j]]
      rrpp.args$ind.i <- x
      if(subTest) rrpp.args$ind_s.i <- ind_s[[j]]
      Yi <- do.call(rrpp, rrpp.args)
      getBetas(Hf, Yi = Yi, pert.rows)
      
    })
  }
  
  if(pbbar)  close(pb)

  cnms <- colnames(Y)
  
  betas.out <- lapply(1:max(1, k), function(j){
    res <- lapply(1:perms, function(jj){
      y <- betas[[jj]][[j]]$B
      colnames(y) <- cnms
      y
    })
    names(res) <- names(ind)
    res
  })
  
  names(betas.out) <-trms
  
  d.out <- lapply(1:max(1, k), function(j){
    res <- sapply(1:perms, function(jj){
      y <- betas[[jj]][[j]]$d
      y
    })
    
    res
  })
  
  if(is.list(d.out)) d.out <- do.call(rbind, d.out)
  try(dimnames(d.out) <- list(unlist(b.names), 
                              names(ind)), 
      silent = TRUE)

  if(is.matrix(d.out) && any(rownames(d.out) == "(Intercept)")){
    rmove <- which(rownames(d.out) == "(Intercept)")
    dnames <- rownames(d.out)[-rmove]
    d.out <- d.out[-rmove, ]
    if(is.vector(d.out)){
      d.out <- matrix(d.out, 1, length(d.out))
      colnames(d.out) <- names(ind)
      rownames(d.out) <- dnames
    }
  }

  list(random.coef = betas.out, 
       random.coef.distances = d.out)
}

# ellipse.points
# A helper function for plotting ellipses from non-parametric CI data
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
  if(is.null(x$Fs)) {
    AN <- getANOVAStats(object, stat = "all")
    x[c("SS", "MS", "Rsq", "Fs", "cohenf")] <-
       AN[ c("SS", "MS", "Rsq", "Fs", "cohenf")]
    AN <- NULL
  }
  
  if(is.null(x$Fs)) x$Fs <- x$F
  df <- x$df
  k <- length(df)-2
  kk <- length(object$LM$term.labels)
  if(k > 0 && k != kk) {
    warning(
      paste(
        "\nThis is not an error!  It is a friendly warning.\n",
        "\nANOVA is missing some terms, likely because",
        "\nsome independent variables were redundant.",
        "\nIf the residual SS is 0, results should not be trusted\n", 
        "\nUse suppressWarnings() to turn off these warnings. \n\n", sep = " "),
      noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
  } 
    
  SS <- x$SS
  MS <- x$MS 
  RSS <-x$RSS
  TSS <- x$TSS
  perms <- object$PermInfo$perms
  pm <- object$PermInfo$perm.method
  if(pm == "RRPP" && object$PermInfo$full.resid)
    pm <- "FMRP"
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
        
        warning(
          paste(
            "\nThis is not an error!  It is a friendly warning.\n",
            "\nCalculating effect size on SS is illogical with GLS.",
            "\nEffect type has been changed to F distributions.\n",
            "\nUse suppressWarnings() to turn off these warnings. \n\n", sep = " "),
          noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
        
        effect.type = "F"
      }
      
      if(effect.type == "MS") {
        warning(
          paste(
            "\nThis is not an error!  It is a friendly warning.\n",
            "\nCalculating effect size on MS is illogical with GLS.",
            "\nEffect type has been changed to F distributions.\n",
            "\nUse suppressWarnings() to turn off these warnings. \n\n", sep = " "),
          noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
        
        effect.type = "F"
      }
    } else est <- "OLS"
    
    ow <- options()$warn
    options(warn = -1)
    
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
      if(!inherits(Z, "matrix")) Z <- matrix(Z, 1, length(Z))
      P.val <- apply(Z, 1, pval) 
      Z <- apply(Z, 1, effect.size)
      } else P.val <- NULL
    
    Residuals <- c(df[k+1], RSS[[1]], RSS[[1]]/df[k+1], 
                   RSS[[1]]/TSS[[1]], rep(NA, 3))
    Total <- c(df[k+2], TSS[[1]], rep(NA, 5))
    tab <- data.frame(Df = df[1:k], SS=SS, MS = MS, Rsq = Rsq, 
                      F = Fs, Z = Z, P.val = P.val)
    tab <- rbind(tab, Residuals = Residuals, Total = Total)
    colnames(tab)[NCOL(tab)] <- paste("Pr(>", effect.type, ")", sep="")
    class(tab) = c("anova", class(tab))
    SS.type <- x$SS.type
    
    if(any(tab$Df == 0)){
      rmove <- which(tab$Df == 0)
      tab <- tab[-rmove,]
    }
    
    options(warn = ow)
    
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
  PermInfo <- getPermInfo(refModel, attribute = "all")
  ind <- PermInfo$perm.schedule
  perms <- length(ind)
  
  if(refModel$LM$gls) {
    X <- refModel$LM$X
    Y <- refModel$LM$Y
    Pcov <- refModel$LM$Pcov
    Cov <- refModel$LM$Cov  
    if(!is.null(Cov) && is.null(Pcov)) {
      Pcov <- Cov.proj(Cov)
    }
    X <- if(!is.null(Pcov)) Pcov %*% X else
      X * sqrt(refModel$LM$weights)
    Y <- if(!is.null(Pcov)) Pcov %*% Y else
      Y * sqrt(refModel$LM$weights)
    
  } else {
    X <- refModel$LM$X
    Y <- refModel$LM$Y
  }
  
  B <- if(refModel$LM$gls) refModel$LM$gls.coefficients else 
    refModel$LM$coefficients
  
  refQR <- getModels(refModel, "qr")
  refQR <- refQR$full[[length(refQR$full)]]
  
  K <- length(lm.list)
  QRlist <- lapply(1:K, function(j){
    m <- lm.list[[j]]
    mQR <- getModels(m, "qr")
    mQR <- mQR$full[[length(mQR$full)]]
    mQR
  })
  
  U <- as.matrix(refQR$Q)
  n <- refModel$LM$n
  p <- refModel$LM$p
  Yh <- as.matrix(fastFit(U, Y, n, p))
  R <- as.matrix(Y) - Yh
  
  
  Ulist <- lapply(1:K, function(j){
    QRlist[[j]]$Q
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
    Pcov <- refModel$LM$Pcov
    Cov <- refModel$LM$Cov  
    if(!is.null(Cov) && is.null(Pcov)) {
      Pcov <- Cov.proj(Cov)
    }
    int <- if(!is.null(Pcov))  Pcov %*% rep(int, n) else
      sqrt(refModel$LM$weights)
  } else int <- rep(int, n)
  
  Qint <- QRforX(int, reduce = FALSE)
  U0 <- Qint$Q
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
  
  dfe <- n - c(getRank(refQR), unlist(lapply(1:K, 
                     function(j) getRank(QRlist[[j]]))))
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

aov.me <- function(object){
  
  perms <- getPermInfo(object, "perms")

  out <- list(table = object$AOV, perm.method = "RRPP", perm.number = perms,
              est.method = "OLS", SS.type = "Within-subject II", effect.type = "SNR",
              call = object$call)
  class(out) <- "anova.lm.rrpp"
  out
}

# getSlopes
# gets the slopes for groups from a lm.rrpp fit
# used in pairwise
getSlopes <- function(fit, x, g){
  
  beta <- fit$LM$random.coef
  k <- length(beta)
  beta <- beta[[k]]
  n <- fit$LM$n
  p <- fit$LM$p
  X <- fit$LM$X
  X <- X[, colnames(X)  %in% rownames(beta[[1]])]
  getFitted <- function(b) X %*% b
  fitted <- lapply(beta, getFitted)
  Xn <- model.matrix(~ g * x + 0)
  Q <- QRforX(Xn)
  H <- tcrossprod(fast.solve(Q$R), Q$Q)
  getCoef <- function(f) H %*% f
  Coef <- lapply(fitted, getCoef)
  group.slopes <- function(B){ # B of form ~ group * x + 0
    gp <- QRforX(model.matrix(~g))$rank
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
# after constraining covariates to mean values
# used in pairwise
getLSmeans <- function(fit, g){
  beta <- fit$LM$random.coef
  k <- length(beta)
  beta <- beta[[k]]
  n <- fit$LM$n
  dat <- fit$LM$data
  covCheck <- sapply(dat, class)
  for(i in 1:length(covCheck)) if(covCheck[i] == "numeric") 
    dat[[i]] <- mean(dat[[i]])
  L <- model.matrix(fit$LM$Terms, data = dat)
  L <- L[, colnames(L)  %in% rownames(beta[[1]])]
  getFitted <- function(b) {
    b <- b[rownames(b) %in% colnames(L),]
    L %*% b
  } 
  fitted <- lapply(beta, getFitted)
  Xn <- model.matrix(~ g + 0)
  Q <- QRforX(Xn, reduce = FALSE)
  H <- tcrossprod(fast.solve(Q$R), Q$Q)
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
  
  y <- sapply(M, function(x) as.vector(as.dist(x)))
  if(is.vector(y)) y <- matrix(y, 1, length(y))
  P <- apply(y, 1, pval)
  D <- as.dist(M[[1]])
  D[1:length(D)] <- P
  P <- as.matrix(D)
  diag(P) <- 1
  P
}

# effect.list
# effect size across a list
# used in pairwise
effect.list <- function(M){
  
  y <- sapply(M, function(x) as.vector(as.dist(x)))
  if(is.vector(y)) y <- matrix(y, 1, length(y))
  Z <- apply(y, 1, effect.size)
  D <- as.dist(M[[1]])
  D[1:length(D)] <- Z
  as.matrix(D)
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

# r.summary.from.list
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
  QR <- QRforX(x)
  H <- tcrossprod(fast.solve(QR$R), QR$Q)
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
  gls <- fit$LM$gls
  n <- fit$LM$n
  p <- fit$LM$p.prime
  X <- as.matrix(fit$LM$X)
  Y <- as.matrix(fit$LM$Y)
  R <- if(gls) fit$LM$gls.residuals else fit$LM$residuals
  PCA <- ordinate(R, tol = tol, rank. = min(c(pc.no, p)))
  rnk <- length(PCA$d)
  w <- fit$LM$weights
  Cov <- fit$LM$Cov
  Pcov <- NULL
  if(gls) {
    Pcov <- try(getModelCov(fit, "Pcov"), 
                silent = TRUE)
    if(inherits(Pcov, "try-error")) Pcov <- NULL
  }
  
  R <- if(gls) {
    if(!is.null(Pcov)) Pcov %*% R else R * sqrt(w)
  } else R
  
  R <- as.matrix(R)
  
  if(!is.null(w)) {
    excl <- w <= 0
    logdetC <- sum(log(w[!excl]))
  } else {
    logdetC <- if(gls) determinant(Cov, logarithm = TRUE)$modulus else 0
  }
  
  if(NCOL(R) > rnk) R <- ordinate(R, rank. = rnk)$x
  Sig <- as.matrix(crossprod(R) / n)
  if(kappa(Sig) > 1e10) Sig <- RiReg(Sig, R)
  logdetSig <- determinant(Sig, logarithm = TRUE)$modulus
  
  ll <- -0.5 * (n * rnk * log(2 * pi) + rnk * logdetC +
                  n * logdetSig + n * rnk) 
  
  list(logL = ll, rank = rnk)
  
}

cov.trace <- function(fit) {
  n <- fit$LM$n
  p <- fit$LM$p.prime
  if(fit$LM$gls){
    Pcov <- try(getModelCov(fit, "Pcov"), 
               silent = TRUE)
    if(inherits(Pcov, "try-error"))
      Pcov <- NULL
    
    Sig <- if(!is.null(Pcov))  crossprod(Pcov %*% fit$LM$gls.residuals)/n else
      crossprod(fit$LM$gls.residuals * sqrt(fit$LM$weights))/n
    
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
      nfit <- lm.fit(as.matrix(Xg * sqrt(w)), as.matrix(Yg * sqrt(w)))
      fitted <- Yg * sqrt(w) - res * sqrt(w)
      nY <- (fitted + nfit$residuals)/sqrt(w)
    } else if(!is.null(Pcov)) {
      PY <- Pcov %*% Yg
      nfit <- lm.fit(as.matrix(Pcov %*% Xg), PY)
      fitted <- PY - as.matrix(Pcov %*% res)
      nY <- fast.solve(Pcov) %*% (fitted + nfit$residuals)
    } else {
      nfit <- lm.fit(as.matrix(Xg), as.matrix(Yg))
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
  if(!is.null(Pcov)) V <- crossprod(Pcov %*% Yc)/n
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

# looPCOne
# finds one PC vector via cross-validation
# used in looCV

looPCOne <-function(fit, n.ind) {
  X <- fit$LM$X
  Y <- fit$LM$Y
  x <- X[n.ind,]
  y <- Y[n.ind,]
  X <- X[-n.ind,]
  Y <- Y[-n.ind,]
  gls <- fit$LM$gls
  Pcov <- if(gls) fit$LM$Pcov[-n.ind, -n.ind] else NULL
  QR <- if(gls) QRforX(Pcov %*% X) else QRforX(X)
  H <- tcrossprod(fast.solve(QR$R), QR$Q)
  B <- if(gls)  H %*% (Pcov %*% Y) else
    H %*% Y
  S <- svd(X %*% B)
  y %*% S$v
}

# looPCOne
# finds each PC vector via cross-validation for a model fit
# used in looCV
looPCAll<-function(fit, ...) {
  n <- fit$LM$n
  gls <- fit$LM$gls
  
  ord.args <- list(...)
  ord.args$Y <- fit$LM$Y
  QR <- getModels(fit, "qr")
  QR <- QR$full[[length(QR$full)]]
  ord.args$A <- tcrossprod(QR$Q)
  if(is.null(ord.args$tol)) ord.args$tol <- 1e-6
  
  ord <- do.call(ordinate, ord.args)
  k <- 1:length(ord$d)
  
  res <- t(sapply(1:n, function(j) looPCOne(fit, j)))
  dimnames(res) <- if(gls) dimnames(fit$LM$gls.fitted) else
    dimnames(fit$LM$fitted)
  d <- svd(fit$LM$fitted)$d
  res <- res[,k]
  
  res <- ordinate(res, ord$x, rank. = max(k))
  list(raw = ord, cv = res)
}
