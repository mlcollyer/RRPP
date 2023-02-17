#' predict for lm.rrpp model fits
#'
#' @description Computes predicted values from an \code{\link{lm.rrpp}} 
#' model fit, using bootstrapped residuals
#' to generate confidence intervals.  (Residuals are the residuals of 
#' the lm.rppp fit, not its null model.  The bootstrap
#' procedure resamples residual vectors with replacement.)
#' The bootstrap permutations use the same number of iterations and 
#' seed as used
#' in the \code{\link{lm.rrpp}} model fit. A \code{\link{predict.lm.rrpp}} 
#' object can be plotted using various options.
#' See \code{\link{plot.predict.lm.rrpp}}.
#' 
#' Note that if data offsets are used (if the offset argument is used 
#' when fitting a \code{\link{lm.rrpp}} model),
#' they are ignored for estimating coefficients over iterations.  
#' Offsets are subtracted from data in \code{\link[stats]{lm}} and 
#' added to predicted values in \code{\link[stats]{predict.lm}}, 
#' effectively adjusted the intercept and then un-adjusting
#' it for predictions.  This causes problems if the newdata have a 
#' different number of observations than the original
#' model fit.
#' 
#'
#' @param object Object from \code{\link{lm.rrpp}}.
#' @param newdata Data frame of either class \code{\link{data.frame}} 
#' or \code{\link{rrpp.data.frame}}.  If null,
#' the data frame from the lm.rrpp fit will be used, effectively calculating 
#' all fitted values and
#' their confidence intervals.  If a numeric variable is missing from newdata, 
#' an attempt to average the values
#' will be made in prediction; i.e., least squares means for factor levels 
#' can be found.  All factors used in the
#' \code{\link{lm.rrpp}} fit should be represented in the newdata data frame, 
#' with appropriate factor levels.
#' @param confidence The desired confidence interval level for prediction.
#' @param block An optional factor for blocks within which to restrict resampling
#' permutations.
#' @param ... Other arguments (currently none)
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' \dontrun{
#' # See examples for lm.rrpp to see how predict.lm.rrpp works in conjunction
#' # with other functions
#' 
#' data(Pupfish)
#' 
#' # CS is centroid (fish) size
#' fit <- lm.rrpp(coords ~ log(CS)  + Sex*Pop, 
#' SS.type = "I", data = Pupfish, iter = 999) 
#'
#' # Predictions (holding alternative effects constant)
#' 
#' shapeDF <- expand.grid(Sex = levels(Pupfish$Sex), Pop = levels(Pupfish$Pop))
#' rownames(shapeDF) <- paste(shapeDF$Sex, shapeDF$Pop, sep = ".")
#' shapeDF
#' 
#' shapePreds <- predict(fit, shapeDF)
#' summary(shapePreds)
#' summary(shapePreds, PC = TRUE)
#' 
#' shapePreds99 <- predict(fit, shapeDF, confidence = 0.99)
#' summary(shapePreds99, PC = TRUE)
#' 
#' # Plot prediction
#' 
#' plot(shapePreds, PC = TRUE)
#' plot(shapePreds, PC = TRUE, ellipse = TRUE)
#' plot(shapePreds99, PC = TRUE)
#' plot(shapePreds99, PC = TRUE, ellipse = TRUE)
#' }
#' 
predict.lm.rrpp <- function(object, newdata = NULL, block = NULL,
                            confidence = 0.95, ...) {
  if(!inherits(object, "lm.rrpp")) stop("Object is not class lm.rrpp")
  Terms <- object$LM$Terms
  trms <- object$LM$term.labels
  k <- length(trms)
  TT <- delete.response(Terms)
  
  if(!missing(newdata) || !is.null(newdata)) {
    
    if(!inherits(newdata, "data.frame") && !inherits(newdata, "rrpp.data.frame"))
      stop("newdata must be an object of class data.frame or rrpp.data.frame")
    if(confidence < 0 || confidence > 1) 
      stop("Confidence level must be between 0 and 1")
    
    if(inherits(newdata, "rrpp.data.frame")) {
      vars <- all.vars(TT)
      refd <- names(newdata) %in% vars
      if(!all(refd))
        stop("\nOne or more variables in newdata are not model terms.\n",
             call. = FALSE)
      if(length(refd) != length(newdata))
        stop("\nIt's not possible to coerce the rrpp.data.frame into a 
             useable data frame for prediction.\n",
             call. = FALSE)
      
      class(newdata) <- "list"
    }
    
  }
  
  if(!is.null(newdata) && is.list(newdata)) newdata <- as.data.frame(newdata)
  
  if(missing(newdata) || is.null(newdata)) {
    nX <- object$LM$X
  } else {
    nX <- matrix(colMeans(object$LM$X), NROW(newdata), NCOL(object$LM$X), 
                 byrow = TRUE)
    colnames(nX) <- colnames(object$LM$X)
    rownames(nX) <- rownames(newdata)
    
    o.names <- all.vars(TT)
    n.names <- names(newdata)
    tm <- match(o.names, n.names)
    if(all(is.na(tm)))
      stop("\nVariables in newdata do not match variables used in lm.rrpp fit",
           call. = FALSE)
    if(any(is.na(tm))) {
      
      warning(
        paste(
          "\nThis is not an error!  It is a friendly warning.\n",
          "\nNot all variables in model accounted for in newdata.",
          "\nMissing variables will be averaged from observed data for prediction.\n",
          "\nUse options(warn = -1) to turn off these warnings. \n\n", sep = " "),
        noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
      
    }
    
    nform <- formula(TT[which(!is.na(tm))])
    vars <- colnames(nX)[colnames(nX) %in% colnames(mX)]
    if(length(vars) == 0)
      stop("\nVariables in newdata do not match variables used in lm.rrpp fit",
           call. = FALSE)
    nX[, vars] <- mX[, vars]
  }
  
  o <- object$LM$offset
  if(!is.null(o)) offst = TRUE else offst = FALSE
  
  n <- NROW(nX)
  p <- NCOL(object$LM$Y)
  
  PI <- object$PermInfo$perm.schedule
  seed <- attr(PI, "seed")
  perms <- length(PI)
  indb <- boot.index(length(PI[[1]]), perms -1, block, seed)
  k <- length(object$Models$full)
  betas <- beta.boot.iter(object, indb)
  
  predM <- function(b) as.matrix(nX[, rownames(b)] %*% b)
  preds <- lapply(betas, predM)
  
  alpha = 1 - confidence
  lcl <- function(x) quantile(x, prob = alpha/2)
  ucl <- function(x) quantile(x, prob = 1 - alpha/2)
  meanV <- Reduce("+", preds)/length(preds)
  
  getCL <- function(preds, CL) {
    x <- unlist(preds)
    ind <- rep(1:length(meanV), perms)
    as.vector(by(x, ind, CL))
  }
  
  LCL <- getCL(preds, lcl)
  UCL <- getCL(preds, ucl)
  
  LCL <- matrix(LCL, NROW(meanV), NCOL(meanV))
  UCL <- matrix(UCL, NROW(meanV), NCOL(meanV))
  
  dimnames(LCL) <- dimnames(UCL) <- dimnames(meanV)
  
  if(!is.null(dim(LCL))) {
    if(dim(LCL)[2] == n) LCL <- t(LCL)
  }
  if(!is.null(dim(UCL))) {
    if(dim(UCL)[2] == n) UCL <- t(UCL)
  }
  
  if(p == 1){
    pc.preds <- preds
    pcLCL <- LCL
    pcUCL <- UCL
    pc.meanV <- Reduce("+", pc.preds)/length(pc.preds)
    pca <- NULL
  } else {
    pca <- prcomp(preds[[1]])
    d <- length(which(zapsmall(pca$sdev) > 0))
    R <- pca$rotation[,1:d]
    rotate <- function(x) center(x) %*% R
    pc.preds <- lapply(preds, rotate)
    pc.meanV <- Reduce("+", pc.preds)/length(pc.preds)
    pcLCL <- as.matrix(sapply(1:n, function(j){
      x <- sapply(1:length(pc.preds), function(jj){
        pc.preds[[jj]][j,]
      })
      if(!is.matrix(x)) z <- lcl(x) else z <- apply(x, 1, lcl)
      z
    }))
    
    pcUCL <- as.matrix(sapply(1:n, function(j){
      x <- sapply(1:length(pc.preds), function(jj){
        pc.preds[[jj]][j,]
      })
      if(!is.matrix(x)) z <- ucl(x) else z <- apply(x, 1, ucl)
      z
    }))
    
    if(!is.null(dim(pcLCL))) {
      if(dim(pcLCL)[2] == n) pcLCL <- t(pcLCL)
    }
    if(!is.null(dim(pcUCL))) {
      if(dim(pcUCL)[2] == n) pcUCL <- t(pcUCL)
    }
  }
  
  data.name = as.character(object$call[[2]][[2]])
  rn <- rownames(nX)
  if(is.null(rownames(nX))) rn <- as.character(1:NROW(nX))
  colnames(meanV) <- colnames(LCL) <- colnames(UCL) <- colnames(object$LM$Y)
  rownames(meanV) <- rownames(pc.meanV) <- rownames(LCL) <- rownames(UCL) <- 
    rownames(pcLCL) <- rownames(pcUCL) <- rn
  out <- list(mean = meanV, lcl = LCL, ucl = UCL, 
              pc.mean = pc.meanV, pc.lcl = pcLCL, pc.ucl = pcUCL,
              confidence = confidence, 
              data.name = data.name, random.predicted = preds,
              random.predicted.pc = pc.preds, pca=pca)
  class(out) <- "predict.lm.rrpp"
  out
}
