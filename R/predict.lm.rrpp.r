#' predict for lm.rrpp model fits
#'
#' @description Computes predicted values from an \code{\link{lm.rrpp}} model fit, using bootstrapped residuals
#' to generate confidence intervals.  (Residuals are the residuals of the lm.rppp fit, not its null model.  The bootstrap
#' procedure resamples residual vectors with replacement.)
#' The bootstrap permutations use the same number of iterations and seed as used
#' in the \code{\link{lm.rrpp}} model fit. A \code{\link{predict.lm.rrpp}} object can be plotted using various options.
#' See \code{\link{plot.predict.lm.rrpp}}.
#'
#' @param object Object from \code{\link{lm.rrpp}}.
#' @param newdata Data frame of either class \code{\link{data.frame}} or \code{\link{rrpp.data.frame}}.  If null,
#' the data frame from the lm.rrpp fit will be used, effectively calculating all fitted values and
#' their confidence intervals.
#' @param confidence The desired confidence interval level for prediction.
#' @param ... Other arguments (currently none)
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' # See examples for lm.rrpp to see how predict.lm.rrpp works in conjunction
#' # with other functions
#' 
#' data(Pupfish)
#' names(Pupfish)
#' Pupfish$logSize <- log(Pupfish$CS) # better to not have functions in formulas
#'
#' fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", data = Pupfish) 
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
predict.lm.rrpp <- function(object, newdata, confidence = 0.95, ...) {
  if(!inherits(object, "lm.rrpp")) stop("Object is not class lm.rrpp")
  Terms <- object$LM$Terms
  if (missing(newdata) || is.null(newdata)) {
    newdata <- model.frame(Terms, data = object$LM$data)
    full.predict <- TRUE
  } else full.predict = FALSE
  if(!inherits(newdata, "data.frame") && !inherits(newdata, "rrpp.data.frame"))
    stop("newdata must be an object of class data.frame or rrpp.data.frame")
  if(confidence < 0 || confidence > 1) stop("Confidence level must be between 0 and 1")
  Y <- object$LM$Y
  if(full.predict){
    nX <- object$LM$X
    n <- NROW(nX)
  } else {
    fl <- object$LM$term.labels
    tl <- intersect(fl, names(object$LM$data))
    pl <- match(tl, names(newdata))
    if(all(is.na(pl))) stop("No variables in newdata match variables in lm.rrpp fit")
    gl <- names(newdata)[na.omit(pl)]
    if(length(gl) != length(tl)) {
      cat("\nWarning: Not all variables in model accounted for in newdata.")
      cat("\nMissing variables will be averaged from observed data for prediction.\n\n")
    }
    nd <- as.data.frame(newdata[match(gl, names(newdata))])
    if(any(is.na(pl))){
      add.l <- tl[is.na(match(tl, names(newdata)))]
      nda <- sapply(1:length(add.l), function(j){
        x <- rep(0, NROW(nd))
        x
      })
      colnames(nda) <- add.l
      nd <- data.frame(cbind(nd, nda))
      nd <- nd[, match(tl, names(nd))]
    }
    od <- object$LM$data
    yl <- setdiff(names(od), tl)
    df.add <- names(od[yl])
    nd.names <- names(nd)
    for(i in 1:length(df.add)) nd <- cbind(nd, 0)
    names(nd) <- c(nd.names, df.add)
    nX <- model.matrix(Terms, data = nd) 
    oX <- object$LM$X
    n <- NROW(nX)
    nX <- nX[,intersect(colnames(nX), colnames(oX))]
    cnm <- match(colnames(oX), colnames(nX))
    if(any(is.na(cnm))){
      df.add <- setdiff(colnames(oX), colnames(nX))
      nX.names <- colnames(nX)
      for(i in 1:length(df.add)) nX <- cbind(nX, 0)
      colnames(nX) <- c(nX.names, df.add)
      isall0 <- function(x) all(x == 0)
      nX.check <- apply(nX, 2, isall0)
      for(i in 1:length(nX.check)) 
        if(isTRUE(nX.check[[i]])) nX[,i] <- mean(oX[i,])
    }
  }
  PI <- object$PermInfo$perm.schedule
  seed <- attr(PI, "seed")
  perms <- length(PI)
  indb <- boot.index(length(PI[[1]]), perms -1, seed)
  k <- length(object$LM$term.labels)
  res <- object$LM$wResiduals
  fitted <- object$LM$wFitted
  X <- object$LM$X * sqrt(object$LM$weights)
  Q <- object$LM$wQR
  H <- tcrossprod(solve(qr.R(Q)), qr.Q(Q))
  if(object$LM$gls) {
    P <- object$LM$Pcov
    PY <- crossprod(P, object$LM$Y * sqrt(object$LM$weights))
    PX <- crossprod(P, X)
    glsFit <- lm.fit(PX, PY)
    fitted <- as.matrix(glsFit$fitted.values)
    res <- as.matrix(glsFit$residuals)
    Q <- qr(PX)
    H <- tcrossprod(solve(qr.R(Q)), qr.Q(Q))
  }
  n <- NROW(nX)
  beta.args <- list(f = fitted, r = res, h = H, ind.i = NULL)
  coefs <- lapply(1:length(indb), function(j){
    x <- indb[[j]]
    beta.args$ind.i <- x
    do.call(beta.boot, beta.args)
  })
  predM <- function(b) as.matrix(nX %*% b)
  preds <- lapply(coefs, predM)
  alpha = 1 - confidence
  lcl <- function(x) quantile(x, prob = alpha/2)
  ucl <- function(x) quantile(x, prob = 1 - alpha/2)
  meanV <- Reduce("+", preds)/length(preds)
  preds.array <- simplify2array(preds)
  LCL <- t(apply(preds.array, c(1,2), lcl))
  UCL <- t(apply(preds.array, c(1,2), ucl))
  if(!is.null(dim(LCL))) {
    if(dim(LCL)[2] == n) LCL <- t(LCL)
  }
  if(!is.null(dim(UCL))) {
    if(dim(UCL)[2] == n) UCL <- t(UCL)
  }
  
  if(dim(Y)[1] == 1){
    pc.preds <- preds
    pcLCL <- LCL
    pcUCL <- UCL
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
  
  data.name = as.character(object$call[[2]])[[2]]
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
