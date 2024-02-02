# Important function for performing lm.rrpp or lm.rrpp.ws
# All arguments must be vetted by primary functions.

.lm.rrpp <- function(f1, subjects = NULL, 
                     sub.var.name = NULL,
                     iter, turbo, seed = NULL, 
                     int.first = FALSE,
                     RRPP, full.resid, block = NULL,
                     SS.type, data, Cov,
                     delta, gamma, 
                     print.progress, Parallel, 
                     verbose, ...) {
  
  Parallel.args <- Parallel.setup(Parallel)
  
  STerm <- sub.lev <- NULL
  subTest <- !is.null(subjects)
  L <- c(as.list(environment()), list(...))
  names(L)[which(names(L) == "f1")] <- "formula"
  
  if(int.first) ko = TRUE else ko = FALSE
  SS.type <- L$SS.type
  full.resid <- L$full.resid
  
  if(full.resid && SS.type != "III"){
    SS.type = "III"
    warning(
      paste(
        "\nThis is not an error!  It is a friendly warning.\n",
        "\nBecause a permutation of full model residuals was chosen,",
        "\nSS.type is being forced to be III, as this is the only applicable",
        "\nestimation method when using full model residuals as exchangeable units",
        "\nunder the null hypotheses of model effects.  Additionally, all random",
        "\nANOVA statistics will have the form, |random.stat - observed.stat|.\n",
        "\nUse options(warn = -1) to turn off these warnings. \n\n", sep = " "),
        noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
  }
  
  dots <- list(...)
  if(length(dots) > 0) {
    w <- dots$weights
    o <- dots$offset
  } else w <- o <- NULL
  
  if(!is.null(w) && !is.null(Cov)) {
    w <- NULL
    warning(
    paste(
      "\nThis is not an error!  It is a friendly warning.\n",
      "\nIt is not possible to use both a Cov matrix and weights.",
      "\nBoth are inputs to perform generalized least squares estimation of coefficients,",
      "\nbut at present only one covariance matrix can be used.",
      "\nAs a result, weights are set to NULL.",
      "\nYou could consider adjusting your Cov matrix; e.g., Cov <- Cov * 1/weights,",
      "\nmaking sure the Cov matrix and weights are ordered consistently.\n",
      "\nUse options(warn = -1) to turn off these warnings. \n\n", sep = " "),
    noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
  }
  
  Terms <- D <- NULL
  
  if(print.progress) {
    cat("\nPlease be aware that printing progress slows down the analysis (perhaps slightly).\n")
    cat("\nPreliminary Model Fit...\n")
  } 
  
  if(inherits(f1, "lm")) {
    exchange.args <- f1[attributes(f1)$names %in% 
                          c("terms", "offset", "weights")]
    exchange.args$model <- f1$model
    exchange.args$Y <- Y <- as.matrix(f1$model[[1]])
    exchange.args$tol <- f1$qr$tol
    exchange.args$SS.type <- SS.type
    names(exchange.args)[which(names(exchange.args) == "terms")] <- "Terms"
    if("weights" %in% names(exchange.args))
      names(exchange.args)[which(names(exchange.args) == "weights")] <- "w"
    Terms <- f1$terms
  }
  
  if(inherits(f1, "formula")) {
    exchange.args <- lm.args.from.formula(L)
    if(!is.null(exchange.args$D)) D <- exchange.args$D
    exchange.args <- exchange.args[c("Terms", "Y", "model")]
    exchange.args$tol <- 1e-7
    exchange.args$SS.type <- SS.type
    Terms <- exchange.args$Terms
    Y <- as.matrix(exchange.args$Y)
  }
  
  id <- get.names(Y)
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  if(is.null(id)) {
    id <- 1:n
    Y <- add.names(Y, id)
  }
  
  attr(Terms, ".Environment") <- NULL
  term.labels <- attr(Terms, "term.labels")
  k <- length(term.labels)
  
  offst <- (!is.null(o)) 
  weighted <- (!is.null(w))
  exchange.args <- c(exchange.args, list(w = w, offset = o))
  
  if(!inherits(f1, c("lm", "formula")))
    stop("\nf1 must be either a formula or class lm objects.\n",
         call. = FALSE)
  
  gls <- FALSE
  ols <- TRUE
  
  if(!is.null(w)) {
    if(NROW(w) != n)
      stop("The number of weights does not match the number of observations.  This could be because of missing data.\n",
           call. = FALSE)
    gls <- TRUE
    ols <- FALSE
  }
  
  if(subTest){
    STerm <- which(term.labels == sub.var.name)
    subTest <- !(length(STerm) == 0)
  }

  if(subTest) {
    subjects <- as.factor(subjects)
    sub.lev <- levels(subjects)
  } 
  
  if(length(subTest) == 0) {
    subTest <- FALSE
    STerm <- sub.lev <- NULL
  }
  
  if(!subTest && !is.null(sub.var.name)) {
    STerm <- sub.lev <- NULL
  }
  
  Pcov <- CovEx <- NULL
  
  if(!is.null(Cov) && !subTest) {
    
    Cov.name <- deparse(substitute(Cov))
    Cov.match <- match(Cov.name, names(data))
    if(all(is.na(Cov.match))) Cov <- Cov else Cov <- data[[Cov.match]]
    if(is.null(rownames(Cov))) rownames(Cov) <- 1:n
    if(is.null(colnames(Cov))) colnames(Cov) <- 1:n
    if(!inherits(Cov, "Matrix") && !inherits(Cov, "matrix")) 
      stop("The covariance matrix must be a matrix.")
    if(!is.null(id) && !is.null(rownames(Cov))) {
      if(length(setdiff(id, rownames(Cov))) > 0)
        stop("Data names and covariance matrix names do not match.\n", call. = FALSE)
      Cov <- Cov[id, id]
      
    }
    
    Pcov <- Cov.proj(Cov, id = id)
    gls <- TRUE
    ols <- FALSE
  } 
  
  if(!is.null(Cov) && subTest) {
    
    cov.lev <- levels(as.factor(dimnames(Cov)[[1]]))
    if(!all(sub.lev %in% cov.lev) || !all(cov.lev %in% sub.lev))
      stop("\nThere is a mismatch between subject levels and covariance levels\n
         Make sure there is direct correspondence between subjects in the model\n
         and subjects in the covariance matrix.  Check spelling of names.\n", 
           call. = FALSE)
    Xsub <- model.matrix(~subjects + 0)
    colnames(Xsub) <- sub.lev
    Xs <- Matrix(Xsub, sparse = TRUE)
    nn <- nrow(Cov)
    if(nn == nrow(Xsub)) cov.type <- 1 else 
      if(nn == length(sub.lev)) cov.type <- 2 else
        stop(paste("There is an irreconcilable covariance matrix, in terms of",
                   "subject number or subject factor levels.  Please make sure that",
                   "the number of Cov rows or columns match the number of subject levels",
                   "or the number of observations.", sep = " "), call. = FALSE)
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
      delta[which(zapsmall(delta) == 0)] <- 0.001
      
      n.list <- as.vector(by(subjects, subjects, length))
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
    Pcov <- Cov.proj(CovEx)
    gls <- TRUE
    ols <- FALSE
    
  }
  

  X.args <- exchange.args[c("Terms", "Y", "SS.type", "tol", "model")]
  X.args$subjects.term <- STerm
  Xs <- do.call(getXs, X.args)
  X <- Xs$Xfs[[length(Xs$Xfs)]]
  
  Qs <- lapply(1:2, function(j){
    X.j <- Xs[[j]]
    kk <- length(X.j)
    res <- lapply(1:kk, function(jj){
      x <- X.j[[jj]]
      if(!is.null(Pcov)) x <- Pcov %*% x
      if(!is.null(w)) x <- x * sqrt(w)
      qr(x)
    })
    res
  })
  
  Qs.sparse <- lapply(1:2, function(j){
    X.j <- Xs[[j]]
    kk <- length(X.j)
    res <- lapply(1:kk, function(jj){
      x <- X.j[[jj]]
      if(!is.null(Pcov)) x <- Pcov %*% x
      if(!is.null(w)) x <- x * sqrt(w)
      x.sparse <- Matrix(round(as.matrix(x), 15), sparse = TRUE)
      qr(x.sparse)
    })
    res
  })
  
  names(Qs) <- names(Qs.sparse) <- c("reduced", "full")
  
  ind <- perm.index(n, iter = iter, block = block, seed = seed)
  perms <- iter + 1
  
  if(subTest) {
    ind <- perm.index(n, iter = iter, block = subjects, seed = seed)
    ind_s <- perm.index(n, iter = iter, block = block, seed = seed)
  } else ind_s <- NULL

  checkers.args <- list(Y = Y, Qs = Qs, Qs.sparse = Qs.sparse, Xs = Xs,
                        turbo = turbo, Terms = Terms, Pcov = Pcov, w = w)
  cks <- do.call(checkers, checkers.args)

  Qs <- Qs.sparse <- checkers.args <- NULL
  
  TY <- if(!is.null(Pcov)) Pcov %*% Y else if(!is.null(w)) Y * sqrt(w) else Y
  
  PCA <- p > (n - 1) 
  
  if(PCA) {
    TYp <- ordinate(TY, tol = 1e-10)$x
    p.prime <- ncol(TYp)
  } else {
    TYp <- TY
    p.prime <- p
  }
  
  Ur <- cks$Ur
  kk <- length(Ur)
  
  if(RRPP) {
    if(full.resid) {
      Uf <- cks$Uf
      FR <- lapply(1:max(1, kk), function(j){
        fitted <- as.matrix(fastFit(Uf[[j]], TY, n , p))
        residuals <- as.matrix(TY - fitted)
        out <- list(fitted = fitted, residuals = residuals)
      })
      
      if(subTest) {
        FR[[STerm]]$fitted <- as.matrix(fastFit(Uf[[STerm]], TY, n , p))
        FR[[STerm]]$residuals <- as.matrix(TY - FR[[STerm]]$fitted)
      }
      
      Uf <- NULL
      
    } else {
      
      FR <-lapply(1:max(1, kk), function(j){
        fitted <- as.matrix(fastFit(Ur[[j]], TY, n , p))
        residuals <- as.matrix(TY - fitted)
        out <- list(fitted = fitted, residuals = residuals)
      })
      
      if(subTest) {
        FR[[STerm]]$fitted <- as.matrix(fastFit(Ur[[STerm]], TY, n , p))
        FR[[STerm]]$residuals <- as.matrix(TY - FR[[STerm]]$fitted)
      }
      
    }

  } else {
    FR <- lapply(1:max(1, kk), function(j){
      fitted <-  matrix(0, n, p)
      residuals <- as.matrix(TY)
      list(fitted = fitted, residuals = residuals)
    })
  }
  
  cks$FR <- FR
  
  if(!turbo) {
    cks$offset <- o
    cks$Y <- TY
    cks$trms <- term.labels
    
    beta.args <- list(checkrs = cks, ind = ind,
                      ind_s = ind_s,
                      subTest = subTest,
                      STerm = STerm,
                      Parallel.args = Parallel.args, 
                      print.progress = print.progress)
    
    betas <- do.call(beta.iter, beta.args)
    random.coef <- betas$random.coef
    random.coef.distances <- betas$random.coef.distances
    betas <- beta.args <- NULL
  } else random.coef <- random.coef.distances <- NULL
  
  cks$Y <- TYp
  
  if(PCA){
    
    if(RRPP) {
      if(full.resid) {
        Uf <- cks$Uf
        FR <- lapply(1:max(1, k), function(j){
          fitted <- as.matrix(fastFit(Uf[[j]], TYp, n , p.prime))
          residuals <- as.matrix(TYp - fitted)
          list(fitted = fitted, residuals = residuals)
        })
        
        if(subTest) {
          FR[[STerm]]$fitted <- as.matrix(fastFit(Uf[[STerm]], TYp, n , p))
          FR[[STerm]]$residuals <- as.matrix(TYp - FR[[STerm]]$fitted)
        }
        
        Uf <- NULL
        
      } else {
        FR <- lapply(1:max(1, k), function(j){
          fitted <- as.matrix(fastFit(Ur[[j]], TYp, n , p.prime))
          residuals <- as.matrix(TYp - fitted)
          list(fitted = fitted, residuals = residuals)
        })
        
        if(subTest) {
          FR[[STerm]]$fitted <- as.matrix(fastFit(Ur[[STerm]], TYp, n , p))
          FR[[STerm]]$residuals <- as.matrix(TYp - FR[[STerm]]$fitted)
        }
      }
      
      
      } else {
      FR <- lapply(1:max(1, k), function(j){
        fitted <- matrix(0, n, p.prime)
        residuals <- as.matrix(TYp)
        list(fitted = fitted, residuals = residuals)
      })
      }
    
    cks$FR <- FR
  }
  
  SS.args <- list(checkrs = cks, ind = ind, 
                  ind_s = ind_s, 
                  subTest = subTest, STerm = STerm,
                  print.progress = print.progress,
                  Parallel.args = Parallel.args)
  FR <- NULL
  SS <- do.call(SS.iter, SS.args)
  cks$SS.type <- SS.type

  ANOVA <- anova_parts(cks, SS, full.resid)
  
  if(!verbose){
    ANOVA$MS <- ANOVA$Fs <- ANOVA$cohenf <- NULL
  }
  
  SS.args <- NULL
  
  obs.fit <- lm.rrpp.fit(X, Y, Pcov = Pcov, w = w, offset = o, 
                         tol = exchange.args$tol)
  
  QR <- obs.fit$qr
  U <- qr.Q(QR)
  R <- suppressWarnings(qr.R(QR))
  Hb <- as.matrix(tcrossprod(fast.solve(R), U))
  if(!identical(colnames(X), colnames(R)))
    Hb <- Hb[colnames(X),]
  coefficients <- as.matrix(Hb %*% TY)
  if(is.null(rownames(coefficients)))  
    rownames(coefficients) <- rownames(Hb)
  R <- U <- QR <- NULL

  QR <- if(gls && !is.null(Pcov)) qr(Pcov %*% X) else
    if(gls && !is.null(w)) qr(X * sqrt(w)) else qr(X)
  
  LM <- list(form = formula(Terms), 
             coefficients = coefficients,
             ols = ols,
             gls = gls,
             Y = Y,  
             X = X, 
             n = n, p = p, p.prime = p.prime,
             QR = QR,
             Terms = Terms, term.labels = term.labels,
             data = exchange.args$model,
             random.coef = if(verbose) random.coef else NULL,
             random.coef.distances = if(verbose) 
               random.coef.distances else NULL
  )
  
  attr(LM$Terms,".Environment") <- NULL
  attr(LM$form,".Environment") <- NULL
  if(verbose && !turbo) {
    environment(LM$random.coef) <- 
      environment(LM$random.coef.distances) <- NULL
  }
  
  LM$weights <- w
  LM$offset <- o
  
  if(gls) {
    names(LM)[[2]] <- "gls.coefficients"
    LM$gls.fitted <- as.matrix(obs.fit$fitted.values)
    LM$gls.residuals <- as.matrix(obs.fit$residuals)
    rownames(LM$gls.fitted) <- rownames(LM$gls.residuals) <- id
    LM$gls.mean <- if(NCOL(LM$gls.fitted) > 1) colMeans(LM$gls.fitted) else
      mean(LM$gls.fitted)
  } else {
    LM$fitted <- as.matrix(obs.fit$fitted.values)
    LM$residuals <- as.matrix(obs.fit$residuals)
    rownames(LM$fitted) <- rownames(LM$residuals) <- id
    LM$mean <- if(NCOL(LM$fitted) > 1) colMeans(LM$fitted) else
      mean(LM$fitted)
  }
  
  if(!is.null(Cov)) {
    LM$Cov <- if(!is.null(CovEx)) CovEx else Cov
    LM$Pcov <- if(verbose) Pcov else NULL
    rm(Cov, Pcov, CovEx)
  }
  
  
  PermInfo <- list(perms = perms,
                   perm.method = ifelse(RRPP==TRUE,"RRPP", "FRPP"), 
                   full.resid = full.resid, block = block,
                   perm.schedule = ind, perm.seed = seed)
  if(!verbose) PermInfo$perm.schedule <- NULL
  rm(ind)
  rm(ind_s)
  
  out <- list(call = match.call(), 
              LM = LM, ANOVA = ANOVA, PermInfo = PermInfo, turbo = turbo)
  environment(out$PermInfo$perm.seed) <- 
    environment(out$PermInfo$perm.sschedule) <- NULL
  
  if(k == 0 && print.progress)
    cat("\nNo terms for ANOVA; only RSS calculated in each permutation\n")
  
  if(!is.null(exchange.args$D)) {
    qrf <- LM$QR
    D.coef <- qr.coef(qrf, D)
    out$LM$dist.coefficients <- D.coef
  }
  
  if(verbose) {
    
    Models <-lapply(1:2, function(j){
      res <- lapply(1:max(1, kk), function(jj){
        X <- Xs[[j]][[jj]]
        qr <- cks$QR[[j]][[jj]]
        list(X = X, qr = qr)
      })
      names(res) <- cks$realized.trms
      res
    })
    names(Models) <- c("reduced", "full")
    
    Model.Terms <- .getTerms(fit = NULL, Terms = Terms, SS.type = SS.type)
    
    for(i in 1:2){
      for(j in 1:max(1, kk)){
        Models[[i]][[j]]$terms <- Model.Terms[[i]][[j]]
      }
    }
    
  } else Models <- NULL

  out$Models <- Models
  out$verbose <- verbose
  
  out$subjects <- subjects
  out$subjects.var <- sub.var.name
  out$subTest <- subTest
  
  class(out) = "lm.rrpp"
  
  out
  
}