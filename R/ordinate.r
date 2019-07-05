
ordinate <- function(Y, A = NULL, Cov = NULL, scale. = FALSE, 
                     tol = NULL, rank. = NULL, newdata = NULL) {
  Y <- try(as.matrix(Y), silent = TRUE)
  if(inherits(Y, "try-error"))
    stop("Y must be a matrix or data frame.\n", call. = FALSE)

  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  I <- diag(n)
  if(is.null(A)) A <- I
  if(!is.matrix(A))
    stop("A must be a matrix with the same number of rows as Y\n", 
         call. = FALSE)
  if(NROW(A) != n)
    stop("A must be a matrix with the same number of rows as Y\n", 
         call. = FALSE)
  if(is.null(rownames(A))) rownames(A) <- rownames(Y) else {
    nnames <- length(intersect(rownames(Y), rownames(A)))
    if(nnames > n)
      stop("The row names of A are not the same as the row names of Y\n",
           call. = FALSE)
    if(isSymmetric(A)) A <- A[rownames(Y), rownames(Y)] else
      A <- A[rownames(Y),]
  }
  X <- matrix(1, n)
  rownames(X) <- rownames(Y)
  if(!is.null(Cov)) Pcov <- Cov.proj(Cov, rownames(Y))
  if(is.null(Cov)) H <- tcrossprod(X)/n else 
    H <- X %*% solve(crossprod(crossprod(Pcov, X))) %*% crossprod(crossprod(Pcov, X), Pcov)
  
  Z <- scale(as.matrix((I - H) %*% Y), center = FALSE, scale = scale.)
  sc <- attr(Z, "scaled:scale")
  if (any(sc == 0)) 
    stop("cannot rescale a constant/zero column to unit variance")
  
  k <- if (!is.null(rank.)) {
    stopifnot(length(rank.) == 1, is.finite(rank.), as.integer(rank.) > 
                0)
    min(as.integer(rank.), n, p)
  } else min(n, p)
  
  s <-  if(!is.null(Cov)) svd(crossprod(A, Pcov %*% Z), 
                              nu = 0, nv = k) else svd(crossprod(A, Z), 
                                                       nu = 0, nv = k)
   
  j <- seq_len(k)
  s$v <- s$v[,j]
  x <- Z %*% s$v
  s$d <- apply(x, 2, sd)
  if (!is.null(tol)) {
    rank <- sum(s$d > (s$d[1L] * tol))
    if (rank < k) {
      j <- seq_len(k <- rank)
      s$v <- s$v[, j, drop = FALSE]
      s$d <- s$d[j]
      x <- x[j]
    }
  }
  dimnames(s$v) <- list(colnames(Z), paste0("PC", j))
  
  r <- list(d = s$d^2, sdev = s$d, rot = s$v, 
            center = colMeans(H %*% Y), 
            scale = if (is.null(sc)) FALSE else sc,
            x = x)

  if(!is.null(newdata)){
    xn <- as.matrix(newdata)
    xn <- scale(xn, center = r$center, scale = r$scale)
    if(NCOL(Z) != NCOL(xn)) stop("Different number of variables in newdata\n", call. = FALSE) 
    r$xn <- xn %*% s$v
  }

  class(r) <- "ordinate"
  r
  
}