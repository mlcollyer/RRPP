#'  Quantify and compare shape change trajectories 
#'
#'  Function estimates attributes of multivariate trajectories 
#'
#'  The function quantifies multivariate trajectories from a set 
#'  of observations, and assesses variation 
#'  in attributes of the trajectories via RRPP. A trajectory is defined 
#'  by a sequence 
#'  of points in the data space. These trajectories can be quantified for 
#'  various attributes (their size, orientation, 
#'  and shape), and comparisons of these attribute enable the statistical 
#'  comparison of shape change 
#'  trajectories (Collyer and Adams 2007; Adams and Collyer 2007; 
#'  Adams and Collyer 2009; Turner et al. 2010; Collyer and Adams 2013). 
#'  
#'  This function is a modified version of \code{\link{pairwise}}, retaining 
#'  the least squares (LS) means as trajectory points.
#'  Analysis starts with a \code{\link{lm.rrpp}} fit (but a procD.lm fit from 
#'  geomorph can also be used).  LS means are calculated using a grouping
#'  variable.  Data can be trajectories, as a start(sensu Adams and Cerney 2007), 
#'  or trajectories can be calculated from data using a factorial model 
#'  (in which case
#'  trajectory points are defined by factor levels).  
#'
#'  This function produces statistics that can be summarized with the 
#'  \code{\link{summary.trajectory.analysis}} function.  The summaries
#'  are consistent with those in the \code{\link{summary.pairwise}} 
#'  function, pertaining to trajectory attributes including,
#'  magnitude difference (MD), the difference in path lengths of trajectories; 
#'  trajectory correlations (TC), better
#'  thought of as angular differences between trajectory principal axes; and if 
#'  trajectories have three or more points,
#'  shape difference (SD), the square root of summed squared point differences, 
#'  after scaling, centering, and rotating trajectories.  The SD is
#'  the "Procrustes" distance between trajectories (Adams and Collyer 2009), 
#'  much the same way as the shape difference between anatomical landmark
#'  configurations in geometric morphometrics.  If attribute = "TC" is chosen 
#'  for the summary, then the angle type ("rad" or "deg",
#'  can be chosen for either radians and degrees, respectively, to return 
#'  angles between principal axes.)
#'  
#'  Plotting can be performed with \code{\link{plot.trajectory.analysis}} and 
#'  \code{\link{add.trajectories}}.  The former
#'  plots all principal component scores for the data, and allows point-by-point 
#'  control of plot parameters.  The later
#'  adds trajectories points and lines, with constrained control.  By saving the 
#'  plot.trajectory.analysis
#'  object, plotting information can be retained and advanced plotting can be 
#'  performed.  See examples below.
#'  
#' @param fit A linear model fit using \code{\link{lm.rrpp}}.
#' @param fit.null An alternative linear model fit to use as a null model for 
#' RRPP, if the null model
#' of the fit is not desired.  Note, if RRPP = FALSE (FRPP rather than RRPP), 
#' then the null model has only an intercept.
#' If the null model is uncertain, using \code{\link{reveal.model.designs}} 
#' will help elucidate the inherent null model used.
#' @param groups A factor or vector coercible to factor that defines trajectories.
#' @param traj.pts Either a single value or a vector coercible to factor to 
#' define trajectory points.  If only a single value, 
#' it is assumed that the data are already in the form, 
#' y1p1, y2p1, y3p1, ...., y2p2, y2p2, y3p2, ..., yjp1, yjp2, yjp3, ..., yjpk, 
#' for j variables comprising k trajectory points;
#' i.e., traj.pts = k.  If a factor, then a group * traj.pt factorial model 
#' is assumed, where traj.pts defines the levels for points within groups.
#' @param pca A logical value to optionally project group:point means onto 
#' principal components (perform PCA on a covariance matrix of the means)
#' This option only applies to factorial designs (traj.pts is a factor).
#' @param print.progress A logical value to indicate whether a progress bar 
#' should be printed to the screen.  
#' This is helpful for long-running analyses.

#' @export
#' @keywords analysis
#' @author Dean Adams and Michael Collyer
#' @return An object of class "trajectory.analysis" returns a list of the 
#' following:
#'   \item{LS.means}{LS.means from pairwise function.}
#'   \item{trajectories}{Trajectories from every permutation.}
#'   \item{PD}{Path distances of trajectories from every permutation.}
#'   \item{MD}{Magnitude differences between trajectories from every permutation.}
#'   \item{TC}{Trajectory correlations from every permutation.}
#'   \item{SD}{Trajectory shape differences from every permutation.}
#' @references Adams, D. C., and M. M. Cerney. 2007. Quantifying biomechanical 
#' motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. The analysis of 
#' character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework 
#' for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. Analysis of two-state 
#' multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @references Collyer, M. L., and D. C. Adams. 2013. Phenotypic trajectory 
#' analysis: comparison of shape change patterns 
#' in evolution and ecology. Hystrix 24: 75-83.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for 
#' analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' 
#' @examples 
#' ### Analysis of sexual dimorphism vectors (factorial approach)
#' data(Pupfish)
#' fit <- lm.rrpp(coords ~ Pop * Sex, data = Pupfish, iter = 199)
#' reveal.model.designs(fit)
#' TA <- trajectory.analysis(fit, groups = Pupfish$Pop, 
#' traj.pts = Pupfish$Sex, print.progress = FALSE)
#' 
#' # Magnitude difference (absolute difference between path distances)
#' summary(TA, attribute = "MD") 
#' 
#' # Correlations (angles) between trajectories
#' summary(TA, attribute = "TC", angle.type = "deg") 
#' 
#' # No shape differences between vectors
#' summary(TA, attribute = "SD") 
#' 
#' # Retain results
#' TA.summary <- summary(TA, attribute = "MD")
#' TA.summary$summary.table
#' 
#' # Plot results
#' TP <- plot(TA, pch = as.numeric(Pupfish$Pop) + 20, bg = as.numeric(Pupfish$Sex),
#' cex = 0.7, col = "gray")
#' add.trajectories(TP, traj.pch = c(21, 22), start.bg = 1, end.bg = 2)
#' legend("topright", levels(Pupfish$Pop), pch =  c(21, 22), pt.bg = 1)
#' 
#' ### Analysis when data are already trajectories (motion paths)
#' 
#' # data are planar Cartesian coordinates (x, y) across 5 points (10 variables)
#' data(motionpaths)
#' fit <- lm.rrpp(trajectories ~ groups, data = motionpaths, iter = 199)
#' TA <- trajectory.analysis(fit, groups = motionpaths$groups, traj.pts = 5)
#' 
#' # Magnitude difference (absolute difference between path distances)
#' summary(TA, attribute = "MD") 
#' 
#' # Correlations (angles) between trajectories
#' summary(TA, attribute = "TC", angle.type = "deg") 
#' 
#' # Shape differences between trajectories 
#' summary(TA, attribute = "SD") 
#' 
#' TP <- plot(TA, pch = 21, bg = as.numeric(motionpaths$groups),
#' cex = 0.7, col = "gray")
#' add.trajectories(TP, traj.pch = 21, traj.bg = 1:4)

trajectory.analysis <- function(fit, fit.null = NULL, groups, 
                                    traj.pts, pca = TRUE, print.progress = FALSE){
  
  n <- fit$LM$n
  perms <- fit$PermInfo$perms
  
  if(is.numeric(traj.pts) && length(traj.pts) == 1) {
    g2 <- NULL
    tp <- traj.pts
  }
  
  if(is.vector(traj.pts) || is.factor(traj.pts)) {
    if(length(traj.pts) == 1) {
      g2 <- NULL
      tp <- traj.pts
    } else {
      g2 <- traj.pts
      tp <- NULL
      p <- levels(g2)
    }
  }
  g1 <- as.factor(groups)
  
  if(!is.null(g2) && NCOL(g2) > 1) 
    stop("traj.pts can be either a single value or a factor, not a matrix.\n",
                                        call. = FALSE)
  if(NCOL(g1) > 1) stop("Groups must be a single factor.\n",
                        call. = FALSE)
  if(length(g1) != n)
    stop("The number of observations does not match group length.\n",
                          call. = FALSE)
  
  if(is.null(tp)) groups <- interaction(g1, g2, lex.order = TRUE) else 
    groups <- g1
  groups <- factor(groups)
  gp.rep <- by(groups, groups, length)
  if(!all(gp.rep > 1)) 
    stop("Not every trajectory point has replication (more than one observation).\n",
                            call. = FALSE)
  
  PW <- pairwise(fit, fit.null, groups, covariate = NULL, 
                     print.progress = print.progress)
  
  means <- PW$LS.means
  if(is.null(tp) && pca) {
    PCA <- if(fit$LM$gls) prcomp(fit$LM$gls.fitted) else prcomp(fit$LM$fitted)
    rot <- PCA$rotation
    Y.cent <- matrix(colMeans(fit$LM$Y), NROW(means[[1]]), 
                     ncol(means[[1]]), byrow = TRUE)
    means <- lapply(means, function(x) (x - Y.cent) %*% rot)
  } else PCA <- NULL
  
  
  if(!is.null(tp)){
    p <- ncol(means[[1]])/tp
    if(p != floor(p)) 
      stop("The number of variables divided by the number of trajectory points is not an integer")
    gl <- levels(g1)
    
    trajectories <- lapply(1:perms, function(j){
      m <- means[[j]]
      tj <- lapply(1:length(gl), function(jj){
        mm <- m[jj,]
        matrix(mm, traj.pts, p, byrow = TRUE)
      })
      names(tj) <- gl
      tj
    })
    
    line.trajectories <- means
  }
  
  if(is.null(tp)) {
    gl <- levels(g1)
    traj.list <- list()
    nt <- nrow(means[[1]])
    tpts <-  nt / length(gl)
    start <- seq(1, nt - tpts + 1, tpts)
    end <- seq(tpts, nt, tpts)
    for(i in 1:length(gl)) traj.list[[i]] <- seq(start[i], end[i], 1)
    names(traj.list) <- gl
    
    trajectories <- lapply(1:perms, function(j){
      m <- means[[j]]
      tj <- lapply(1:length(gl), function(jj){
        m[traj.list[[jj]],]
      })
      names(tj) <- names(traj.list)
      tj
    })
    
    line.trajectories <- lapply(1:perms, function(j) {
      m <- trajectories[[j]]
      res <- sapply(1:length(gl), function(jj) {
        mm <- m[[jj]]
        as.vector(t(mm))
      })
      res <- t(res)
      rownames(res) <- gl
      res
    })
  }
  
  # Pairwise differences in path length
  
  PD <- lapply(1:perms, function(j){
    tj <- trajectories[[j]]
    ts <- trajsize(tj)
    names(ts) <- gl
    ts
  })
  
  MD <- lapply(PD, function(x) as.matrix(dist(x)))
  
  # Pairwise correlations
  
  tn <- length(gl)
  Tcor <- if(fit$LM$p > 1) lapply(1:perms, function(j){
    x <- trajectories[[j]]
    to <- trajorient(x, tn)
    dimnames(to) <- list(gl, gl)
    to
  }) else NULL

  # Pairwise shape differences
  
  p <- NROW(trajectories[[1]][[1]])
  
  if(p > 2) {
    
    if(print.progress) {
      cat("Beginning trajectory shape analysis.\nThis could take long for highly multivariate data.\n")
      pb <- txtProgressBar(min = 0, max = perms, initial = 0, style=3)
    }
    
    SD <- lapply(1:perms, function(j) {
      x <- trajectories[[j]]
      ts <- trajshape(x)
      dimnames(ts) <- list(gl, gl)
      if(print.progress) setTxtProgressBar(pb,j)
      ts
    })
    
    if(print.progress) close(pb)
    
  } else SD <- NULL
  
  
  names(trajectories) <- names(MD) <- 
    names(PD) <- c("obs", paste("iter", 1:(perms - 1), sep = "."))
  if(!is.null(Tcor)) names(Tcor) <- names(PD)
  if(!is.null(SD)) names(SD) <- names(PD)
  
  # output
  out <- list(LS.means = means, trajectories = trajectories, PD = PD, 
              MD = MD, TC = Tcor, SD = SD, pca = PCA, 
              fit = fit, n.trajectories = if(is.matrix(trajectories[[1]])) 1 else length(trajectories[[1]]),
              n.points = p,
              type = if(is.null(g2)) "single.factor" else "factorial")
  
  class(out) <- "trajectory.analysis"
  out
  
}

