#' Pairwise comparisons of model effects 
#'
#' Function generates pairwise statistics for comparing model fits and
#' returns important statistics for hypothesis tests.
#'
#' The function statistically compares the effect sizes of two or more models fit
#' and evaluated using RRPP.  Input for the function is a list of fitted models
#' of the class 'model.comparison', whose options included type = 'Z' and 
#' verbose = TRUE when the models were compared with that function.
#' 
#' A two-sample test is performed on each pair of models, comparing the strength 
#' of model fits to one another (Collyer and Adams 2025). This might be used to 
#' compare the strength of fit of the data to differing statistical models (as in 
#' model selection) or for comparing the fit across differing 
#' datasets for the same model to determine whether the strength of fit in one dataset
#' is greater than that found in another (see Collyer and Adams 2025). In the latter
#' case, one is advised to include a vector containing the sample sizes of each
#' dataset, so that two-sample tests may account for differences in sample size. 
#' 
#' @param ... Either an object of class \code{\link{model.comparison}}, 
#' or several objects of class \code{\link{lm.rrpp}}.  If the former, arguments of
#' type = 'Z' and verbose = TRUE are required.  If the latter, a model comparsion analysis
#' will first be performed with these arguments.
#' @param nsamp An optional vector containing the sample sizes for each model fit
#' @param two.tailed A logical value to indicate whether a two-tailed test (typical and default) should be performed.
#' @param predictor An optional vector to be passed to \code{\link{model.comparison}}, if used.
#' @param tol An optional value to be passed to \code{\link{model.comparison}}, if used.
#' @param pc.no An optional value to be passed to \code{\link{model.comparison}}, if used.
#' @param gls.null An optional logical value to be passed to \code{\link{model.comparison}}, if used.
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' 
#' @return A list containing the following
#' \item{sample.z}{A vector of model effect sizes.}
#' \item{pairwise.z}{A matrix of pairwise test statistics comparing model effect sizes.}
#' \item{pairwise.P}{A matrix of pairwise significance levels.}
#' \item{tails}{Number of tails used for P-value calculation.}
#' @references Collyer, M.L., and D.C. Adams. 2025. Permutational Biometry. 
#'   Volume 1: Univariate Data. Iowa State University Digital Press. (Forthcoming).
#' @examples
#' \dontrun{
#'  data(Pupfish)
#'  Pupfish$logSize <- log(Pupfish$CS)
#'  fit1 <- lm.rrpp(coords ~ logSize, data = Pupfish,
#'  print.progress = FALSE)
#'  fit2 <- lm.rrpp(coords ~ Pop, data = Pupfish,
#'  print.progress = FALSE)
#'  fit3 <- lm.rrpp(coords ~ Sex, data = Pupfish, 
#'  print.progress = FALSE)
#'  fit6 <- lm.rrpp(coords ~ logSize + Sex * Pop, data = Pupfish, 
#'  print.progress = FALSE)
#'  Mod.C <- model.comparison(fit1, fit2, fit3, fit6,
#'  pc.no = 4, type = "Z", verbose = TRUE)
#'  res <- pairwise.model.Z(Mod.C)
#'  summary(res, stats.table = TRUE)
#'  summary(res, stats.table = FALSE)
#'  }

pairwise.model.Z <- function(..., 
                             nsamp = NULL, 
                             two.tailed = TRUE,
                             predictor = NULL,
                             tol = NULL,
                             pc.no = NULL,
                             gls.null = FALSE
                             ){
  dots <- list(...)
  if(length(dots) == 1){
    if(!inherits(dots[[1]], "model.comparison"))
      stop("\nNot a class model.comparison object.\n", call. = FALSE)
    MC <- dots[[1]]
  } else {
    MC <- model.comparison(..., type = "Z", verbose = TRUE,
                           predictor = predictor, pc.no = pc.no,
                           tol = tol, gls.null = gls.null)
  }
  rm(dots)
  if(is.null(MC$random.logL))
    stop("\nThe model.comparsion analysis requires type = 'Z' and verose = TRUE.\n",
         call. = FALSE)
  
  dists <- MC$random.logL
  k <- length(dists)
  list.names <- names(dists)
  if(is.null(list.names)) list.names <- seq(1:k)
  sdn <- function(x) sqrt(sum((x - mean(x))^2) / length(x))
  k.combn <- combn(k, 2, simplify = FALSE)
  bct <- lapply(dists, function(x) box.cox(x)$transformed)
  tails <- if(two.tailed) 2 else 1
  list.drs <- sapply(1:k, function(j) bct[[j]][1] - mean(bct[[j]])) 
  list.sds <- sapply(1:k, function(j) sdn(bct[[j]]))
  list.zs <- sapply(1:k, function(j) effect.size(dists[[j]], center=TRUE))
  z12 <- sapply(1:length(k.combn), function(j){
    x <- k.combn[[j]]
    a <- x[1]
    b <- x[2]
    r1 <- list.drs[a]
    r2 <- list.drs[b] 
    sd1 <- list.sds[a]
    sd2 <- list.sds[b]
    if (is.null(nsamp)) {n1 <- n2 <- 1}
    if (!is.null(nsamp)) {n1 <- nsamp[a]
      n2 <- nsamp[b]}
    (r1 - r2)/sqrt(sd1^2 * (2 * n1 / ((n1 + n2))) + 
                     sd2^2 * (2 * n2 / ((n1 + n2))))
  }) 
    z12.p <- sapply(1:length(z12), function(j) 
      pnorm(abs(z12[[j]]), lower.tail = FALSE) * tails)
  d <- rep(0, k)
  names(d) <- list.names
  D <-dist(d)
  z12.pw <- p12.pw <- D
  for(i in 1:length(z12)) z12.pw[i] <-z12[i]
  for(i in 1:length(z12)) p12.pw[i] <-z12.p[i]
  names(list.zs) <- names(list.sds) <-list.names
  pairwise.z <- as.matrix(z12.pw)
  pairwise.P <- as.matrix(p12.pw)
  diag(pairwise.P) <- 1
  
  out <- list(sample.z = list.zs,
              pairwise.z = abs(pairwise.z),
              pairwise.P = pairwise.P, 
              tails = tails)
  
  attr(out, "class") <- "pairwise.model.Z"
  out
}