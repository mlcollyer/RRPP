#' Reveal the inter-subject variability from a measurement error analysis
#'
#' Function produces both a list of inter-subject Euclidean distance matrices, 
#' based on replicate measurements of the same subjects, and one matrix that 
#' summarizes the variability among the inter-subject distances, across subjects.  
#' This function can be considered a tool for the evaluation of subject 
#' estimate precision.  The function, \code{\link{plot.interSubVar}} can produce a 
#' heat map of the inter-subject variability.
#'
#' 
#' @param ME A measurement error object
#' @param type A value to indicate the type of variability (statistic)
#' to measure, which can be
#' one of range (the maximum value minus the minimum value, not the two values),
#' standard deviation (sd), variance (var), or coefficient of variation (cv).  No 
#' attempt is made to assure the distribution of values is appropriate for the 
#' statistics.  For example, if only two replicates are available, using sd, var, or
#' cv might not be wise.  Or if the replicated values are exact, cv 
#' will be NA (and other stats will be 0).  Choice of statistic should consider
#' the distribution of values.
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @return An object of class \code{interSubVar} is a list containing the 
#' following
#' \item{var.map}{A distance matrix object with values that map the variability
#' statistic used for inter-subject Euclidean distances.}
#' \item{distance.mats}{The inter-subject distance matrices for every replicate.}
#' \item{subject.order}{A vector of subject levels in the order that was used to
#' guarantee consistent sorting across distance matrices.}
#' \item{var.map}{The variability type (statistic) that was used.}
#' @examples
#' # Measurement error analysis on simulated data of fish shapes
#' 
#' data(fishy)
#' 
#' # Analysis unconcerned with groups 
#' 
#' ME1 <- measurement.error(
#'   Y = "coords",
#'   subjects = "subj",
#'   replicates = "reps",
#'   data = fishy)
#' 
#' anova(ME1)
#' ICCstats(ME1, subjects = "Subjects", with_in = "Systematic ME")
#' plot(ME1)
#' 
#' # Analysis concerned with groups 
#' 
#' ME2 <- measurement.error(
#'   Y = "coords",
#'   subjects = "subj",
#'   replicates = "reps",
#'   groups = "groups",
#'   data = fishy)
#' 
#' anova(ME2)
#' ICCstats(ME2, subjects = "Subjects", 
#'   with_in = "Systematic ME", groups = "groups")
#' P <- plot(ME2)
#' focusMEonSubjects(P, subjects = 18:20, shadow = TRUE)
#'  
interSubVar <- function(ME, type = c("range", "sd", "var", "cv")){
  type <- match.arg(type)
  if(!type %in% c("range", "sd", "var", "cv")) type <- "range"
  
  statfx <- if(type == "cv") function(x) {sd(x) / mean(x) * 100} else
    if(type == "sd") function(x) sd(x) else if(type == "var") function(x)
      var(x) else function(x) {max(x) - min(x)}
  
  Y <- ME$LM$data$Y
  subj <- ME$LM$data$subjects
  reps <- ME$LM$data$replicates
  rep.levels <- levels(reps)
  subj.levels <- levels(subj)
  
  tb <- table(subj, reps)
  if(any(tb == 0))
    stop(paste("The ME design is not balanced and not all inter-subject distances",
               "can be estimated for every replicate.\n", sep = " "), call. = FALSE)
  
  result <- lapply(1:length(rep.levels), function(j){
    keep <- which(reps == rep.levels[j])
    Yj <- Y[keep, ]
    subj.j <- subj[keep]
    rownames(Yj) <- subj.j
    Yj <- Yj[subj.levels, ]
    dist(Yj)
  })
  
  names(result) <- rep.levels
  
  result2 <- sapply(result, as.vector)
  result3 <- apply(result2, 1, statfx)
  
  d <- result[[1]]
  d[1:length(d)] <- result3
  
  out <- list(var.map = d,
              distance.mats = result,
              subject.order = attr(d, "Labels"),
              type = type)
  
  class(out) <- "interSubVar"
  out
  
}