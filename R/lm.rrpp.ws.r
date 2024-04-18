#' Linear Model Evaluation with RRPP performed within subjects
#'
#' Function performs a linear model fit over many random permutations of 
#' data, using a randomized residual permutation procedure restricted to
#' subjects.
#'
#' The function fits a linear model using ordinary least squares (OLS) or 
#' generalized least squares (GLS) estimation of coefficients over any 
#' number of random permutations of the data, but the permutations are mostly 
#' restricted to occur with subject blocks for any model terms other than subjects.  
#' All functionality should resemble that of \code{\link{lm.rrpp}}.  However, 
#' an argument for research subjects is also required.  The purpose of this function 
#' is to account for the non-independence among observations of research subjects 
#' (due to sampling within subjects), while also allowing for the non-independence 
#' among subjects to be considered (Adams and Collyer, submitted).  
#' 
#' By comparison, the covariance matrix option in \code{\link{lm.rrpp}} must have a
#' one-to-one match to observations, which can be matched by the row names of the data.
#' In this function, the covariance matrix can be the same one used in \code{\link{lm.rrpp}}
#' but the number of observations can be greater.  For example, if subjects are
#' species or some other level of taxonomic organization, data can comprise measurements
#' on individuals.  Users have the option to expand the covariance matrix for subjects
#' or input one they have generated.
#' 
#' Irrespective of covariance matrix type, the row names of the data matrix must match the
#' subjects.  This step assures that the analysis can proceed in \code{\link{lm.rrpp}}.  It
#' is also best to make sure to use an \code{\link{rrpp.data.frame}}, so that the subjects
#' can be a name in that data frame.  For example, if research subjects are species and
#' data (observations) are collected from individuals within species, then a procedure like 
#' the following should produce results:
#' 
#' rownames(Y) <- species
#' 
#' rdf <- rrpp.data.frame(Y = Y, subjects = species, x = x)
#' 
#' fit <- lm.rrpp.ws(Y ~ species * x, subject = species, data = rdf, Cov = myCov, ...)
#' 
#' where ... means other arguments.  The covariances in the the Covariance matrix can be 
#' sorted by the subjects factor but data will not be sorted.  Therefore, names matching
#' the subjects is essential.  Additionally, subjects must be a factor in the data frame
#' or a factor in the global environment.  It cannot be part of a list.  Something like
#' subjects <- mylist$species will not work.  Assuring that data and subjects are in the 
#' same \code{\link{rrpp.data.frame}} object as data is the best way to avoid errors.
#' 
#' Most attributes for this analysis are explained with \code{\link{lm.rrpp}}.  
#' The notable different attributes for this function are that: (1) a covariance 
#' matrix for the non-independence of subjects can be either a symmetric matrix 
#' that matches in dimensions the number of subjects or the number of observations; 
#' (2) a parameter (delta) that can range between near 0 and 1 to calibrate the 
#' covariances between observations of different subjects; and (3) a 
#' parameter (gamma) that is either 1 (equal) or the square-root of the subject 
#' sample size (sample) to calibrate the covariances among observations 
#' within subjects.  If delta = 0, it is expected that the covariance between
#' individual observations, between subjects, is the same as expected from the
#' covariance matrix, as if observations were the single observations made on subjects.
#' As delta approaches 1, the observations become more independent, as if it is 
#' expected that the many observations would not be expected to be
#' as correlated as if from one observation.  Increasing delta might be useful, if, 
#' for example, many individuals are sampled within species, from different locations,
#' different age groups, etc.  Alternatively, the sample size (n_i) for subject i 
#' can also influence the trend of inter-subject covariances.  If more individual
#' observations are sampled, the correlation between subjects might be favored
#' to be smaller compared to fewer observations.  The covariances can be adjusted
#' to allow for greater independence among observations to be assumed for larger samples.
#'
#' A design matrix, \bold{X}, is constructed with 0s and 1s to indicate subjects association,
#' and it is used to expand the covariance matrix (\bold{C}) by \bold{XCt(X)}, where \bold{t(X)}
#' is the matrix transpose.  The parameters in \bold{X} are multiplied by exp(-delta * gamma)
#' to scale the covariances.  (If delta = 0 and gamma = 1, they are unscaled.)  
#' 
#' These options for scaling covariances could be important for data with 
#' hierarchical organization.
#' For example, data sampled from multiple species with expected covariances 
#' among species based on phylogenetic distances, might be expected to not covary as
#' strongly if sampling encounters other strata like population, sex, and age.
#' An a priori expectation is that covariances among observations would be expected to
#' be smaller than between species, if only one observation per species were made.
#' 
#' If one wishes to have better control over between-subject and within-subject
#' covariances, based on either a model or empirical knowledge, a covariance matrix should 
#' be generated prior to analysis.  One can input a covariance matrix with dimensions 
#' the same as \bold{XCt(X)}, if they prefer to define covariances in an alternative way.
#' A function to generate such matrices based on separate
#' inter-subject and intra-subject covariance matrices is forthcoming.
#' 
#' IMPORTANT.  It is assumed that either the levels of the covariance matrix (if 
#' subject by subject) match the subject levels in the subject argument, or that
#' the order of the covariance matrix (if observation by observation) matches the order
#' of the observations in the data.  
#' No attempt is made to reorder a covariance matrix by observations
#' and row-names of data are not used to re-order the covariance matrix.  If the
#' covariance matrix is small (same in dimension as the number of subject levels), the function
#' will compile a large covariance matrix that is correct in terms of order, but
#' this is based on the subjects argument, only.
#' 
#' The covariance matrix is important for describing the expected covariances
#' among observations, especially knowing observations between and within subjects 
#' are not independent.  However, the randomization of residuals in a permutation
#' procedure (RRPP) is also important for testing inter-subject and 
#' intra-subject effects.  There are two RRPP philosophies used.  If the
#' variable for subjects is part of the formula, the subject effect is evaluated with
#' type III sums of squares and cross-products (estimates SSCPs between a model with all
#' terms and a model lacking subject term), and RRPP performed for all residuals of the
#' reduced model.  Effects for all other terms are evaluated with type II SSCPs and RRPP
#' restricted to randomization of reduced model residuals, within subject blocks.  This
#' assures that subject effects are held constant across permutations, so that intra-subject
#' effects are not confounded by inter-subject effects.
#' 
#' 
#' More details will be made and examples provided after publication of articles introducing 
#' the novel RRPP approach.
#' 
#' 
#' The \code{\link{lm.rrpp}} arguments not available for this function include: 
#' full.resid, block, and SS.type.  These arguments are fixed because of
#' the within-subject blocking for tests, plus the requirement for type II SS
#' for within-subject effects.
#' 
#' 
#' @param f1 A formula for the linear model (e.g., y~x1+x2).  
#' @param subjects A variable that can be found in the data frame indicating the research subjects
#' for the analysis.  This variable must be in the data frame.  Is can be either numeric 
#' (if its slot in the data frame is known) or a character, e.g., "sub_id".  It is imperative that
#' it is ordered the same as the data but that the data do not have row names the same as subjects.
#' For example, the subjects variable in the data frame might be sub_id: sub1, sub1, sub1, sub2,
#' sub2, sub2, ... and the row names of the data might be obs1, obs2, obs3, obs4, obs5, obs6, ...  
#' The data do not need to have row names but the subjects variable has to be provided.
#' @param iter Number of iterations for significance testing
#' @param turbo A logical value that if TRUE, suppresses coefficient estimation 
#' in every random permutation.  This will affect subsequent analyses that 
#' require random coefficients (see \code{\link{coef.lm.rrpp}})
#' but might be useful for large data sets for which only ANOVA is needed.
#' @param seed An optional argument for setting the seed for random 
#' permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found 
#' for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  
#' One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param int.first A logical value to indicate if interactions of first 
#' main effects should precede subsequent main effects
#' @param RRPP A logical value indicating whether residual randomization 
#' should be used for significance testing
#' @param Cov An optional argument for including a covariance matrix 
#' to address the non-independence
#' of error in the estimation of coefficients (via GLS).  If included, 
#' any weights are ignored.  This matrix must match in dimensions either
#' the number of subject levels or the number of observations.
#' @param delta A within-subject scaling parameter for covariances, ranging from 
#' 0 to 1.  If delta = 0, a sight value (0.001) is added to assure variances of the 
#' covariance matrix are 0.1 percent larger than covariances.
#' @param gamma A sample-size scaling parameter that is adjusted to be 1 ("equal")
#' scaling or the square-root of the sample size for subject observations ("sample").
#' @param data A data frame for the function environment, see 
#' \code{\link{rrpp.data.frame}}.  A data frame is required for this analysis.
#' @param print.progress A logical value to indicate whether a progress 
#' bar should be printed to the screen.
#' This is helpful for long-running analyses.
#' @param verbose A logical value to indicate if all possible output from an analysis 
#' should be retained. Generally this should be FALSE, unless one wishes to extract, 
#' e.g., all possible terms, model matrices, QR decomposition, or random permutation 
#' schemes.
#' @param Parallel Either a logical value to indicate whether parallel processing 
#' should be used, a numeric value to indicate the number of cores to use, or a predefined
#' socket cluster.  This argument defines parallel processing via the \code{parallel} library. 
#' If TRUE, this argument invokes forking or socket cluster assignment of all processor cores, 
#' except one.  If FALSE, only one core is used. A numeric value directs the number of cores to 
#' use, but one core will always be spared.  If a predefined socket cluster (Windows) is provided,
#' the cluster information will be passed to \code{parallel}.
#' @param ... Arguments typically used in \code{\link{lm}}, such as 
#' weights or offset, passed on to
#' \code{lm.rrpp} for estimation of coefficients.  If both weights and 
#' a covariance matrix are included,
#' weights are ignored (since inverses of weights are the diagonal elements 
#' of weight matrix, used in lieu
#' of a covariance matrix.)
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{lm.rrpp.ws} is a list containing the 
#' following
#' \item{call}{The matched call.}
#' \item{LM}{Linear Model objects, including data (Y), coefficients, 
#' design matrix (X), sample size
#' (n), number of dependent variables (p), dimension of data space (p.prime),
#' QR decomposition of the design matrix, fitted values, residuals,
#' weights, offset, model terms, data (model) frame, random coefficients 
#' (through permutations),
#' random vector distances for coefficients (through permutations), 
#' whether OLS or GLS was performed, 
#' and the mean for OLS and/or GLS methods. Note that the data returned 
#' resemble a model frame rather than 
#' a data frame; i.e., it contains the values used in analysis, which 
#' might have been transformed according to 
#' the formula.  The response variables are always labeled Y.1, Y.2, ..., 
#' in this frame.}
#' \item{ANOVA}{Analysis of variance objects, including the SS type, 
#' random SS outcomes, random MS outcomes,
#' random R-squared outcomes, random F outcomes, random Cohen's f-squared 
#' outcomes, P-values based on random F
#' outcomes, effect sizes for random outcomes, sample size (n), number of 
#' variables (p), and degrees of freedom for
#' model terms (df).  These objects are used to construct ANOVA tables.}
#' \item{PermInfo}{Permutation procedure information, including the number 
#' of permutations (perms), The method
#' of residual randomization (perm.method), and each permutation's sampling 
#' frame (perm.schedule), which
#' is a list of reordered sequences of 1:n, for how residuals were 
#' randomized.}
#' \item{Models}{Reduced and full model fits for every possible model 
#' combination, based on terms
#' of the entire model, plus the method of SS estimation.}
#' @seealso \code{\link{lm.rrpp}}; \code{\link{measurement.error}}
#' @references Adams, D.C and M.L Collyer. (submitted) Extended phylogenetic regression models for 
#' comparing within-species patterns across the Tree of Life. 
#' Methods in Ecology and Evolution
#' @references ter Braak, C.J.F. 1992. Permutation versus bootstrap significance tests in 
#' multiple regression and ANOVA. pp .79–86 In Bootstrapping and Related Techniques. eds K-H. Jockel, 
#' G. Rothe & W. Sendler.Springer-Verlag, Berlin.
#' \code{\link[stats]{lm}} for more on linear model fits.
#' @examples 
#' \dontrun{
#' data(fishy)
#' 
#' suppressWarnings(fit <- lm.rrpp.ws(coords ~ subj + groups * reps,
#'   subjects = "subj", 
#'   data = fishy))
#' 
#' anova(fit)
#' }

lm.rrpp.ws <- function(f1, subjects, 
                       iter = 999, turbo = FALSE, 
                       seed = NULL, int.first = FALSE,
                       RRPP = TRUE, 
                       data, Cov = NULL,
                       delta = 0.001, 
                       gamma = c("sample", "equal"),
                       print.progress = FALSE, 
                       verbose = FALSE,
                       Parallel = FALSE, ...) {
  
  L <- c(as.list(environment()), list(...))
  if(is.null(L$data))
    stop("\nWithin-subjects analysis requires a data frame, continaing the subjects factor.\n",
         call. = FALSE)
  if(is.numeric(subjects) && length(subjects) > 1)
    stop("\nThe subjects argument should be a single numeric or character values.\n",
         call. = FALSE)
  if(is.numeric(subjects)) subj.no <- subjects else
    subj.no <- which(names(L$data) %in% subjects)
  
  if(length(subj.no) == 0)
    stop("\nSubjects factor not found in data frame.\n",
         call. = FALSE)
  
  if(is.null(L$data))
    stop(paste("\nA data frame is required for this analysis,",
               "and it must contain the subjects factor.\n", sep = " "),
         call. = FALSE)
  L$sub.var.name <- names(L$data)[[subj.no]]
  L$subjects <- as.factor(L$data[[subj.no]])
  L$SS.type <- "IIws"
  L$full.resid <- FALSE
  L$gamma <- match.arg(gamma)
  L$block <- NULL
  if(length(L$gamma) > 1) L$gamma <- "equal"
  out <- do.call(.lm.rrpp, L)
  out$call <- match.call()
  out$ANOVA$SS.type <- "Within-subject II"
  class(out) <- c("lm.rrpp.ws", class(out))
  out

}

