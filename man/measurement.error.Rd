% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measurement.error.r
\name{measurement.error}
\alias{measurement.error}
\title{Evaluation of measurement error for two or more multivariate measurements, 
for common research subjects}
\usage{
measurement.error(
  Y,
  subj,
  reps,
  groups = NULL,
  iter = 999,
  seed = NULL,
  multivariate = FALSE,
  use.PCs = TRUE,
  tol = 0.001,
  Parallel = TRUE,
  print.progress = FALSE
)
}
\arguments{
\item{Y}{A matrix (n x p) of data for n observations and p variables.}

\item{subj}{A vector or factor of research subjects (each subject should occur twice or more).
The length of the vector must equal the number of observations and will be coerced into a factor.}

\item{reps}{A vector or factor for replicate measurements for research subjects.
The length of the vector must equal the number of observations and will be coerced into a factor.}

\item{groups}{An optional vector, coercible to factor, to be included in the linear model
(as an interaction with replicates)..
This would be of interest if one were concerned with systematic ME occurring perhaps differently among 
certain strata within the data.  For example, systematic ME because of an observer bias might
only be observed with females or males.}

\item{iter}{Number of iterations for significance testing}

\item{seed}{An optional argument for setting the seed for random 
permutations of the resampling procedure.
If left NULL (the default), the exact same P-values will be found 
for repeated runs of the analysis (with the same number of iterations).
If seed = "random", a random seed will be used, and P-values will vary.  
One can also specify an integer for specific seed values,
which might be of interest for advanced users.}

\item{multivariate}{Logical value for whether to include multivariate analyses.  Intraclass correlation 
matrices and relative eigenanalysis are based on products of sums of squares and cross-products (SSCP)
matrices, some of which must be inverted and potentially require
 significant computation time.  If FALSE, only statistics based on dispersion of values are calculated.}

\item{use.PCs}{A logical argument for whether to use the principal components of the data.
This might be helpful for relative eigenanalysis, and if p > n, 
in which case inverting singular covariance matrices would not be possible.}

\item{tol}{A value indicating the magnitude below which 
components should be omitted., if use.PCs is TRUE. (Components are omitted if their 
standard deviations are less than or equal to tol times the 
standard deviation of the first component.)  See \code{\link{ordinate}} for more details.}

\item{Parallel}{The same argument as in \code{\link{lm.rrpp}} to govern parallel processing (
either a logical vale -- TRUE or FALSE -- or the number of threaded cores to use).  See \code{\link{lm.rrpp}} 
for additional details.}

\item{print.progress}{A logical value to indicate whether a progress 
bar should be printed to the screen.}
}
\value{
Objects of class "measurement.error" return a list of the following:
 \item{AOV}{Analysis of variance to test for systematic error, based on dispersion of values.}
 \item{mAOV}{Multivariate AOV based on product of the inverse of the random component (SSCP) of ME
 times the systematic component of ME.}
 \item{icc}{The intraclass correlation coefficient (ICC) based on the dispersion of values.}
 \item{mult.icc.eigs}{The eigenvalues of ICC matrices, culled to principal dimensions with positive eigenvalues.}
 \item{SSCP}{The sums of squares and cross-products matrices for model effects.}
 \item{SSCP.ME.product}{The products of the inverse of the random ME SSCP and the SSCP matrices
 for systematic ME,.  These are the same matrix products used for eigenanalysis.  
 This is the observed matrix.}
 \item{SSCP.ME.product.std}{A list of the symmetric forms of standradized SSCP.ME.products 
 that yield orthogonal eigenvectors.}
 \item{all.stats}{All SS, MS, eigen values, etc., from the RRPP analyses performed.  This is the same
 as the output found in an \code{\link{lm.rrpp}} object, updated with \code{\link{manova.update}}.}
}
\description{
Function performs analyses concerned with the repeatability (reliability) of multivariate data 
(measurements) collected from the same research subjects.  Although there is no
requirement for repeated measurements on all research subjects, the analysis assumes
that multiple observations are made.
}
\details{
This function performs analyses as described in XXX to assess systematic and random components of 
measurement error (ME).  It basically performs ANOVA with RRPP,
but with different restricted randomization strategies.  The reliability of research subject variation 
can be considered by restricting randomization within replicates; the consistency of replicate measures
can be considered by restricting randomization within subjects.
Inter-subject variation remains constant across all random permutations within subjects and 
inter-replicate variation remains constant across all random permutations within replicates.  Type II
sums of squares and cross-products (SSCP) are calculated to assure conditional estimation.

The results include univariate-like statistics based on dispersion of values and
eigenanalysis performed on a signal to noise matrix product of SSCP matrices 
(sensu Bookstein and Mitteroecker, 2014) 
including the inverse of the random component of ME and the systematic
component of ME.  The multivriate test is a permutational form of multivariate ANOVA (MANOVA).

The first eigenvalue (Roy's maximum root) is used as a test statistic 
for tests using MANOVA.
Intraclass correlation coefficients (ICC) are also calculated, both based on dispersion of values and 
covariance matrices, as descriptive statistics.  Multivariate generalizations of the statistics
described by Liljequist et al. (2019) are also used, along with eigenanalysis.  
Three statistics describe the ICC for the population,
agreement of measurements among subjects, and consistency between measurements.  The last statistic does not 
necessarily measure the sameness between measurements but the consistency of change between measurements,
which might be indicative of a systematic measurement error.  If groups are used, these three statistics are 
repeated, using the SSCP for groups-adjusted data.  This approach accounts for group differences,
which would avoid large subject variation compared to measurement error inflating ICC values.  If there are 
inherently disparate groups from which subjects are sampled, this approach can elucidate better agreement and 
consistency in light of group differences.
}
\examples{
# TBD
}
\references{
yet to be determined.

Bookstein, F. L., & Mitteroecker, P. (2014). Comparing covariance matrices by relative eigenanalysis, 
with applications to organismal biology. Evolutionary biology, 41(2), 336-350.

Liljequist, D., Elfving, B., & Skavberg Roaldsen, K. (2019). Intraclass correlation–A discussion 
and demonstration of basic features. PloS one, 14(7), e0219854.
}
\seealso{
\code{\link{lm.rrpp}}, \code{\link{manova.update}}
}
\author{
Michael Collyer
}
\keyword{analysis}