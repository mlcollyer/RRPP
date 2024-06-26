% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRPP.utils.r
\name{summary.pairwise}
\alias{summary.pairwise}
\title{Print/Summary Function for RRPP}
\usage{
\method{summary}{pairwise}(
  object,
  stat.table = TRUE,
  test.type = c("dist", "stdist", "mdist", "VC", "DL", "var"),
  angle.type = c("rad", "deg"),
  confidence = 0.95,
  show.vectors = FALSE,
  ...
)
}
\arguments{
\item{object}{Object from \code{\link{pairwise}}}

\item{stat.table}{Logical argument for whether results should be 
returned in one table 
(if TRUE) or separate pairwise tables (if FALSE)}

\item{test.type}{the type of statistic to test.  See below
should be used in the test.}

\item{angle.type}{If test.type = "VC", whether angle results are 
expressed in radians or degrees.}

\item{confidence}{Confidence level to use for upper confidence 
limit; default = 0.95 (alpha = 0.05)}

\item{show.vectors}{Logical value to indicate whether vectors 
should be printed.}

\item{...}{Other arguments passed onto pairwise}
}
\description{
See \code{\link{pairwise}} for further description.
}
\details{
The following summarize the test that can be performed: 

\describe{
\item{\bold{Distance between vectors, "dist"}}{ Vectors for LS means or 
slopes originate at the origin and point to some location, having both a 
magnitude
and direction.  A distance between two vectors is the inner-product of of 
the vector difference, i.e., the distance between their endpoints.  For
LS means, this distance is the difference between means.  For multivariate 
slope vectors, this is the difference in location between estimated change 
for the dependent variables, per one-unit change of the covariate considered.
For univariate slopes, this is the absolute difference between slopes.}
\item{\bold{Standardized distance between vectors, "stdist"}}{Same as the 
distance between vectors, but distances are divided by the model standard error 
(square-root of the trace of the residual covariance matrix). Pairwise tests
with this statistic should be consistent with ANOVA results.}
\item{\bold{Mahalanobis distance between vectors, "mdist"}}{Similar to the
standardized distance between vectors but the inverse of the 
residual covariance matrix is used in calculation of the distance, rather than
dividing the Euclidean distance between means and dividing by the standard error. 
If the residual covariance matrix is singular, Mahalanobis distances will not be 
calculated.  Pairwise tests with this statistic should be consistent with 
MANOVA results.}
\item{\bold{Vector correlation, "VC"}}{ If LS mean or slope vectors are 
scaled to unit size, the vector correlation is the inner-product of the 
scaled vectors.
The arccosine (acos) of this value is the angle between vectors, which 
can be expressed in radians or degrees.  Vector correlation indicates 
the similarity of 
vector orientation, independent of vector length.}
\item{\bold{Difference in vector lengths, "DL"}}{  If the length of a 
vector is an important attribute -- e.g., the amount of multivariate 
change per one-unit
change in a covariate -- then the absolute value of the difference in 
vector lengths is a practical statistic to compare vector lengths.  
Let d1 and
d2 be the distances (length) of vectors.  Then |d1 - d2| is a statistic 
that compares their lengths.}
\item{\bold{Variance, "var"}}{  Vectors of residuals from a linear 
model indicate can express the distances of observed values from 
fitted values.  Mean
squared distances of values (variance), by group, can be used to 
measure the amount of dispersion around estimated values for groups.  
Absolute
differences between variances are used as test statistics to compare 
mean dispersion of values among groups.  Variance degrees of freedom 
equal n, 
the group size, rather than n-1, as the purpose is to compare mean 
dispersion 
in the sample.  (Additionally, tests with one subject in a group 
are possible, or at least not a hindrance to the analysis.)}
}

The argument, \code{test.type} is used to select one of the tests 
above.  See \code{\link{pairwise}} for examples.

 \subsection{Notes for RRPP 0.6.2 and subsequent versions}{ 
 In previous versions of pairwise, \code{\link{summary.pairwise}} had three 
 test types: "dist", "VC", and "var".  When one chose "dist", for LS mean 
 vectors, the statistic was the inner-product of the vector difference.  
 For slope vectors, "dist" returned the absolute value  of the difference 
 between vector lengths, which is "DL" in 0.6.2 and subsequent versions.  This
 update uses the same calculation, irrespective of vector types.  Generally,
 "DL" is the same as a contrast in rates for slope vectors, but might not have
 much meaning for LS means.  Likewise, "dist" is the distance between vector
 endpoints, which might make more sense for LS means than slope vectors.  
 Nevertheless, the user has more control over these decisions with version 0.6.2
 and subsequent versions.
}

 \subsection{Notes for RRPP 2.0.4 and subsequent versions}{ 
 The test types, standardized distance between vectors, "stdist", and Mahalanobis distances
 between vectors were added.  The former
 simply divides the distance between vectors by model standard error (square-root 
 of the trace of the residual covariance matrix).  This is a multivariate generalization
 of a t-statistic for the difference between means.  It is not the same as Hotelling 
 squared-T-statistic, which requires incorporating the inverse of the residual covariance 
 matrix in the calculation of the distance, a calculation that also requires a non-singular
 covariance matrix.  However, the Mahalanobis distances are similar (and proportional) 
 to the Hotelling squared-T-statistic. Pairwise tests using Mahalanobis distances represent
 a non-parametric analog to the parametric Hotelling squared-T test.  Both tests should be 
 better for GLS model fits compared
 to Euclidean distances between means, as the total sums of squares are more likely to vary
 across random permutations. In general, if ANOVA is 
 performed a pairwise test with "stdist" is a good association; if MANOVA is performed,
 a pairwise test with "mdist" is a good association.
}
}
\author{
Michael Collyer
}
\keyword{utilities}
