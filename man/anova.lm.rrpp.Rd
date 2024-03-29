% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.lm.rrpp.r
\name{anova.lm.rrpp}
\alias{anova.lm.rrpp}
\title{ANOVA for lm.rrpp model fits}
\usage{
\method{anova}{lm.rrpp}(
  object,
  ...,
  effect.type = c("F", "cohenf", "SS", "MS", "Rsq"),
  error = NULL,
  print.progress = TRUE
)
}
\arguments{
\item{object}{Object from \code{\link{lm.rrpp}}}

\item{...}{Additional lm.rrpp model fits or other arguments passed to anova.}

\item{effect.type}{One of "F", "cohenf", "SS", "MS", "Rsq" to choose from 
which distribution of statistics to calculate effect sizes (Z).  
See \code{\link{lm.rrpp}}.}

\item{error}{An optional character string to define MS error term for 
calculation of F values. See \code{\link{lm.rrpp}} for examples.}

\item{print.progress}{A logical argument if multiple models are used and 
one wishes to view progress for sums of squares (SS) calculations.}
}
\description{
Computes an analysis of variance (ANOVA) table using 
distributions of random statistics from \code{\link{lm.rrpp}}.  
ANOVA can be performed on one model or multiple models.  
If the latter, the first model is considered a null model for 
comparison to other models.  The ANOVA is functionally similar to a 
non-parametric likelihood ratio test for all null-full model comparisons
Residuals from the null model will be used to generate random pseudo-values 
via RRPP for evaluation of subsequent models. The permutation schedule from 
the null model will be used for random permutations.
This function does not correct for improper null models.  One must assure 
that the null model is nested within the other models.  Illogical results 
can be generated if this is not the case.
}
\examples{
\dontrun{
# See examples for lm.rrpp to see how anova.lm.rrpp works in conjunction
# with other functions

data(Pupfish)
names(Pupfish)
Pupfish$logSize <- log(Pupfish$CS) # better to not have functions in formulas

# Single-Model ANOVA

fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", 
data = Pupfish, print.progress = FALSE, iter = 999) 
anova(fit)
anova(fit, effect.type = "MS")
anova(fit, effect.type = "Rsq")
anova(fit, effect.type = "cohenf")

# Multi-Model ANOVA (like a Likelihood Ratio Test)
fit.size <- lm.rrpp(coords ~ logSize, SS.type = "I", data = Pupfish, 
print.progress = FALSE, iter = 999) 
fit.sex <- lm.rrpp(coords ~ logSize + Sex, SS.type = "I", data = Pupfish, 
print.progress = FALSE, iter = 999) 
fit.pop <- lm.rrpp(coords ~ logSize + Pop, SS.type = "I", data = Pupfish, 
print.progress = FALSE, iter = 999) 
anova(fit.size, fit.sex, fit.pop, 
print.progress = FALSE) # compares two models to the first

# see lm.rrpp examples for mixed model ANOVA example and how to vary SS type
}

}
\author{
Michael Collyer
}
\keyword{utilities}
