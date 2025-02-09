# CHANGES IN RRPP 2.1.1
+ New function, `print.summary.betaTest`, to complement `summary.betaTest`.

### BUG FIXES
+ Updated, `summary.betaTest`, to provide invisible output, so that the summary 
table could be extracted without screen printing.
+ Added error catch for `.lm.rrpp` for mismatched numbers of subjects and data.

# CHANGES IN RRPP 2.1.0

### NEW FEATURES
+ New function, `kcomp`, which performs a K-component analysis.
+ New function, `betaTest`, which performs tests of coefficients vectors,
with enhanced flexibility compared to `coef.lm.rrpp`.

# CHANGES IN RRPP 2.0.4

### NEW FEATURES
+ New function, `mahal_dist`, which emulates the `dist` function but allows a covariance matrix to be used for generalized distances.
+ Added two tests to the `pairwise` suite of functions: `stdist` and `mdist`, for standardized and Malahanobis distances, respectively.
+ Added argument, `groups.first`, to `measurement.error` for ordering terms having subjects and groups in ANOVA.

### BUG FIXES
+ Fixed the missing re-order of data for `anc.BM`.
+ Added factor coercion to `measurement.error`.
+ Forced coercion of matrices for `getLSmeans`.
+ Adjusted sensitivity for `QRforX` for determining 0 values.
+ Fixed missing `Pcov` issue in `manova.update`.

# CHANGES IN RRPP 2.0.3

### NEW FEATURES
+ New function, `QRforX`, which produces consistent QR decomposition results despite differences between`base::qr` and `Matrix::qr` functions.

### BUG FIXES
+ Fixed an important issue with `removeRedundant` and `getRank` to properly remove axes based on a `base::qr` pivot strategy, but still avoiding the making of large, dense matrices before using `base::qr`.

### OTHER CHANGES
+ Expanded function tests for CRAN checks.
+ Updated `R`-dependency to >= 4.4.0

# CHANGES IN RRPP 2.0.2

### BUG FIXES
+ Fixed an important issue with `checkers` function that assumed sparse-matrix (`Matrix`) format of an object without coercing this format.

# CHANGES IN RRPP 2.0.1

### NEW FEATURES
+ Updated all QR decompositions that forced dense-matrix methods when not needed.
+ Added functions to remove design matrix parameter redundancies without using QR decomposition.

### BUG FIXES
+ Fixed issues with model matrices for new data in `predict.lm.rrpp` and made the code more universal for different model formulae.
+ Fixed standard error method for intercepts in `.coef.lm.rrpp`.
+ Fixed some circularity issues with `.lm.rrpp` associated with covariance matrices (causing non-adjusted matrices).
+ Fixed some issues with multiple functions that misused the logical `subTest` output from `lm.rrpp.ws`.
+ Fixed issues with `anova.multi.model` and `logL` related to missing `Pcov` when GLS models are used.

# CHANGES IN RRPP 2.0.0

### NEW FEATURES
* Added arguments for restricted resampling to `lm.rrpp` and `predict.lm.rrpp`, and permutation of full model residuals (with restrictions) to `lm.rrpp`.
* Added new function, `measurement.error`.
* Added new function, `lm.rrpp.ws`.
* Added `verbose` arguments to most analytical functions to reduce dense results.
* Added utility functions (`getANOVAStats`, `getPermInfo`, `getTerms`, and `getModels`) to accommodate `verbose` function choices, plus allow users to easily extract objects as output.
* Added new function, `ICCstats`, for calculating ICC statistics for `lm.rrpp` fits.
* Added S3 generic functions for `measurement.error`, but also added a plot utility function, `focusMEonSubjects`, which builds off of `plot.measurement.error`, to focus attention on specific plot details.
* Added function, `interSubVar`, a utility function to use along with `measurement.error`,
to visualize how variation between replicate measures, within subjects, can affect inter-subject variation.

### BUG FIXES
*  Fixed SSCP issue in `summary.lm.rrpp` that produces matrices of 0s.  The issue was inherent to `lm.rrpp`, from erroneous assignment of models, only for output (not for mechanics).
* Attempted fix of `predict.lm.rrpp` for QR-truncated design matrices. 
* Fixed illogical `print.progess` conditions in `pairwise`.
* Fixed potential problem with reliance on `Matrix::qrQ`, using rather, `Matrix::qr.Q`.
* Fixed GLS estimation issue with `pairwise`.

### OTHER CHANGES
* argument `SE` added to `coef.lm.rrpp` for calculating bootstrap standard error.

# CHANGES IN RRPP 1.3.1

### NEW FEATURES
*  Added ability to use either Pearson or normalized residuals in residual plots for `plot.lm.rrpp`
* Updated arguments for `plot.ordinate` to better work with `plot.default`.

### BUG FIXES
*  Fixed S4 issue with GLS for `plot.lm.rrpp` diagnostic plots.

### OTHER CHANGES

# CHANGES IN RRPP 1.3.0

### NEW FEATURES
* New S3 function, `logLik.lm.rrpp` to obtain log-likelihood from an `lm.rrpp` object.
* New function, `scaleCov` to scale covariance matrices with linear or exponential scalars.

### BUG FIXES
*  Various lingering bugs from version 1.2 associated with sparse matrix calculation issues were removed.

### OTHER CHANGES
* Updated `model.comparison` to include Z scores calculated from log-likelihoods.

# CHANGES IN RRPP 1.2.3

### BUG FIXES
* Fixed a bug in `predict.lm.rrpp` for new data frames with only one observation.
* Fixed a bug in `SS.iter.main` that incidentally wrapped RSS.model by rows rather than columns.
* Fixed a bug in `summary.lm.rrpp` that did not properly index a matrix of `RSS.model`.

# CHANGES IN RRPP 1.2.2

### BUG FIXES
*  Fixed a typo in `plot.lm.rrpp` for diagnostic plots, which forced an error.

# CHANGES IN RRPP 1.2.1

### NEW FEATURES
* Updated many functions to be able to use `Matrix` class matrices for more efficient computation, when needed.
* Updated `checkers` function to use better algorithms to switch among different class matrices, to better save memory and increase computation time efficiency.

### BUG FIXES
*  Added catch to `anc.BM` for singleton nodes.

### OTHER CHANGES
* Updated `lm.rrpp` to have less detritus during use.  Also adjusted/updated support functions to work with updates.

# CHANGES IN RRPP 1.1.2

### BUG FIXES
* Fixed a few issues with dropped names when working with `data.frame` objects.

# CHANGES IN RRPP 1.1.1

### BUG FIXES
* Fixed issue in `lm.args.from.formula` for intercept only models with covariance matrices.
* Fixed issue with non-full rank design matrices and dropped terms in `lm.rrpp` support functions.
* Fixed issue in `logL` support function for non-full rank design matrices.

# CHANGES IN RRPP 1.1.0

### NEW FEATURES
* `na.omit.rrpp.data.frame` added for handling missing data.
* `looCV` function added as diagnostic tool.

### OTHER CHANGES
* `coef.lm.rrpp` updated to provide results based on `SS.type`, rather than
type III SS, only.
* Complete scanning and updating of code to identify and fix names lost due to matrix calculations or vector/matrix conversions in `R`.

### BUG FIXES
* Fixed missing coefficients for GLS estimation in `coef.lm.rrpp`
* Fixed missing data issue in `lm.rrpp` (passed onto `lm`).
* Fixed some `as.matrix` names dropping in support code.

# CHANGES IN RRPP 1.0.0 (Major Release)

### NEW FEATURES
* Computation of sums of squares in `lm.rrpp` using `Matrix` package and sparse matrix calculations.  This speeds up computation time and requires far less memory allocation.
* The `lm.rrpp` function now has a an argument, `turbo`, which can suppress calculating coefficients in random permutations, if unneeded, which can speed up analysis of large data sets.
* More options for parallel processing.
* New function, `convert2ggplot`, for coercing `RRPP` plots into ggplot objects.

# CHANGES IN RRPP 0.6.2 (Patch Release)

### NEW FEATURES
* Added option to add abscissa to `plot.predict.lm.rrpp`.
* Added Box-Cox transformation to `effect.size`.
* Added option to flip axes in `plot.ordinate`.
* Updated `predict.lm.rrpp` so that functions in formulae are permissible.

### OTHER CHANGES 
* Adapted `summary.pairwise` to perform degree transformations rather than `print.summary.pairwise`, so that objects saved are the same as objects printed.
* Removed error trap from pairwise for n = 1 groups.

### BUG FIXES
* Fixed class check on `model.comparison`.
* Fixed `manova.update` pc dimension issue (output)
* Fixed `xlim` and `ylim` to be adjustable in `plot.ordinate`.
* Added transform logical output to ordinate.
* Fixed intercept only model issue in lm.rrpp.

# CHANGES IN RRPP 0.6.1 (Patch Release)

### NEW FEATURES
### OTHER CHANGES
### BUG FIXES
* Refined `manova.update` to have more efficient code and better notes.
* Fixed a bug in `summary.manova.lm.rrpp` that mixed up rows and columns of a matrix of random stats.
* Fixed a bug that eliminated row names for `lm.rrpp` data when converting a vector of data to a matrix.
* Added a step to coerce single-column matrices to vectors in `rrpp.data.frame` to prevent downstream issues.
* Made it possible to use `trajectory.analysis` with univariate response data, omitting vector correlations output.

# CHANGES IN RRPP 0.6.0 (Minor Release)

### NEW FEATURES
* `ordinate` function.
* `summary.ordinate` and `plot.ordinate` S3 functions
* `add.tree` function (for plotting with `plot.ordinate`)

### OTHER CHANGES
* Updated support functions for `lm.rrpp` to provide better flexibility for different formulas.  
* Changed `$LM$data` in `lm.rrpp` to be a model frame rather than a data frame, consistent with `$model` from `lm`.

### BUG FIXES 
* Fixed issue with transformation of dependent variable -- e.g., log(y) -- in formula for lm.rrpp. The bug was a failure to transform the variable(s).

# CHANGES IN RRPP 0.5.2 (Patch Release)

### BUG FIXES 
* Fixed issue with GLS estimation of coefficients in pairwise (used GLS residuals instead of
transformed residuals, by mistake).

# CHANGES IN RRPP 0.5.1 (Patch Release)

### BUG FIXES 
* Fixed data.frame error in lm.rrpp for within-formula transformations (like log or poly).
* Arranged effect.type argument options in lm.rrpp to match anova.lm.rrpp).

# CHANGES IN RRPP 0.5.0 (Minor Release)

### NEW FEATURES
* `prep.lda` A new function to generate arguments for `lda` in the `MASS` library.

### OTHER CHANGES
* `classify` deprecated (in favor of `prep.lda`)
* Changed the output structure for lm.rrpp objects.  Either OLS or GLS statistics are returned, not both.  There is no longer differentiation between GLS and weighted stats, as weighted stats are GLS estimated stats.  The difference is obvious with the inclusion of weights or Cov matrix.
* No longer a `verbose` option with `manova.update`.  The function was optimized to provide verbose output without having to slow down computation time.

### BUG FIXES 
* Fixed some issues with linear model weights and offset, throughout many functions.
* Sweeping changes to be CRAN 4.0 compliant



# CHANGES IN RRPP 0.4.3 (Patch Release)

### NEW FEATURES

### OTHER CHANGES

### BUG FIXES 
* Fixed a parallel processing issue with internal function, `SS.iter` (produced incorrect RSS).
* Fixed multi-model anova TSS issue.
* Fixed an issue with the centering of trajectories in `trajectoty.analysis`.
* Fixed an issue for calculating fitted values for GLS in the `lm.rrpp` function.


------

# CHANGES IN RRPP 0.4.2 (Patch Release)

### NEW FEATURES

### OTHER CHANGES
* Added `tol` and `pc.no` arguments to `model.comparison` (were fixed before) so that users have more control of the analysis.

### BUG FIXES 
* Tweaked logL tolerance calculation to be consistent with `prcomp`.
* Tweaked support code for `lm.rrpp` to work better with missing data frames.
* Fixed absolute eigenvalue issue with Cov.proj.
* Added forgotten code from last update to fix non-full rank design matrices.
* Fixed `trajectory.analysis` traj.list issue, to not use grep for sorting trajectories.  (Now lexical ordering of interactions is used.)
* Changed `det` to `determinant` in all needing functions, to use modulus for near-singular matrices
* Fixed bug in `plot.lm.rrpp` (code lines out of order)

------

# CHANGES IN RRPP VERSION 0.4.1 (Patch release)

### NEW FEATURES

### OTHER CHANGES
* Updated source in Description.

### BUG FIXES 
* Fixed univariate data issue with `model.comparison`
* Fixed xlab flexibility issue for regression plots in `plot.procD.lm`
    
------

# CHANGES IN RRPP VERSION 0.4.0 (Minor release)

### NEW FEATURES
* `manova.update` function
* `trajectory.analysis` function
* `reveal.model.designs` function
* New vignette for `ANOVA versus MANOVA in RRPP`

### OTHER CHANGES
* Added a pairwise variance comparison for the `pairwise` function.

### BUG FIXES 
* Tuned F-stat calculations to allow for model-specific residual variances, for multiple terms.
* Updated `procD.lm` to better work with data in the global environment rather than a data frame. 

------

# CHANGES IN RRPP VERSION 0.3.0 (Minor release)

### NEW FEATURES
* Added `print.summary.pairwise`.
* Added `model.comparison` function.
* Added classify function.

### OTHER CHANGES

* Added an update to allow classify to work on univariate data.

### BUG FIXES 

* Tuned criterion for assessing whether generalized inverse should be performed.
* Fixed bug with GLS variance estimation in pairwise.
* Fixed some issues with univariate data for `classify`.
* Fixed some issues with univariate data for `pairwise`.
* Fixed the logL function within model.comparisons for GLS determinants 
    (was returning 0).
* Fixed some issues with the aov.multimodel subfunction of `anova.lm.rrpp`, related to GLS permutations and intercept only models.  
* Added random SS output to aov.multimodel subfunction of `anova.lm.rrpp`, so that it can be called by other functions/packages.
    
------

# CHANGES IN RRPP VERSION 0.2.0 (Minor release)

### NEW FEATURES

* `pairwise function`: allows pairwise comparison of means or slopes for a `lm.rrpp fit`.
* A vignette for using RRPP, which is the same as Appendix S2 in Collyer and Adams 
    (2008). RRPP: An R package for fitting linear models to high-dimensional data using 
    residual randomization. Methods in Ecology and Evolution.  (submitted)

### OTHER CHANGES

* Added multi-model inference capability to `anova.lm.rrpp`
    
### BUG FIXES 

* Fixed issue for `coef.lm.rrpp` tests when type II or type III SS is chosen, to make sure that appropriate coefficients are used.
 
------

# CHANGES IN RRPP VERSION 0.1.0 (Major release)

###	New Release!

###	New features

*	`anova.lm.rrpp.r`
*   `coef.lm.rrpp.r`
*	`lm.rrpp.r`
*	`predict.lm.rrpp.r`
*	`RRPP.support.code.r`
*	`RRPP.utils.r`

* Added a `NEWS.md` file to track changes to the package.
