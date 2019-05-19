# RRPP VERSION 0.4.2 (Patch release)

### NEW FEATURES

### OTHER CHANGES
* Added `tol` and `pc.no` arguments to `model.comparison` (were fixed before) so that users have more control of the analysis.

### BUG FIXES 
* Tweaked logL tolerance (0.001) calculation to be consistent with prcomp.
* Tweaked support code for lm.rrpp to work better with missing data frames.
* Fixed eigen absolute eigenvalue issue with Cov.proj.
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
* Updated `procD.lm` to better work with data in the gloabl environment rather than a data frame. 

------

# CHANGES IN RRPP VERSION 0.3.0 (Minor release)

### NEW FEATURES
* Added `print.summary.pairwise`.
* Added `model.comparison` function.
* Added classify function.

### OTHER CHANGES

* Added an update to allow classify to work on univeraite data.

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

* Added multi-model inference capabaility to `anova.lm.rrpp`
    
### BUG FIXES 

* Fixed issue for `coef.lm.rrpp` tests when type II or type III SS is chosen, to make sure that appropriate coefficients are used.
 
------

# CHANGES IN RRPP VERSION 0.1.0 (Major release)

###	New Release!

###	New features

*	`anova.lm.rrpp.r`
* `coef.lm.rrpp.r`
*	`lm.rrpp.r`
*	`predict.lm.rrpp.r`
*	`RRPP.support.code.r`
*	`RRPP.utils.r`
