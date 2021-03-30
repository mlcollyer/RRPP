## Resubmission
<<<<<<< HEAD
This is a patch release, version 0.6.2, that fixes a few bugs in version 0.6.1 (which was a patch release).   The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this patch release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.
=======
This is a major release, 1.0.0, that incorporates a `Matrix` package dependency, increasing speed and reducing memory allocation for large data sets.   Additionally, a new function, `convert2ggplot`, makes it possible to convert package plot objects to `ggplot` objects.  This function forces a dependency on `ggplot2` and its packages dependencies, therein.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

`RRPP 1.0.0` and `geomorph 4.0.0` are to be concurrently updated and released.

A previous attempt to submit this package to `CRAN`, requiring compilation and dependency on `Rcpp` and `RcppArmadillo` produced an installation error despite passing all tests, previously.  This version has removed any need for compilation.
>>>>>>> develop

## Test environments
* local OS X install, R 4.0.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

