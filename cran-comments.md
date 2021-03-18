## Resubmission
This is a major release, 1.0.0, that incorporates new `Rcpp` and `RcppArmadillo` functions for calculating sums of squares via `C++`, increasing speed and reducing memory allocation for large data sets.   Additionally, a new function, `convert2ggplot`, makes it possible to convert package plot objects to `ggplot` objects.  This function forces a dependency on `ggplot2` and that packages dependencies, therein.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

`RRPP 1.0.0` and `geomorph 4.0.0` are to be concurrently updated and released.

## Test environments
* local OS X install, R 4.0.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

