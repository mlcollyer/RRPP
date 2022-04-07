## Resubmission
This is a patch release, 1.3.0, that fixes a few small bugs caused by using `Matrix` class sparse matrices, plus adds two new functions: `logLik.lm.rrpp` and `scaleCov`.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.1.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There were two notes on the R-hub Builder for Windows Server that could not be resolved: 

