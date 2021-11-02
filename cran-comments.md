## Resubmission
This is a patch release, 1.1.2, that includes bug fixes for bugs incurred with minor release 1.1.0.  One function, `rrpp.data.frame` also now has an additional and optional argument for including observation names to assure functionality in functions that use `rrpp.data.frame`.  Functions were also tested to be compatible with R 4.1.2.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.1.2
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

