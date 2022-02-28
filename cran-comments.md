## Resubmission
This is a minor release, 1.2.1, that mostly improves computational efficiency with functions that identify when a switch to `Matrix` class sparse matrices can be made.  There was also one small bug fix.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.1.2
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There was a repeated note about a potential invalid URL.  The URL (journal article) is valid but service might be temporarily unavailable.

