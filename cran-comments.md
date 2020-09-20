## Resubmission
This is a patch release, version 0.6.2, that fixes a few bugs in version 0.6.1 (which was a patch release).   The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this patch release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.0.2
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There was a NOTE:

  checking for future file timestamps ... NOTE
  unable to verify current time
  
Research on this issue suggests that devtools::check is using an external web resource that is no longer active.