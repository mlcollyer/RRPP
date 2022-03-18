## Resubmission
This is a patch release, 1.2.3, that fixes three small bugs, which were inadvertently induced by minor release 1.2.1.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.1.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There was a repeated note on the R-hub Builder for Windows Server that could not be resolved: 
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

