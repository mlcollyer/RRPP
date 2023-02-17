## Resubmission
This is a patch release, 1.3.2, that adds a new function plus S3 generic functions.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.2.1
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There were two notes that could not be resolved, one on the R-hub Builder for the Windows Server:

'lastMiKTeXException'

and one on the Fedora Linux server:

Skipping checking HTML validation: no command 'tidy' found

