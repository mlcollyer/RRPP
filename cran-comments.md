## Resubmission
This is a major release, 2.0.0, which adds several new functions, options for existing functions, and some function efficiencies.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release and any necessary changes will be concomitantly submitted. 

## Test environments
* local OS X install, R 4.3.2
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There were three notes that could not be resolved.  These notes have appeared on previous submissions, as well.

Two on the R-hub Builder for the Windows Server:

* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ''NULL''

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

One on both the Fedora Linux and Ubuntu Linux servers:

Skipping checking HTML validation: no command 'tidy' found

