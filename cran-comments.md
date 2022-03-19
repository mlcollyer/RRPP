## Resubmission
This is a patch release, 1.2.3, that fixes three small bugs, which were inadvertently induced by minor release 1.2.1.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release.  All other dependencies should be unaffected, functionally, as no attribute or functional changes have been introduced.

## Test environments
* local OS X install, R 4.1.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. There were two notes on the R-hub Builder for Windows Server that could not be resolved: 

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

The above note has appeared on this and the Rhub check for the previous release, 1.2.2.  It was not an issue on the previous release,

* checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                       user system elapsed
  lm.rrpp             17.69   2.52   49.39
  predict.lm.rrpp     11.45   1.51   27.91
  pairwise             7.37   0.80   20.89
  trajectory.analysis  6.20   0.97   13.25
  coef.lm.rrpp         4.42   1.02   13.02
  manova.update        3.88   0.27    9.22
  convert2ggplot       2.94   0.47    7.28
  
After receiving prep errors after three attempts, this note started to appear in multiple attempts, but does not seem to make much sense.   None of these examples have ever taken much time, previously, and these functions were not updated.  I believe this is a continuing issue inherent with the Windows server, as suggested by the prep errors that stopped without any package changes.
  
