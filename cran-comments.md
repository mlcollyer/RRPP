## Resubmission

Initial submission was flagged with a note:

Found the following Rd file(s) with Rd \link{} targets missing package
anchors:
  convert2ggplot.Rd: ggplot
Please provide package anchors for all Rd \link{} targets not in the
package itself and the base packages.

These links were updated to include package anchors, and this is a re-submission of 2.0.3.

This is a patch release, 2.0.3, which adds a function, `QRforX`, which assures consistent QR results, whether obtained from `base::qr` or `Matrix::qr`.  Some additional tests have also been included.  The package is also updated to require R >= 4.4.0, to be consistent with `Matrix` dependency.

The package, `geomorph`, is developed by the same authors as `RRPP`.  The reverse dependency has been addressed during development of this release and any necessary changes will be concomitantly submitted. 

## Test environments
* local OS X install, R 4.4.0
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.