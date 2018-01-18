## Test environments
* local OS X install, R 3.4.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michael Collyer <mlcollyer@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Evlauation (3:21)
  radomization (8:95)

** running examples for arch 'i386' ... [19s] NOTE
Examples with CPU or elapsed time > 10s
         user system elapsed
lm.rrpp 11.25   0.15   11.64

** running examples for arch 'x64' ... [22s] NOTE
Examples with CPU or elapsed time > 10s
         user system elapsed
lm.rrpp 13.56   0.22    13.9
  
## Updates
\
I have fixed the typos and removed on-screen printing from examples, plus reduced iterations to shorten elapsed time.
