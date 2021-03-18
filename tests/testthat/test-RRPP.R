### -----------------------------------------------------------------------
###
### All parallel processing options must be run outside of CRAN check
###
### All tests that require dependencies should be run outside of CRAN check
###
### -----------------------------------------------------------------------


### lm.rrpp ---------------------------------------------------------------

test_that("fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                                print.progress = FALSE, iter = 3))
})

test_that("fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3))
})

test_that("fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS, data = Pupfish, 
                  print.progress = FALSE, iter = 3))
})

test_that("fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  print.progress = FALSE, iter = 3))
})

test_that("fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "II", print.progress = FALSE, iter = 3))
})

test_that("fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "III", print.progress = FALSE, iter = 3))
})

test_that("fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  turbo = TRUE, print.progress = FALSE, iter = 3))
})


test_that("fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  print.progress = FALSE, iter = 3))
})

test_that("fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  int.first = TRUE,
                  print.progress = FALSE, iter = 3))
})

# needed: Cov application

### coef.lm.rrpp ----------------------------------------------------------

test_that("coef.fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3)))
})

test_that("coef.fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3)))
})

test_that("coef.fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS, data = Pupfish, 
                  print.progress = FALSE, iter = 3)))
})

test_that("coef.fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  print.progress = FALSE, iter = 3)))
})

test_that("coef.fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "II", print.progress = FALSE, iter = 3)))
})

test_that("coef.fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "III", print.progress = FALSE, iter = 3)))
})

test_that("coef.fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  turbo = TRUE, print.progress = FALSE, iter = 3)))
})


test_that("coef.fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  print.progress = FALSE, iter = 3)))
})

test_that("coef.fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  int.first = TRUE,
                  print.progress = FALSE, iter = 3)))
})


test_that("coef.t.fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS, data = Pupfish, 
                       print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "II", print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "III", print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       turbo = TRUE, print.progress = FALSE, iter = 3),
          test = TRUE))
})



test_that("coef.t.fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       print.progress = FALSE, iter = 3),
          test = TRUE))
})

test_that("coef.t.fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       int.first = TRUE,
                       print.progress = FALSE, iter = 3),
          test = TRUE))
})

# needed: Cov application

### anova.lm.rrpp ----------------------------------------------------------


test_that("anova.fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3)))
})

test_that("anova.fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3)))
})

test_that("anova.fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS, data = Pupfish, 
                       print.progress = FALSE, iter = 3)))
})

test_that("anova.fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       print.progress = FALSE, iter = 3)))
})

test_that("anova.fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "II", print.progress = FALSE, iter = 3)))
})

test_that("anova.fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "III", print.progress = FALSE, iter = 3)))
})

test_that("anova.fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       turbo = TRUE, print.progress = FALSE, iter = 3)))
})


test_that("anova.fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       print.progress = FALSE, iter = 3)))
})

test_that("anova.fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       int.first = TRUE,
                       print.progress = FALSE, iter = 3)))
})

test_that("anova.fit11.b.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3),
                effect.type = "cohenf"))
})

test_that("anova.fit11.c.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3),
                effect.type = "Rsq"))
})

test_that("anova.fit11.d.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3),
                effect.type = "MS"))
})

test_that("anova.fit11.e.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3),
                error = c("Pop:Sex", "Pop:Sex", "Residuals")))
})

# needed: Cov application

### predict.lm.rrpp -------------------------------------------------------
### manova.update ---------------------------------------------------------
### pairwise --------------------------------------------------------------
### trajectory.analysis ---------------------------------------------------
### ordinate --------------------------------------------------------------
### plot functions --------------------------------------------------------
### -----------------------------------------------------------------------