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
                                print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "II", print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "III", print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  turbo = TRUE, print.progress = FALSE, iter = 3, verbose = TRUE))
})


test_that("fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  int.first = TRUE,
                  print.progress = FALSE, iter = 3, verbose = TRUE))
})

test_that("fit12.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, iter = 3, verbose = TRUE))
})


test_that("fit13.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(lm.rrpp(TailLength ~ 1, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE))
})

test_that("fit14.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE))
})

test_that("fit15.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(lm.rrpp(cbind(TailLength, BodyWidth) ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE))
})

test_that("fit16.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(lm.rrpp(TailLength ~ SVL + 0, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE))
})

test_that("fit17.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(lm.rrpp(cbind(TailLength, BodyWidth) ~ SVL + 0, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE))
})

### coef.lm.rrpp ----------------------------------------------------------

test_that("coef.fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("coef.fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ 1, data = Pupfish, 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS, data = Pupfish, 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "II", 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "III", 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  turbo = TRUE, 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})


test_that("coef.fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  int.first = TRUE,
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})


test_that("coef.fit12.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  iter = 3, verbose = TRUE)))
})


test_that("coef.fit13.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(TailLength ~ 1, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit14.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)))
})

test_that("coef.fit15.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(cbind(TailLength, BodyWidth) ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)))
})

test_that("coef.t.fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "II", print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "III", print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       turbo = TRUE, print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})



test_that("coef.t.fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(coef(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       int.first = TRUE,
                       print.progress = FALSE, iter = 3, verbose = TRUE),
          test = TRUE))
})

test_that("coef.t.fit12.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, iter = 3, verbose = TRUE),
               test = TRUE))
})


test_that("coef.t.fit13.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(TailLength ~ 1, data = PlethMorph, 
                       print.progress = FALSE, 
                       Cov = PlethMorph$PhyCov,
                       iter = 3, verbose = TRUE),
               test = TRUE))
})

test_that("coef.t.fit14.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                       print.progress = FALSE, 
                       Cov = PlethMorph$PhyCov,
                       iter = 3, verbose = TRUE),
               test = TRUE))
})

test_that("coef.t.fit15.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(coef(lm.rrpp(cbind(TailLength, BodyWidth) ~ SVL, data = PlethMorph, 
                       print.progress = FALSE, 
                       Cov = PlethMorph$PhyCov,
                       iter = 3, verbose = TRUE),
               test = TRUE))
})

# needed: Cov application

### betaTest---------------------------------------------------------------

test_that("beta01.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta02.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ 1, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
})

test_that("beta03.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta04.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta05.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "II", print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta06.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  SS.type = "III", print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta07.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                  turbo = TRUE, print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})


test_that("beta10.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta11.works", {
  library(RRPP)
  data("Pupfish")
  mod <- lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                  int.first = TRUE,
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta12.works", {
  library(RRPP)
  data("PlethMorph")
  mod <- lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})


test_that("beta13.works", {
  library(RRPP)
  data("PlethMorph")
  mod <- lm.rrpp(TailLength ~ 1, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))

})

test_that("beta14.works", {
  library(RRPP)
  data("PlethMorph")
  mod <- lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta15.works", {
  library(RRPP)
  data("PlethMorph")
  mod <- lm.rrpp(cbind(TailLength, BodyWidth) ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta16.works", {
  library(RRPP)
  data("PlethMorph")
  mod <- lm.rrpp(TailLength ~ SVL + 0, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

test_that("beta17.works", {
  library(RRPP)
  data("PlethMorph")
  mod <- lm.rrpp(cbind(TailLength, BodyWidth) ~ SVL + 0, data = PlethMorph, 
                  print.progress = FALSE, 
                  Cov = PlethMorph$PhyCov,
                  iter = 3, verbose = TRUE)
  
  succeed(betaTest(mod))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p)))
  succeed(betaTest(mod, Beta = rep(1, mod$LM$p), include.md = TRUE))
})

### anova.lm.rrpp ----------------------------------------------------------


test_that("anova.fit01.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords[,1] ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit02.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ 1, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit03.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit04.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit05.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "II", print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit06.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       SS.type = "III", print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit07.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS + Pop, data = Pupfish, 
                       turbo = TRUE, print.progress = FALSE, iter = 3, verbose = TRUE)))
})


test_that("anova.fit10.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit11.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                       int.first = TRUE,
                       print.progress = FALSE, iter = 3, verbose = TRUE)))
})

test_that("anova.fit11.b.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3, verbose = TRUE),
                effect.type = "cohenf"))
})

test_that("anova.fit11.c.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3, verbose = TRUE),
                effect.type = "Rsq"))
})

test_that("anova.fit11.d.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ CS * Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3, verbose = TRUE),
                effect.type = "MS"))
})

test_that("anova.fit11.e.works", {
  library(RRPP)
  data("Pupfish")
  succeed(anova(lm.rrpp(coords ~ Pop * Sex, data = Pupfish, 
                        int.first = TRUE,
                        print.progress = FALSE, iter = 3, verbose = TRUE),
                error = c("Pop:Sex", "Pop:Sex", "Residuals")))
})


test_that("anova.fit12.works", {
  library(RRPP)
  data("PlethMorph")
  succeed(anova(lm.rrpp(TailLength ~ SVL, data = PlethMorph, 
                  print.progress = FALSE, iter = 3, verbose = TRUE)))
})


### predict.lm.rrpp -------------------------------------------------------


test_that("predict1.works", {
  library(RRPP)
  data("PupfishHeads")
  suppressWarnings(fit <- lm.rrpp(log(headSize) ~ sex + locality/year, SS.type = "III", 
                 data = PupfishHeads, print.progress = FALSE, iter = 3, verbose = TRUE))
  succeed(predict(fit))
})


### manova.update ---------------------------------------------------------

test_that("manova1.works", {
  library(RRPP)
  data("Pupfish")
  fit <- lm.rrpp(coords ~ CS + Sex, SS.type = "I", 
          data = Pupfish, print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(manova.update(fit, print.progress = FALSE))
})

### pairwise --------------------------------------------------------------

test_that("pairwise1.works", {
  library(RRPP)
  data("Pupfish")
  fit <- lm.rrpp(coords ~ CS + Sex, SS.type = "I", 
          data = Pupfish, print.progress = FALSE, iter = 3, verbose = TRUE)
  succeed(pairwise(fit, groups = Pupfish$Sex, print.progress = FALSE))
})


### trajectory.analysis ---------------------------------------------------

test_that("trajectory.analysis1.works", {
  library(RRPP)
  data("Pupfish")
  fit <- lm.rrpp(coords ~ Pop * Sex, data = Pupfish, print.progress = FALSE,
                 iter = 3, verbose = TRUE)
  succeed(trajectory.analysis(fit, groups = Pupfish$Pop, 
                              traj.pts = Pupfish$Sex, print.progress = FALSE))
})

test_that("trajectory.analysis2.works", {
  library(RRPP)
  data("motionpaths")
  fit <- lm.rrpp(trajectories ~ groups, data = motionpaths, iter = 3, verbose = TRUE)
  succeed(trajectory.analysis(fit, groups = motionpaths$groups, traj.pts = 5))
})

### ordinate --------------------------------------------------------------


test_that("ordinate1.works", {
  library(RRPP)
  data("PlethMorph")
  Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
                                  "Snout.eye", "BodyWidth", 
                                  "Forelimb", "Hindlimb")])
  PlethMorph$Y <- as.matrix(Y)
  R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
               iter = 0, print.progress = FALSE)$LM$residuals
  succeed(ordinate(R, scale. = FALSE))
})

test_that("ordinate2.works", {
  library(RRPP)
  data("PlethMorph")
  Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
                                  "Snout.eye", "BodyWidth", 
                                  "Forelimb", "Hindlimb")])
  PlethMorph$Y <- as.matrix(Y)
  R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
               iter = 0, print.progress = FALSE)$LM$residuals
  succeed(ordinate(R, scale. = TRUE))
})

test_that("ordinate3.works", {
  library(RRPP)
  data("PlethMorph")
  Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
                                  "Snout.eye", "BodyWidth", 
                                  "Forelimb", "Hindlimb")])
  PlethMorph$Y <- as.matrix(Y)
  R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
               iter = 0, print.progress = FALSE)$LM$residuals
  succeed(ordinate(R, scale. = TRUE, 
                   transform. = FALSE, 
                   Cov = PlethMorph$PhyCov))
})

test_that("ordinate4.works", {
  library(RRPP)
  data("PlethMorph")
  Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
                                  "Snout.eye", "BodyWidth", 
                                  "Forelimb", "Hindlimb")])
  PlethMorph$Y <- as.matrix(Y)
  R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
               iter = 0, print.progress = FALSE)$LM$residuals
  succeed(ordinate(R, scale. = TRUE, 
                   transform. = TRUE, 
                   Cov = PlethMorph$PhyCov))
})

test_that("ordinate5.works", {
  library(RRPP)
  data("PlethMorph")
  Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
                                  "Snout.eye", "BodyWidth", 
                                  "Forelimb", "Hindlimb")])
  PlethMorph$Y<- as.matrix(Y)
  R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
               iter = 0, print.progress = FALSE)$LM$residuals
  succeed(ordinate(R, A = PlethMorph$PhyCov, scale. = TRUE))
})

test_that("ordinate6.works", {
  library(RRPP)
  data("PlethMorph")
  Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
                                  "Snout.eye", "BodyWidth", 
                                  "Forelimb", "Hindlimb")])
  PlethMorph$Y <- as.matrix(Y)
  R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
               iter = 0, print.progress = FALSE)$LM$residuals
  succeed(ordinate(R, A = PlethMorph$PhyCov, 
                   scale. = TRUE,
                   transform. = FALSE, 
                   Cov = PlethMorph$PhyCov))
})

### measurement.error ----------------------------------------------------


test_that("measurement.error1.works", {
  library(RRPP)
  data("fishy")
  succeed(ME1 <- measurement.error(Y = "coords", subjects = "subj",
                           replicates = "reps", data = fishy, 
                           iter = 3, Parallel = FALSE,
                           verbose = TRUE))
})

test_that("measurement.error2.works", {
  library(RRPP)
  data("fishy")
  succeed(ME2 <- measurement.error(Y = "coords", subjects = "subj",
                                   replicates = "reps", data = fishy,
                                   groups = "groups",
                                   iter = 3, Parallel = FALSE,
                                   verbose = TRUE))
})

test_that("measurement.error3.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                                   replicates = "reps", data = fishy,
                                   groups = "groups",
                                   iter = 3, Parallel = FALSE,
                                   verbose = TRUE)
  succeed(anova(ME2))
})

test_that("ICCstats1.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                                   replicates = "reps", data = fishy,
                                   groups = "groups",
                                   iter = 3, Parallel = FALSE,
                           verbose = TRUE)
  succeed(ICCstats(ME2, subjects = "Subjects"))
})

test_that("ICCstats2.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                           replicates = "reps", data = fishy,
                           groups = "groups",
                           iter =3, Parallel = FALSE,
                           verbose = TRUE)
  succeed(ICCstats(ME2, subjects = "Subjects", 
                   with_in = "Systematic ME"))
})

test_that("ICCstats3.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                           replicates = "reps", data = fishy,
                           groups = "groups",
                           iter = 3, Parallel = FALSE,
                           verbose = TRUE)
  succeed(ICCstats(ME2, subjects = "Subjects", 
                   with_in = c("Systematic ME", "Systematic ME:Groups")))
})

### getModels---------------------------------------------------------------

test_that("getModels.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                           replicates = "reps", data = fishy,
                           groups = "groups",
                           iter = 3, Parallel = FALSE,
                           verbose = TRUE)
  succeed(getModels(ME2, "all"))
})

### getANOVAstats ------------------------------------------------------------

test_that("getANOVAStats.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                           replicates = "reps", data = fishy,
                           groups = "groups",
                           iter = 3, Parallel = FALSE,
                           verbose = TRUE)
  succeed(getANOVAStats(ME2, "all"))
})

### getPermInfo ----------------------------------------------------------

test_that("getANOVAStats.works", {
  library(RRPP)
  data("fishy")
  ME2 <- measurement.error(Y = "coords", subjects = "subj",
                           replicates = "reps", data = fishy,
                           groups = "groups",
                           iter = 3, Parallel = FALSE,
                           verbose = TRUE)
  succeed(getPermInfo(ME2, "all"))
})


### getModelCov---------------------------------------------------------------

test_that("getModelCov1.works", {
  library(RRPP)
  data("PlethMorph")
  fitGLS <- lm.rrpp(TailLength ~ SVL, 
                    data = PlethMorph, 
                    Cov = PlethMorph$PhyCov, 
                    print.progress = FALSE, iter = 3)
  
  succeed(getModelCov(fitGLS, "Cov"))
})

test_that("getModelCov2.works", {
  library(RRPP)
  data("PlethMorph")
  fitGLS <- lm.rrpp(TailLength ~ SVL, 
                    data = PlethMorph, 
                    Cov = PlethMorph$PhyCov, 
                    print.progress = FALSE, iter = 3)
  
  succeed(getModelCov(fitGLS, "Pcov"))
})


### getResCov ---------------------------------------------------------------

test_that("getResCov1.works", {
  library(RRPP)
  data("PlethMorph")
  fitGLS <- lm.rrpp(TailLength ~ SVL, 
                    data = PlethMorph, 
                    Cov = PlethMorph$PhyCov, 
                    print.progress = FALSE, iter = 3)
  
  succeed(getResCov(fitGLS, useDf = TRUE, standardize = FALSE))
})

test_that("getResCov2.works", {
  library(RRPP)
  data("PlethMorph")
  fitGLS <- lm.rrpp(TailLength ~ SVL, 
                    data = PlethMorph, 
                    Cov = PlethMorph$PhyCov, 
                    print.progress = FALSE, iter = 3)
  
  succeed(getResCov(fitGLS, useDf = TRUE, standardize = TRUE))
})

test_that("getResCov3.works", {
  library(RRPP)
  data("PlethMorph")
  fitGLS <- lm.rrpp(TailLength ~ SVL, 
                    data = PlethMorph, 
                    Cov = PlethMorph$PhyCov, 
                    print.progress = FALSE, iter = 3)
  
  succeed(getResCov(fitGLS, useDf = FALSE, standardize = FALSE))
})

test_that("getResCov4.works", {
  library(RRPP)
  data("PlethMorph")
  fitGLS <- lm.rrpp(TailLength ~ SVL, 
                    data = PlethMorph, 
                    Cov = PlethMorph$PhyCov, 
                    print.progress = FALSE, iter = 3)
  
  succeed(getResCov(fitGLS, useDf = FALSE, standardize = TRUE))
})