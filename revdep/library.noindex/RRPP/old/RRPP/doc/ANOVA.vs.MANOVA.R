## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(RRPP)
data(Pupfish)
fit <- lm.rrpp(coords ~ Sex*Pop, SS.type = "I", 
               data = Pupfish, print.progress = FALSE) 
attributes(fit)
attributes(fit$ANOVA)
attributes(getANOVAStats(fit, stat = "all"))

## -----------------------------------------------------------------------------
anova(fit)

## -----------------------------------------------------------------------------
fitm <- manova.update(fit, print.progress = FALSE, tol = 0)
attributes(fitm)
attributes(fitm$MANOVA)

## -----------------------------------------------------------------------------
summary(fitm, test = "Roy")
summary(fitm, test = "Pillai")


## -----------------------------------------------------------------------------
fitm <- manova.update(fit, print.progress = FALSE, tol = 0.001)
anova(fit)
summary(fitm)

## -----------------------------------------------------------------------------
fitm <- manova.update(fit, print.progress = FALSE, PC.no  = 10)
anova(fit)
summary(fitm)

