pkgname <- "geomorph"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('geomorph')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("arrayspecs")
### * arrayspecs

flush(stderr()); flush(stdout())

### Name: arrayspecs
### Title: Convert landmark data matrix into array (p x k x n)
### Aliases: arrayspecs
### Keywords: utilities

### ** Examples

x<-matrix(rnorm(18),nrow=3)  # Random triangles (all coordinates on same 
 # row for each triangle)
arrayspecs(x,3,2) 
 
x2<-matrix(rnorm(18),ncol=2) # Random triangles (each landmark on its 
# own row)
arrayspecs(x2,3,2)



cleanEx()
nameEx("bilat.symmetry")
### * bilat.symmetry

flush(stderr()); flush(stdout())

### Name: bilat.symmetry
### Title: Analysis of bilateral symmetry
### Aliases: bilat.symmetry
### Keywords: analysis

### ** Examples

#Example of matching symmetry
# NOT RUN
# data(mosquito)
# gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, 
# ind=mosquito$ind, 
# side=mosquito$side,
# replicate=mosquito$replicate)
# mosquito.sym <- bilat.symmetry(A = wingshape, ind = ind, side = side,
# replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 149, 
# data = gdf)
# summary(mosquito.sym)
# plot(mosquito.sym, warpgrids = TRUE)
# mosquito.sym$shape.anova # extract just the anova table on shape

# Previous example, performing GPA first
# Y.gpa <- gpagen(mosquito$wingshape)
# mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = ind, side = side,
# replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 149, 
# data = gdf)
# summary(mosquito.sym2)
# summary(mosquito.sym) # same results

#Example of object symmetry

# data(lizards)
# gdf <- geomorph.data.frame(shape = lizards$coords, 
# ind = lizards$ind, 
# replicate = lizards$rep)
# liz.sym <- bilat.symmetry(A = shape, ind = ind, rep = rep, 
# object.sym = TRUE, 
# land.pairs = lizards$lm.pairs, data = gdf, RRPP = TRUE, iter = 149)
# summary(liz.sym)

# Example of object symmetry in 3D and including semilandmarks

# data(scallops)
# gdf <- geomorph.data.frame(shape = scallops$coorddata, 
# ind = scallops$ind)
# scallop.sym <- bilat.symmetry(A = shape, ind = ind, 
# object.sym = TRUE, 
# curves= scallops$curvslide, surfaces = scallops$surfslide,
# land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 149)
# summary(scallop.sym)
# NOTE one can also: plot(scallop.sym, warpgrids = TRUE, mesh = NULL)
# NOTE one can also: scallop.sym$data.type # recall the symmetry type



cleanEx()
nameEx("combine.subsets")
### * combine.subsets

flush(stderr()); flush(stdout())

### Name: combine.subsets
### Title: Combine separate landmark configurations
### Aliases: combine.subsets
### Keywords: utilities

### ** Examples

# NOT RUN
# data(larvalMorph) 
# head.gpa <- gpagen(larvalMorph$headcoords, 
#   curves = larvalMorph$head.sliders)
# tail.gpa <- gpagen(larvalMorph$tailcoords, 
 # curves = larvalMorph$tail.sliders)

# Combine original data without GPA (plot to see relative size of  
# heads and tails)

 # all.lm <- combine.subsets(head = larvalMorph$headcoords,
 # tail = larvalMorph$tailcoords, gpa = FALSE, CS.sets = NULL)
 # plotAllSpecimens((all.lm$coords))
 
 # Combine with GPA and relative centroid size
 
# comb.lm <- combine.subsets(head = head.gpa, tail = tail.gpa, gpa = TRUE)
# summary(comb.lm)

# (configurations are actual relative size)
# comb.lm$coords[,,1]

# Plot all specimens and just first specimen and color code landmarks 
# par(mfrow = c(1,2))
# plotAllSpecimens(comb.lm$coords)
# plot(comb.lm$coords[,,1], pch = 21, bg = c(rep(1,26), 
# rep(2,64)), asp = 1)

# Override relative centroid size

# comb.lm <- combine.subsets(head = head.gpa$coords, 
# tail = tail.gpa$coords, gpa = FALSE, CS.sets = NULL)
# par(mfrow = c(1,2))
# plotAllSpecimens(comb.lm$coords)
# plot(comb.lm$coords[,,1], pch = 21, bg = c(rep(1,26), 
# rep(2,64)), asp = 1)

# Note the head is as large as the tail, which is quite unnatural.

## Normalizing centroid size

# comb.lm <- combine.subsets(head = head.gpa, 
# tail = tail.gpa, gpa = TRUE, norm.CS = TRUE)
# summary(comb.lm)
# par(mfrow = c(1,2))
# plotAllSpecimens(comb.lm$coords)
# plot(comb.lm$coords[,,1], pch = 21, bg = c(rep(1,26), 
# rep(2,64)), asp = 1)
# par(mfrow = c(1,1))

# Note that the head is too large, compared to a real specimen.  
# This option focuses on average distance of points to centroid, 
# but ignores the number of landmarks.  
# Consequently,the density of landmarks in the head and tail are 
# irrelevant and the head size is inflated because of the fewer 
# landmarks in the configuration.

## Weighting centroid size

# comb.lm <- combine.subsets(head = head.gpa, 
# tail = tail.gpa, gpa = TRUE, norm.CS = FALSE, weights = c(0.3, 0.7))
# summary(comb.lm)
# par(mfrow = c(1,2))
# plotAllSpecimens(comb.lm$coords)
# plot(comb.lm$coords[,,1], pch = 21, bg = c(rep(1,26), 
# rep(2,64)), asp = 1)
# par(mfrow = c(1,1))

# Note that the head is way too small, compared to a real specimen.  
# This option allows one to dictate the relative sizes of subsets
#  as portions of the combined set.  An option like this should be 
# used with caution, but can help overcome issues caused by landmark 
# density.



cleanEx()
nameEx("compare.CR")
### * compare.CR

flush(stderr()); flush(stdout())

### Name: compare.CR
### Title: Comparisons of Effect Sizes from Modularity Analyses
### Aliases: compare.CR
### Keywords: analysis

### ** Examples


#NOT RUN
# Example 1: Compare modular signal across datasets
 
# data(pupfish) 
# Y.gpa<-gpagen(pupfish$coords, print.progress = FALSE)    #GPA-alignment   
 
## landmarks on the body and operculum
# land.gps<-rep('a',56); land.gps[39:48]<-'b'

# group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
# levels(group)

# coords.gp <- coords.subset(Y.gpa$coords, group)

# modul.tests <- Map(function(x) modularity.test(x, land.gps,iter=999, 
# print.progress = FALSE), coords.gp) 
         
# the map function performs the integration test on each 3D array 
# in the lists provided

#  modul.tests$Marsh.F
#  modul.tests$Marsh.M
#  modul.tests$Sinkhole.F
#  modul.tests$Sinkhole.M

# group.Z <- compare.CR(modul.tests, CR.null = FALSE)
# summary(group.Z)

# Example 2: Compare alternative modular hypotheses

# 3 module hypothesis (tail now a module)
# land.gps3 <- rep('a',56); land.gps3[39:48]<-'b'; 
# land.gps3[c(6:9,28:38)] <- 'c' 
   
# 4 module hypothesis (eye now a module)
# land.gps4 <- rep('a',56); land.gps4[39:48]<-'b'; 
# land.gps4[c(6:9,28:38)] <- 'c'; 
#  land.gps4[c(10,49:56)] <- 'd'  

# m3.test <- modularity.test(coords.gp$Marsh.F,land.gps3, iter = 499, 
# print.progress = FALSE)
# m4.test <- modularity.test(coords.gp$Marsh.F,land.gps4, iter = 499, 
# print.progress = FALSE)

# model.Z <- compare.CR(modul.tests$Marsh.F,m3.test,m4.test, 
# CR.null = TRUE)
# summary(model.Z)




cleanEx()
nameEx("compare.ZVrel")
### * compare.ZVrel

flush(stderr()); flush(stdout())

### Name: compare.ZVrel
### Title: Comparisons of Effect Sizes from Overall Integration Analyses
### Aliases: compare.ZVrel
### Keywords: analysis

### ** Examples

 
 data("plethodon")
 Y.gpa <- gpagen(plethodon$land)
 
 coords.gp <- coords.subset(Y.gpa$coords, plethodon$species)
 Vrel.gp <- Map(function(x) integration.Vrel(x), coords.gp) 
 
 out <- compare.ZVrel(Vrel.gp$Jord, Vrel.gp$Teyah)
 
 summary(out)
 



cleanEx()
nameEx("compare.evol.rates")
### * compare.evol.rates

flush(stderr()); flush(stdout())

### Name: compare.evol.rates
### Title: Comparing net rates of shape evolution on phylogenies
### Aliases: compare.evol.rates
### Keywords: analysis

### ** Examples

data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
 gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
 names(gp.end)<-plethspecies$phy$tip

ER<-compare.evol.rates(A=Y.gpa$coords, phy=plethspecies$phy,
  method="simulation",gp=gp.end,iter=999)
summary(ER)
plot(ER)



cleanEx()
nameEx("compare.multi.evol.rates")
### * compare.multi.evol.rates

flush(stderr()); flush(stdout())

### Name: compare.multi.evol.rates
### Title: Comparing net rates of evolution among traits on phylogenies
### Aliases: compare.multi.evol.rates
### Keywords: analysis

### ** Examples


data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
land.gp<-c("A","A","A","A","A","B","B","B","B","B","B")  
    #mandible and cranium subsets

EMR<-compare.multi.evol.rates(A=Y.gpa$coords,gp=land.gp, 
    Subset=TRUE, phy= plethspecies$phy,iter=999)
summary(EMR)
plot(EMR)



cleanEx()
nameEx("compare.physignal.z")
### * compare.physignal.z

flush(stderr()); flush(stdout())

### Name: compare.physignal.z
### Title: Comparisons of Phylogenetic Signal Effect Sizes
### Aliases: compare.physignal.z
### Keywords: analysis

### ** Examples


# Example: Compare phylogenetic signal of head components in Plethodon

data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment

## landmarks of the jaw and cranium
jaw <- 1:5
cranium <- 6:11

PS.jaw <- physignal.z(A = Y.gpa$coords[jaw,,], phy = plethspecies$phy, 
lambda = "front", PAC.no = 7, iter=999)

PS.cranium <- physignal.z(A = Y.gpa$coords[cranium,,], phy = plethspecies$phy, 
lambda = "front", PAC.no = 7, iter=999)

PS.list <-list(PS.jaw, PS.cranium)
names(PS.list) <- c("jaw", "cranium")

PS.Z <- compare.physignal.z(PS.list)
summary(PS.Z)




cleanEx()
nameEx("compare.pls")
### * compare.pls

flush(stderr()); flush(stdout())

### Name: compare.pls
### Title: Comparisons of Effect Sizes from Partial Least Squares
### Aliases: compare.pls
### Keywords: analysis

### ** Examples

# Example of comparative morphological integration between pupfish head 
# and body shapes

 data(pupfish) # GPA previously performed
  
 group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
 levels(group)
  
 tail.LM <- c(1:3, 5:9, 18:38)
 head.LM <- (1:56)[-tail.LM]

 tail.coords <- pupfish$coords[tail.LM,,]
 head.coords <- pupfish$coords[head.LM,,]
 
 # Subset 3D array by group, returning a list of 3D arrays
 tail.coords.gp <- coords.subset(tail.coords, group)
 head.coords.gp <- coords.subset(head.coords, group)

 integ.tests <- Map(function(x,y) integration.test(x, y, iter=499, 
 print.progress = FALSE), head.coords.gp, tail.coords.gp)
# the map function performs the integration test on each 3D array in 
# the lists provided

 integ.tests$Marsh.F
 integ.tests$Marsh.M
 integ.tests$Sinkhole.F
 integ.tests$Sinkhole.M

 group.Z <- compare.pls(integ.tests)
 summary(group.Z)

 # Sexual dimorphism in morphological integration in one population
 # but not the other

 # can also list different PLS analyses, separately

 compare.pls(MF = integ.tests$Marsh.F, MM = integ.tests$Marsh.M)




cleanEx()
nameEx("coords.subset")
### * coords.subset

flush(stderr()); flush(stdout())

### Name: coords.subset
### Title: Subset landmark coordinates via a factor
### Aliases: coords.subset
### Keywords: utilities

### ** Examples

data(pupfish) 
group <- factor(paste(pupfish$Pop, pupfish$Sex))
levels(group)
new.coords <- coords.subset(A = pupfish$coords, group = group)
names(new.coords) # see the list levels
# group shape means
lapply(new.coords, mshape)




cleanEx()
nameEx("define.sliders")
### * define.sliders

flush(stderr()); flush(stdout())

### Name: define.sliders
### Title: Select points to "slide" along curves
### Aliases: define.sliders
### Keywords: utilities

### ** Examples

 
## (not run) Use interactive function in rgl window
 # data(scallops)
 # define.sliders(scallops$coorddata[,,1], nsliders=11,
 #   surfsliders = scallops$surfslide) 
 # here the first specimen is used for plotting purposes only
 
## Examples of AUTO mode 
 ## 1 curve of sliding semilandmark
 # Define sliders for scallopdata
 #sliders = define.sliders(c(5:16,1))

 ## 2 curves of sliding semilandmarks
 # Define sliders for 10 landmarks, where LMs 1, 5, and 10 fixed
 # 2, 3, and 4 are along a curve between 1 and 5
 # and 6, 7, 8, and 9 are along a curve between 5 and 10.
 #sliders = rbind(define.sliders(1:5), define.sliders(5:10)) 



cleanEx()
nameEx("estimate.missing")
### * estimate.missing

flush(stderr()); flush(stdout())

### Name: estimate.missing
### Title: Estimate locations of missing landmarks
### Aliases: estimate.missing
### Keywords: utilities

### ** Examples

data(plethodon)
plethland<-plethodon$land
  plethland[3,,2]<-plethland[8,,2]<-NA  #create missing landmarks
  plethland[3,,5]<-plethland[8,,5]<-plethland[9,,5]<-NA  
  plethland[3,,10]<-NA  
  
estimate.missing(plethland,method="TPS")
estimate.missing(plethland,method="Reg")



cleanEx()
nameEx("fixed.angle")
### * fixed.angle

flush(stderr()); flush(stdout())

### Name: fixed.angle
### Title: Rotate a subset of 2D landmarks to common articulation angle
### Aliases: fixed.angle
### Keywords: utilities

### ** Examples

#Example using Plethodon
#Articulation point is landmark 1, rotate mandibular landmarks (2-5) 
# relative to cranium

data(plethspecies) 
# Using specific points:
newLM1 <- fixed.angle(plethspecies$land,
art.pt=1, angle.pts.1 = 5, 
angle.pts.2 = 6, rot.pts = c(2,3,4,5))
Y.gpa1 <- gpagen(newLM1)
plot(Y.gpa1, mean = FALSE)

# Using centroids from subsets
newLM2 <- fixed.angle(plethspecies$land,art.pt=1, 
angle.pts.1 = c(1, 6:11), 
angle.pts.2 = 2:5, 
rot.pts = NULL, angle = 20, 
degrees = TRUE) # rotated points same as second partition
Y.gpa2 <- gpagen(newLM2)
plot(Y.gpa2, mean = FALSE)




cleanEx()
nameEx("geomorph.data.frame")
### * geomorph.data.frame

flush(stderr()); flush(stdout())

### Name: geomorph.data.frame
### Title: Create a data frame with shape data
### Aliases: geomorph.data.frame
### Keywords: utilities

### ** Examples

data(plethodon) 
Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)
gdf <- geomorph.data.frame(Y.gpa)
attributes(gdf)

gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
site = plethodon$site)
attributes(gdf)

# Using geomorph.data.frame to facilitate analysis
anova(procD.lm(coords ~ Csize + species * site, data = gdf))



cleanEx()
nameEx("globalIntegration")
### * globalIntegration

flush(stderr()); flush(stdout())

### Name: globalIntegration
### Title: Quantify global integration relative to self-similarity
### Aliases: globalIntegration
### Keywords: analysis

### ** Examples

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    

globalIntegration(Y.gpa$coords)



cleanEx()
nameEx("gm.prcomp")
### * gm.prcomp

flush(stderr()); flush(stdout())

### Name: gm.prcomp
### Title: Principal and phylogenetically-aligned components analysis of
###   shape data
### Aliases: gm.prcomp
### Keywords: visualization

### ** Examples

 data(plethspecies) 
 Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
 
 ###  Traditional PCA 
 PCA <- gm.prcomp(Y.gpa$coords)
 summary(PCA)
 plot(PCA, main = "PCA")
 plot(PCA, main = "PCA", flip = 1) # flip the first axis
 plot(PCA, main = "PCA", axis1 = 3, axis2 = 4) # change PCs viewed
 
 ### Phylomorphospace - PCA with phylogeny (result is same as above, 
 ### but with estimated ancestral states projected into plot)
 PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
 summary(PCA.w.phylo)
 plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo")
 
 ### Phylogenetic PCA - PCA based on GLS-centering and projection
 # This is the same as the method described by Revell (2009)
 phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE)
 summary(phylo.PCA)
 plot(phylo.PCA, phylo = TRUE, main = "phylo PCA")
 
 ### Phylogenetic PCA - PCA based on GLS-centering and transformed 
 # projection
 # This produces a PCA independent of phylogeny
 phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
 GLS = TRUE, transform = TRUE)
 summary(phylo.tPCA)
 plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA")
 
 ### PaCA - Alignment of data to physlogenetic signal rather than axis of 
 ### greatest variation, like in PCA
 
 # OLS method (rotation of PCA)
 PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
   align.to.phy = TRUE)
 summary(PaCA.ols)
 plot(PaCA.ols, phylo = TRUE, main = "PaCA using OLS")
 
 # GLS method (rotation of Phylogenetic PCA)
 PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
 align.to.phy = TRUE, GLS = TRUE)
 summary(PaCA.gls)
 plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS")
 
 # GLS method (rotation of Phylogenetic PCA with transformed data)
 PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
 align.to.phy = TRUE, GLS = TRUE, transform = TRUE)
 summary(PaCA.gls)
 plot(PaCA.gls, phylo = TRUE, 
   main = "PaCA using GLS and transformed projection")
 
 
 ### Advanced Plotting
 gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
 par(mar=c(2, 2, 2, 2))
 plot(PaCA.ols, pch=22, cex = 1.5, bg = gps, phylo = TRUE) 
 # Modify options as desired
 #  Add things as desired using standard R plotting
 text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", 
   pos = 4, font = 2)
 text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
 legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))
 
 ### NOT RUN
 ### 3D plot with a phylogeny and time on the z-axis
 # plot(PCA.w.phylo, time.plot = TRUE)
 # plot(PCA.w.phylo, time.plot = TRUE, bg = "red", 
 #    phylo.par = list(tip.labels = TRUE, 
 # tip.txt.cex = 2, edge.color = "blue", edge.width = 2))
 



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("gpagen")
### * gpagen

flush(stderr()); flush(stdout())

### Name: gpagen
### Title: Generalized Procrustes analysis of points, curves, and surfaces
### Aliases: gpagen
### Keywords: analysis

### ** Examples

# Example 1: fixed points only
data(plethodon) 
Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE)
summary(Y.gpa)
plot(Y.gpa)

# Example 2: points and semilandmarks on curves
data(hummingbirds)

###Slider matrix
hummingbirds$curvepts

# Using bending energy for sliding
Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
ProcD = FALSE)   
summary(Y.gpa)
plot(Y.gpa)


# Using Procrustes Distance for sliding
Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
ProcD = TRUE)   
summary(Y.gpa)
plot(Y.gpa)


# Example 3: points, curves and surfaces
data(scallops)

# Using Procrustes Distance for sliding
Y.gpa <- gpagen(A = scallops$coorddata, curves = scallops$curvslide, 
surfaces = scallops$surfslide)
# NOTE can summarize as: summary(Y.gpa)
# NOTE can plot as: plot(Y.gpa) 



cleanEx()
nameEx("gridPar")
### * gridPar

flush(stderr()); flush(stdout())

### Name: gridPar
### Title: Set up parameters for grids, points, and links in
###   plotRefToTarget
### Aliases: gridPar
### Keywords: utilities visualization

### ** Examples

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
ref<-mshape(Y.gpa$coords)
plotRefToTarget(ref,Y.gpa$coords[,,39]) # default settings

# Altering points and links
GP1 <- gridPar(pt.bg = "red", pt.size = 1, link.col="blue", link.lwd=2, 
n.col.cell=50)
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP1, mag=2, 
links=plethodon$links, method="TPS")

# Altering point color
GP2 <- gridPar(pt.bg = "green", pt.size = 1) 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP2, mag=3, 
method="vector")

# Altering ref and target points
GP3 <- gridPar(pt.bg = "blue", pt.size = 1.5, tar.pt.bg = "orange", 
tar.pt.size = 1) 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP3, mag=3, 
method="points")

# Altering outline color
GP4 <- gridPar(tar.out.col = "red", tar.out.cex = 0.3) 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP4, mag=3, 
outline=plethodon$outline, method="TPS")

# Altering text labels
GP5 <- gridPar(txt.pos = 3, txt.col = "red") 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP5, mag=3, 
method="vector", label=TRUE)



cleanEx()
nameEx("integration.Vrel")
### * integration.Vrel

flush(stderr()); flush(stdout())

### Name: integration.Vrel
### Title: Quantify integration in a set of traits
### Aliases: integration.Vrel
### Keywords: analysis

### ** Examples

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
integration.Vrel(Y.gpa$coords)



cleanEx()
nameEx("integration.test")
### * integration.test

flush(stderr()); flush(stdout())

### Name: integration.test
### Title: Quantify morphological integration between modules
### Aliases: integration.test
### Keywords: analysis

### ** Examples

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#landmarks on the skull and mandible assigned to partitions
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
IT <- integration.test(Y.gpa$coords, partition.gp=land.gps, iter=999)
summary(IT) # Test summary
P <- plot(IT) # PLS plot

 # Visualize shape at minimum and maximum PLS scores.
 # This is the challenging way
 
 # Block 1 
 minx <- min(P$plot.args$x)
 maxx <- max(P$plot.args$x)
 preds <- shape.predictor(P$A1, 
 x = P$plot.args$x,
 min = minx, max = maxx)
 plotRefToTarget(mshape(P$A1), preds$min)
 plotRefToTarget(mshape(P$A1), preds$max)

 # Block 2 
 miny <- min(P$plot.args$y)
 maxy <- max(P$plot.args$y)
 preds <- shape.predictor(P$A2, 
 x = P$plot.args$y,
 min = miny, max = maxy)
 plotRefToTarget(mshape(P$A2), preds$min)
 plotRefToTarget(mshape(P$A2), preds$max)
 
 ### Visualize shape variation using picknplot.shape Because picknplot  
 ### requires user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(P)

IT$left.pls.vectors # extracting just the left (first block) 
# singular vectors



cleanEx()
nameEx("interlmkdist")
### * interlmkdist

flush(stderr()); flush(stdout())

### Name: interlmkdist
### Title: Calculate linear distances between landmarks
### Aliases: interlmkdist
### Keywords: utilities

### ** Examples

 
data(plethodon)
# Make a matrix defining three interlandmark distances 
lmks <- matrix(c(8,9,6,12,4,2), ncol=2, byrow=TRUE, 
dimnames = list(c("eyeW", "headL", "mouthL"),c("start", "end")))
# where 8-9 is eye width; 6-12 is head length; 4-2 is mouth length
# or alternatively
lmks <- data.frame(eyeW = c(8,9), headL = c(6,12), mouthL = c(4,2), 
row.names = c("start", "end")) 
A <- plethodon$land
lineardists <- interlmkdist(A, lmks)



cleanEx()
nameEx("make_ggplot")
### * make_ggplot

flush(stderr()); flush(stdout())

### Name: make_ggplot
### Title: Convert geomorph plots to ggplot objects
### Aliases: make_ggplot
### Keywords: utilities

### ** Examples


### PLS Example
# NOT RUN
# data(plethodon) 
# Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
# landmarks on the skull and mandible assigned to partitions
# land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
# IT <- integration.test(Y.gpa$coords, partition.gp=land.gps, iter=999)
# summary(IT) # Test summary
# P <- plot(IT) # PLS plot
# make_ggplot(P) # same plot in ggplot

### Allometry example

# data(plethodon) 
# Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment  

# gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
#                           species = plethodon$species) 

# fit <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=0, 
#                 print.progress = FALSE)

# P <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
#                     pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))

# make_ggplot(P)

### Tangent Space plot

# data(plethspecies) 
# Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment

# PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
# P <- plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo")
# make_ggplot(P)




cleanEx()
nameEx("modularity.test")
### * modularity.test

flush(stderr()); flush(stdout())

### Name: modularity.test
### Title: Evaluate the degree of modular signal in shape data
### Aliases: modularity.test
### Keywords: analysis

### ** Examples

# Not Run
# data(pupfish) 
# Y.gpa<-gpagen(pupfish$coords, print.progress = FALSE)    #GPA-alignment    
 #landmarks on the body and operculum
# land.gps<-rep('a',56); land.gps[39:48]<-'b'

# MT <- modularity.test(Y.gpa$coords,land.gps,CI=FALSE,iter=99)
# summary(MT) # Test summary
# plot(MT) # Histogram of CR sampling distribution 
# Result implies modularity present



cleanEx()
nameEx("morphol.disparity")
### * morphol.disparity

flush(stderr()); flush(stdout())

### Name: morphol.disparity
### Title: Morphological disparity for one or more groups of specimens
### Aliases: morphol.disparity
### Keywords: analysis

### ** Examples

# Not Run
# data(plethodon)
# Y.gpa<-gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment
# gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
# site = plethodon$site)

# Morphological disparity for entire data set
# morphol.disparity(coords ~ 1, groups = NULL, data = gdf, 
# iter = 999, print.progress = FALSE)

# Morphological disparity for entire data set, accounting for allometry
# morphol.disparity(coords ~ Csize, groups= NULL, data = gdf, 
# iter = 999, print.progress = FALSE)

# Morphological disparity without covariates, using overall mean
# morphol.disparity(coords ~ 1, groups= ~ species*site, data = gdf, 
# iter = 999, print.progress = FALSE)

# Morphological partial disparities for overal mean
# morphol.disparity(coords ~ 1, groups= ~ species*site, partial = TRUE, 
# data = gdf, iter = 999, print.progress = FALSE)

# Morphological disparity without covariates, using group means
# morphol.disparity(coords ~ species*site, groups= ~species*site, 
# data = gdf, iter = 999, print.progress = FALSE)

# Morphological disparity of different groups than those 
# described by the linear model
# morphol.disparity(coords ~ Csize + species*site, groups= ~ species, 
# data = gdf, iter = 999, print.progress = FALSE)

# Extracting components
# MD <- morphol.disparity(coords ~ Csize + species*site, groups= ~ species, 
# data = gdf, iter = 999, print.progress = FALSE)
# MD$Procrustes.var # just the Procrustes variances


### Morphol.disparity can be used with previously-defined 
### procD.lm or lm.rrpp class objects

# data(plethspecies)
# Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
# gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
# names(gp.end)<-plethspecies$phy$tip

# gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy, 
# gp.end = gp.end)

# pleth.ols <- procD.lm(coords ~ Csize + gp.end, 
# data = gdf, iter = 999) # ordinary least squares
# pleth.pgls <- procD.pgls(coords ~ Csize + gp.end, phy = phy, 
# data = gdf, iter = 999) # phylogenetic generalized least squares

# summary(pleth.ols)
# summary(pleth.pgls)

# morphol.disparity(f1 = pleth.ols, groups = ~ gp.end, data = gdf, 
# iter = 999, print.progress = FALSE)
# morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end, 
# transform = FALSE, data = gdf, 
# iter = 999, print.progress = FALSE) # disparity in tangent space
# morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end,
# transform = TRUE, data = gdf, 
# iter = 999, print.progress = FALSE) # disparity in transformed space 

# Three plots that correspond to the three tests
# PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
# pPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
# GLS = TRUE, transform = FALSE)
# tpPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
# GLS = TRUE, transform = TRUE)

# par(mfrow = c(1,3))

# Phylomorphospace
# PC.plot <- plot(PCA, pch = 19, phylo = TRUE, main = "PCA-OLS")
# shapeHulls(PC.plot, groups = gp.end)

# Phylo-PCA
# pPC.plot <- plot(pPCA, pch = 19, phylo = TRUE, main = "pPCA - GLS, not transformed")
# shapeHulls(pPC.plot, groups = gp.end)

# Transformed phylo-PCA
# tpPC.plot <- plot(tpPCA, pch = 19, phylo = TRUE, main = "tpPCA - GLS, transformed")
# shapeHulls(tpPC.plot, groups = gp.end)

# par(mfrow = c(1,1))

 ### Variance using RRPP (not necessarily the same as disparity)
 # PW <- pairwise(pleth.ols, groups = gp.end)
 # summary(PW, test.type = 'var')
 # PW2 <- pairwise(pleth.pgls, groups = gp.end)
 # summary(PW2, test.type = 'var')




cleanEx()
nameEx("mshape")
### * mshape

flush(stderr()); flush(stdout())

### Name: mshape
### Title: Estimate mean shape for a set of aligned specimens
### Aliases: mshape
### Keywords: utilities

### ** Examples

data(plethodon)
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
A <- Y.gpa$coords
A[[1]] <- NA # make a missing value, just for example

mshape(Y.gpa$coords)   # mean (consensus) configuration
# mshape(A, na.action = 1) # will return an error
mshape(A, na.action = 2) # returns NA in spot of missing value
mshape(A, na.action = 3) # finds mean values from all possible values




cleanEx()
nameEx("na.omit.geomorph.data.frame")
### * na.omit.geomorph.data.frame

flush(stderr()); flush(stdout())

### Name: na.omit.geomorph.data.frame
### Title: Handle missing values in rrpp.data.frame objects
### Aliases: na.omit.geomorph.data.frame
### Keywords: utilities

### ** Examples

data(plethspecies)
Y.gpa <- gpagen(plethspecies$land, verbose = TRUE)
gdf <- geomorph.data.frame(Y.gpa)
gdf$d <- Y.gpa$procD
gdf$group <- c(rep(1, 4), rep(2, 4), NA) # one unknown group designation
gdf
ndf <- na.omit(gdf)
ndf



cleanEx()
nameEx("phylo.integration")
### * phylo.integration

flush(stderr()); flush(stdout())

### Name: phylo.integration
### Title: Quantify phylogenetic morphological integration between two or
###   more sets of variables under Brownian motion
### Aliases: phylo.integration
### Keywords: analysis

### ** Examples


data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 

IT<- phylo.integration(Y.gpa$coords,partition.gp=land.gps,
  phy=plethspecies$phy,iter=999)
summary(IT) # Test summary
P <- plot(IT) # PLS plot

 # Block 1 
 minx <- min(P$plot.args$x)
 maxx <- max(P$plot.args$x)
 preds <- shape.predictor(P$A1, 
 x = P$plot.args$x,
 min = minx, max = maxx)
 plotRefToTarget(mshape(P$A1), preds$min)
 plotRefToTarget(mshape(P$A1), preds$max)

 # Block 2 
 miny <- min(P$plot.args$y)
 maxy <- max(P$plot.args$y)
 preds <- shape.predictor(P$A2, 
 x = P$plot.args$y,
 min = miny, max = maxy)
 plotRefToTarget(mshape(P$A2), preds$min)
 plotRefToTarget(mshape(P$A2), preds$max)
 
 ### Visualize shape variation using picknplot.shape Because picknplot  
 ### requires user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(P)




cleanEx()
nameEx("phylo.modularity")
### * phylo.modularity

flush(stderr()); flush(stdout())

### Name: phylo.modularity
### Title: Evaluate the degree of phylogenetic modular signal in Procrustes
###   shape variables
### Aliases: phylo.modularity
### Keywords: analysis

### ** Examples

# Not Run
# data(plethspecies)
# Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
# land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 

# MT <- phylo.modularity(Y.gpa$coords, partition.gp=land.gps, 
# phy=plethspecies$phy, 
# CI = FALSE, iter=499)
# summary(MT) # Test summary
# plot(MT) # Histogram of CR sampling distribution 



cleanEx()
nameEx("physignal")
### * physignal

flush(stderr()); flush(stdout())

### Name: physignal
### Title: Assessing phylogenetic signal in Procrustes shape variables
### Aliases: physignal
### Keywords: analysis

### ** Examples

data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    

#Test for phylogenetic signal in shape
PS.shape <- physignal(A=Y.gpa$coords,phy=plethspecies$phy,iter=999)
summary(PS.shape)
plot(PS.shape)
plot(PS.shape$PACA, phylo = TRUE)
PS.shape$K.by.p # Phylogenetic signal profile

#Test for phylogenetic signal in size
PS.size <- physignal(A=Y.gpa$Csize,phy=plethspecies$phy,iter=999)
summary(PS.size)
plot(PS.size)



cleanEx()
nameEx("physignal.z")
### * physignal.z

flush(stderr()); flush(stdout())

### Name: physignal.z
### Title: Assessing phylogenetic signal effect size in Procrustes shape
###   variables
### Aliases: physignal.z
### Keywords: analysis

### ** Examples

data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    

# Test for phylogenetic signal in shape
PS.shape <- physignal.z(A = Y.gpa$coords, phy = plethspecies$phy, 
lambda = "front", iter=999)
summary(PS.shape)

# Problem with ill-conditioned residual covariance matrix; try shaving one dimension

PS.shape <- physignal.z(A = Y.gpa$coords, phy = plethspecies$phy, 
lambda = "front", PAC.no = 7, iter=999)
summary(PS.shape)
plot(PS.shape)
plot(PS.shape$PACA, phylo = TRUE)

# Test for phylogenetic signal in size
PS.size <- physignal.z(A = Y.gpa$Csize, phy = plethspecies$phy, 
lambda = "front", iter=999)
summary(PS.size)
plot(PS.size)



cleanEx()
nameEx("picknplot.shape")
### * picknplot.shape

flush(stderr()); flush(stdout())

### Name: picknplot.shape
### Title: Pick points in geomorph scatterplots to visualize shape
###   variation
### Aliases: picknplot.shape
### Keywords: visualization

### ** Examples


### Because picknplot requires user decisions, the following examples
### are not run (but can be with removal of #s)

# 2d
# data(plethodon) 
# Y.gpa <- gpagen(plethodon$land)
# pleth.pca <- gm.prcomp(Y.gpa$coords)
# pleth.pca.plot <- plot(pleth.pca)
# picknplot.shape(pleth.pca.plot) 
# May change arguments for plotRefToTarget
# picknplot.shape(plot(pleth.pca), method = "points", mag = 3, 
# links=plethodon$links)

# 2d with phylogeny
# data(plethspecies) 
# Y.gpa <- gpagen(plethspecies$land)
# gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
# pleth.phylo <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
# pleth.phylomorphospace <- plot(pleth.phylo, phylo = TRUE, cex = 2, 
# pch = 22, bg = gps, phylo.par = list(edge.color = "blue", 
# edge.width = 2, edge.lty = 2,
# node.pch = 22, node.bg = "black"))
# links.species <- plethodon$links[-11,]
# links.species[11, 1] <- 11
# picknplot.shape(pleth.phylomorphospace, method = "points", 
# links = links.species)

# 2d allometry 
# gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
# species = plethodon$species) 
# fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, 
# print.progress = FALSE)
# Predline
# PA <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
# method = "PredLine", pch = 19)
# picknplot.shape(PA)

# 3d and two-b-pls
# data("scallops")
# Y.gpa <- gpagen(scallops$coorddata, curves = scallops$curvslide, 
#              surfaces = scallops$surfslide)
# PLS <- two.b.pls(Y.gpa$coords, Y.gpa$Csize)
# PLS.plot = plot(PLS)
# picknplot.shape(PLS.plot)



cleanEx()
nameEx("plotAllSpecimens")
### * plotAllSpecimens

flush(stderr()); flush(stdout())

### Name: plotAllSpecimens
### Title: Plot landmark coordinates for all specimens
### Aliases: plotAllSpecimens
### Keywords: visualization

### ** Examples

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment

plotAllSpecimens(Y.gpa$coords,links=plethodon$links)



cleanEx()
nameEx("plotAllometry")
### * plotAllometry

flush(stderr()); flush(stdout())

### Name: plotAllometry
### Title: Plotting to assist visualization of shape-size covariation
###   (allometry)
### Aliases: plotAllometry
### Keywords: utilities

### ** Examples


# Simple allometry
data(plethodon) 
Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment  

gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
species = plethodon$species) 
fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, 
print.progress = FALSE)

# Predline
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
method = "PredLine", pch = 19)

# same as
logSize <- log(gdf$Csize)
plot(fit, type = "regression", reg.type = "PredLine", 
predictor = logSize, pch = 19)

# RegScore
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
method = "RegScore", pch = 19)

# same as
plot(fit, type = "regression", reg.type = "RegScore", 
predictor = logSize, pch = 19)

# CAC
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
method = "CAC", pch = 19)

# same (first plot) as
PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
plot(PLS)

# Group Allometries
fit <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=0, 
print.progress = FALSE)

# CAC (should not change from last time; model change has no effect)
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "CAC", 
pch = 19)

# Predline
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))

# RegScore
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))

# Size-Shape PCA

pc.plot <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
method = "size.shape", 
pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
summary(pc.plot$size.shape.PCA)

# Are species' shape differences just a manifestation of shape allometry?

fit3 <- procD.lm(coords ~ species, data = gdf, iter = 0, 
print.progress = FALSE)
plotAllometry(fit3, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
pch = 19, col = as.numeric(gdf$species))

# No evidence this is the case




cleanEx()
nameEx("plotOutliers")
### * plotOutliers

flush(stderr()); flush(stdout())

### Name: plotOutliers
### Title: Find potential outliers
### Aliases: plotOutliers
### Keywords: utilities

### ** Examples

data(plethodon)
# let's make some outliers
newland <- plethodon$land
newland[c(1,8),,2] <- newland[c(8,1),,2]
newland[c(3,11),,26] <- newland[c(11,3),,2]
Y<- gpagen(newland) # GPA
out <- plotOutliers(Y$coords) # function returns dimnames and address 
# of all specimens ordered
plotOutliers(Y$coords, inspect.outliers = TRUE) # function also produces 
# plots of identified outlier specimens compared to the mean shape

# example with groups
plotOutliers(Y$coords, groups = plethodon$species, 
inspect.outliers = TRUE)
 



cleanEx()
nameEx("plotRefToTarget")
### * plotRefToTarget

flush(stderr()); flush(stdout())

### Name: plotRefToTarget
### Title: Plot shape differences between a reference and target specimen
### Aliases: plotRefToTarget
### Keywords: visualization

### ** Examples

# Not Run
# Two dimensional data
# data(plethodon) 
# Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
# ref<-mshape(Y.gpa$coords)
# plotRefToTarget(ref,Y.gpa$coords[,,39])
# plotRefToTarget(ref,Y.gpa$coords[,,39], mag=2, outline=plethodon$outline)   
#magnify by 2X
# plotRefToTarget(ref,Y.gpa$coords[,,39], method="vector", mag=3)
# plotRefToTarget(ref,Y.gpa$coords[,,39], method="points", 
# outline=plethodon$outline)
# plotRefToTarget(ref,Y.gpa$coords[,,39], method="vector", 
# outline=plethodon$outline, mag=2.5)
# plotRefToTarget(ref,Y.gpa$coords[,,39], 
# gridPars=gridPar(pt.bg = "green", pt.size = 1),
# method="vector",mag=3)

# Three dimensional data
# data(scallops)
# Y.gpa<-gpagen(A=scallops$coorddata, curves=scallops$curvslide, 
# surfaces=scallops$surfslide)
# ref<-mshape(Y.gpa$coords)
# plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
# scallinks <- matrix(c(1,rep(2:16, each=2),1), nrow=16, byrow=TRUE)
# plotRefToTarget(ref,Y.gpa$coords[,,1],
# gridPars=gridPar(tar.pt.bg = "blue", tar.link.col="blue",
# tar.link.lwd=2), method="points", links = scallinks)




cleanEx()
nameEx("plotspec")
### * plotspec

flush(stderr()); flush(stdout())

### Name: plotspec
### Title: Plot 3D specimen, fixed landmarks and surface semilandmarks
### Aliases: plotspec
### Keywords: visualization

### ** Examples


## Not Run
# data(scallopPLY)
# ply <- scallopPLY$ply
# digitdat <- scallopPLY$coords
# plotspec(spec = ply, digitspec = digitdat, fixed = 16, 
# centered = TRUE, fixed.pt.col = "red", 
# fixed.pt.size = 15, col = "blue", size = 5)



cleanEx()
nameEx("procD.lm")
### * procD.lm

flush(stderr()); flush(stdout())

### Name: procD.lm
### Title: Procrustes ANOVA/regression for Procrustes shape variables
### Aliases: procD.lm
### Keywords: analysis

### ** Examples

## Not run: 
##D ### ANOVA example for Goodall's F test (multivariate shape vs. factors)
##D data(plethodon) 
##D Y.gpa <- gpagen(plethodon$land)    #GPA-alignment  
##D gdf <- geomorph.data.frame(Y.gpa, 
##D site = plethodon$site, 
##D species = plethodon$species) # geomorph data frame
##D 
##D fit1 <- procD.lm(coords ~ species * site, 
##D data = gdf, iter = 999, turbo = TRUE,
##D RRPP = FALSE, print.progress = FALSE) # randomize raw values
##D fit2 <- procD.lm(coords ~ species * site, 
##D data = gdf, iter = 999, turbo = TRUE,
##D RRPP = TRUE, print.progress = FALSE) # randomize residuals
##D 
##D summary(fit1)
##D summary(fit2)
##D 
##D # RRPP example applications
##D 
##D coef(fit2)
##D coef(fit2, test = TRUE)
##D anova(fit2) # same as summary
##D anova(fit2, effect.type = "Rsq")
##D # if species and site were modeled as random effects ...
##D anova(fit2, error = c("species:site", "species:site", "Residuals"))  
##D # not run, because it is a large object to print 
##D # remove # to run
##D # predict(fit2) 
##D 
##D # TPS plots for fitted values of some specimens
##D 
##D M <- Y.gpa$consensus
##D plotRefToTarget(M, fit2$GM$fitted[,,1], mag = 3)
##D plotRefToTarget(M, fit2$GM$fitted[,,20], mag = 3)
##D 
##D ### THE FOLLOWING ILLUSTRATES SIMPLER SOLUTIONS FOR 
##D ### THE NOW DEPRECATED advanced.procD.lm FUNCTION AND
##D ### PERFORM ANALYSES ALSO FOUND VIA THE morphol.disparity FUNCTION
##D ### USING THE pairwise FUNCTION
##D 
##D # Comparison of LS means, with log(Csize) as a covariate
##D 
##D # Model fits
##D fit.null <- procD.lm(coords ~ log(Csize) + species + site, data = gdf, 
##D iter = 999, print.progress = FALSE, turbo = TRUE)
##D fit.full <- procD.lm(coords ~ log(Csize) + species * site, data = gdf, 
##D iter = 999, print.progress = FALSE, turbo = TRUE)
##D 
##D # ANOVA, using anova.lm.rrpp function from the RRPP package 
##D # (replaces advanced.procD.lm)
##D anova(fit.null, fit.full, print.progress = FALSE)
##D 
##D # Pairwise tests, using pairwise function from the RRPP package
##D gp <-  interaction(gdf$species, gdf$site)
##D 
##D PW <- pairwise(fit.full, groups = gp, covariate = NULL)
##D 
##D # Pairwise distances between means, summarized two ways 
##D # (replaces advanced.procD.lm):
##D summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE)
##D summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)
##D 
##D # Pairwise comaprisons of group variances, two ways 
##D # (same as morphol.disaprity):
##D summary(PW, test.type = "var", confidence = 0.95, stat.table = TRUE)
##D summary(PW, test.type = "var", confidence = 0.95, stat.table = FALSE)
##D morphol.disparity(fit.full, groups = gp, iter = 999)
##D 
##D ### Regression example
##D data(ratland)
##D rat.gpa<-gpagen(ratland)         #GPA-alignment
##D gdf <- geomorph.data.frame(rat.gpa) # geomorph data frame is easy 
##D # without additional input
##D 
##D fit <- procD.lm(coords ~ Csize, data = gdf, iter = 999, turbo = TRUE,
##D RRPP = TRUE, print.progress = FALSE) 
##D summary(fit)
##D 
##D ### Extracting objects and plotting options
##D # diagnostic plots
##D plot(fit, type = "diagnostics") 
##D # diagnostic plots, including plotOutliers
##D plot(fit, type = "diagnostics", outliers = TRUE) 
##D 
##D # PC plot rotated to major axis of fitted values
##D plot(fit, type = "PC", pch = 19, col = "blue") 
##D # Regression-type plots 
##D 
##D # Use fitted values from the model to make prediction lines
##D plot(fit, type = "regression", 
##D predictor = gdf$Csize, reg.type = "RegScore", 
##D pch = 19, col = "green")
##D 
##D # Uses coefficients from the model to find the projected regression 
##D # scores
##D rat.plot <- plot(fit, type = "regression", 
##D predictor = gdf$Csize, reg.type = "RegScore", 
##D pch = 21, bg = "yellow") 
##D 
##D # TPS grids for min and max scores in previous plot
##D preds <- shape.predictor(fit$GM$fitted, x = rat.plot$RegScore, 
##D                         predmin = min(rat.plot$RegScore), 
##D                         predmax = max(rat.plot$RegScore))
##D M <- rat.gpa$consensus
##D plotRefToTarget(M, preds$predmin, mag=2)
##D plotRefToTarget(M, preds$predmax, mag=2)
##D                         
##D attributes(fit)
##D fit$fitted[1:3, ] # the fitted values (first three specimens)
##D fit$GM$fitted[,, 1:3] # the fitted values as Procrustes 
##D # coordinates (same specimens)
##D 
##D ### THE FOLLOWING ILLUSTRATES SIMPLER SOLUTIONS FOR 
##D ### THE NOW DEFUNCT nested.update FUNCTION USING
##D ### THE anova GENERIC FUNCTION
##D 
##D data("larvalMorph")
##D Y.gpa <- gpagen(larvalMorph$tailcoords, 
##D curves = larvalMorph$tail.sliders,
##D ProcD = TRUE, print.progress = FALSE)
##D gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
##D family = larvalMorph$family)
##D 
##D fit <- procD.lm(coords ~ treatment/family, data = gdf, turbo = TRUE,
##D print.progress = FALSE, iter = 999)
##D anova(fit) # treatment effect not adjusted
##D anova(fit, error = c("treatment:family", "Residuals")) 
##D # treatment effect updated (adjusted)
## End(Not run)






cleanEx()
nameEx("procD.pgls")
### * procD.pgls

flush(stderr()); flush(stdout())

### Name: procD.pgls
### Title: Phylogenetic ANOVA/regression for Procrustes shape variables
### Aliases: procD.pgls
### Keywords: analysis

### ** Examples

### Example of D-PGLS for high-dimensional data 
data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy)

pleth.pgls <- procD.pgls(coords ~ Csize, phy = phy, data = gdf, 
iter = 999)
anova(pleth.pgls)
summary(pleth.pgls)  #similar output

### Working with procD.pgls objects
predict(pleth.pgls)
plot(pleth.pgls, type="regression", reg.type="RegScore", 
predictor = gdf$Csize)
attributes(pleth.pgls) # Note the PGLS object
attributes(pleth.pgls$PGLS) # PGLS details embedded within PGLS object
pleth.pgls$LM$Pcov # the projection matrix derived from the 
# phylogenetic covariance matrix
pleth.pgls$pgls.fitted # the PGLS fitted values 
pleth.pgls$GM$pgls.fitted # The same fitted values, in a 3D array

# Changing lambda value

pleth.pgls2 <- procD.pgls(coords ~ Csize, phy = phy, lambda = 0.5, 
data = gdf, iter = 999)

anova(pleth.pgls)
anova(pleth.pgls2)



cleanEx()
nameEx("read.ply")
### * read.ply

flush(stderr()); flush(stdout())

### Name: read.ply
### Title: Read mesh data (vertices and faces) from ply files
### Aliases: read.ply
### Keywords: IO

### ** Examples

# If the file has no mesh color, or color is undesirable, user can 
# assign this as follows:
# Using the example scallop PLY
data(scallopPLY) 
myply <- scallopPLY$ply
myply$material$color <- "gray" # using color word
myply$material$color <- "#FCE6C9" # using RGB code



cleanEx()
nameEx("readland.shapes")
### * readland.shapes

flush(stderr()); flush(stdout())

### Name: readland.shapes
### Title: Read landmark data from a shapes object (StereoMorph)
### Aliases: readland.shapes
### Keywords: IO

### ** Examples

# A true example is not possible, as digitizing experiences are 
# unique, but here is a general set-up
# myShapes <- readShapes("myDigitizingFile") # data from readShapes 
# from StereoMorph 
# myGMdata <- readland.shapes(myShapes) # just reading in the fixed 
# landmarks
# myGMdata <- readland.shapes(myShapes, 
#      nCurvePts = c(10, 15, 5), 
#      continuous.curve = 2) # fixed landmarks plus curve points 
# for three curves, one closed
# myGPA <- gpagen(myGMdata, ProcD = FALSE) # GPA perfomed with 
# minimized bending energy



cleanEx()
nameEx("rotate.coords")
### * rotate.coords

flush(stderr()); flush(stdout())

### Name: rotate.coords
### Title: Rotate or flip landmark or coordinate configurations
### Aliases: rotate.coords
### Keywords: utilities

### ** Examples

data(plethodon)
Y.gpa <- gpagen(plethodon$land)
plot(Y.gpa)
Y.gpa2 <- rotate.coords(Y.gpa, "flipX")
plot(Y.gpa2)
Y.gpa3 <- rotate.coords(Y.gpa2, "rotateCC")
plot(Y.gpa3)

spec1 <- Y.gpa$coords[,,1]
plot(spec1, asp = 1)
spec1 <- rotate.coords(spec1, "flipY")
plot(spec1, asp = 1)

specs1to3 <- Y.gpa$coords[,,1:3]
plotAllSpecimens(specs1to3)
specs1to3 <- rotate.coords(specs1to3, "rotateC")
plotAllSpecimens(specs1to3)



cleanEx()
nameEx("shape.predictor")
### * shape.predictor

flush(stderr()); flush(stdout())

### Name: shape.predictor
### Title: Shape prediction from numeric predictors
### Aliases: shape.predictor

### ** Examples

# Examples using Plethodon data
# NOT RUN
# data("plethodon")

# Y.gpa <- gpagen(plethodon$land)    #GPA-alignment    
# plot(gm.prcomp(Y.gpa$coords))

# preds <- shape.predictor(Y.gpa$coords, x= NULL, Intercept = FALSE, 
# pred1 = -0.1, pred2 = 0.1) # PC 1 extremes, sort of
# M <- mshape(Y.gpa$coords)
# plotRefToTarget(M, preds$pred1)
# plotRefToTarget(M, preds[[1]]) # same result
# plotRefToTarget(M, preds$pred2)

# PCA <- gm.prcomp(Y.gpa$coords)
# PC <- PCA$x[,1]
# preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
# pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
# plotRefToTarget(M, preds$pred1)
# plotRefToTarget(M, preds$pred2)

# PC <- PCA$x[,1:2]
# user-picked spots can be anything, but it in this case, apparent groups
# preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
#                          pred1 = c(0.045,-0.02), 
#                          pred2 = c(-0.025,0.06), 
#                          pred3 = c(-0.06,-0.04)) 
# plotRefToTarget(M, preds$pred1)
# plotRefToTarget(M, preds$pred2)
# plotRefToTarget(M, preds$pred3)

# allometry example - straight-up allometry

# preds <- shape.predictor(Y.gpa$coords, x= log(Y.gpa$Csize), 
#                          Intercept = TRUE, 
#                          predmin = min(log(Y.gpa$Csize)), 
#                          predmax = max(log(Y.gpa$Csize))) 

# plotRefToTarget(M, preds$predmin, mag=3)
# plotRefToTarget(M, preds$predmax, mag=3)

# allometry example - using RegScore or PredLine via procD.lm

# gdf <- geomorph.data.frame(Y.gpa)
# plethAllometry <- procD.lm(coords ~ log(Csize), data=gdf)
# allom.plot <- plot(plethAllometry, 
# type = "regression", 
# predictor = log(gdf$Csize),
# reg.type ="RegScore") # make sure to have a predictor 

# preds <- shape.predictor(plethAllometry$GM$fitted, 
#                          x= allom.plot$RegScore, Intercept = FALSE, 
#                          predmin = min(allom.plot$RegScore), 
#                          predmax = max(allom.plot$RegScore)) 
# plotRefToTarget(M, preds$predmin, mag=3)
# plotRefToTarget(M, preds$predmax, mag=3)

# allom.plot <- plot(plethAllometry, 
# type = "regression", 
# predictor = log(gdf$Csize),
# reg.type ="PredLine")
# preds <- shape.predictor(plethAllometry$GM$fitted, 
#                          x= allom.plot$PredLine, Intercept = FALSE, 
#                          predmin = min(allom.plot$PredLine), 
#                          predmax = max(allom.plot$PredLine)) 
# plotRefToTarget(M, preds$predmin, mag=3)
# plotRefToTarget(M, preds$predmax, mag=3)

# using factors via PCA

# gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
#         site = plethodon$site)
# pleth <- procD.lm(coords ~ species*site, data=gdf)
# PCA <- prcomp(pleth$fitted)
# plot(PCA$x, asp=1, pch=19)

# means <- unique(round(PCA$x,3))
# means # note: suggests 3 PCs useful enough

# preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x= PCA$x[,1:3],
#                          Intercept = FALSE,
#                          pred1 = means[1,1:3], 
#                          pred2 = means[2,1:3],
#                          pred3 = means[3,1:3], 
#                          pred4 = means[4,1:3])                   
# plotRefToTarget(M, preds$pred1, mag=2)
# plotRefToTarget(M, preds$pred2, mag=2)
# plotRefToTarget(M, preds$pred3, mag=2)
# plotRefToTarget(M, preds$pred4, mag=2)

# Using a design matrix for factors

# X <- pleth$X
# X # includes intercept; remove for better functioning 
# X <- X[,-1]
# symJord <- c(0,1,0) # design for P. Jordani in sympatry
# alloJord <- c(0,0,0) # design for P. Jordani in allopatry
# preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x = X, 
#                          Intercept = TRUE, 
#                          symJord=symJord, alloJord=alloJord)
# plotRefToTarget(M, preds$symJord, mag=2)
# plotRefToTarget(M, preds$alloJord, mag=2)

# PLS Example

# data(plethShapeFood) 
# Y.gpa<-gpagen(plethShapeFood$land)    #GPA-alignment    

# 2B-PLS between head shape and food use data
# PLS <-two.b.pls(A1 = plethShapeFood$food, A2 = Y.gpa$coords, iter=999) 
# summary(PLS)
# plot(PLS)

# preds <- shape.predictor(Y.gpa$coords, plethShapeFood$food, 
#                          Intercept = FALSE,
#                          method = "PLS",
#                          pred1 = 2, pred2 = -4, pred3 = 2.5) 
                        # using PLS plot as a guide
# M <- mshape(Y.gpa$coords)
# plotRefToTarget(M, preds$pred1, mag=2)
# plotRefToTarget(M, preds$pred2, mag=2)
# plotRefToTarget(M, preds$pred3, mag=2)




cleanEx()
nameEx("shapeHulls")
### * shapeHulls

flush(stderr()); flush(stdout())

### Name: shapeHulls
### Title: Update Plots with Convex Hulls for Groups
### Aliases: shapeHulls
### Keywords: utilities

### ** Examples


# Via procD.lm and plot.procD.lm

data("pupfish")
gdf <- geomorph.data.frame(coords = pupfish$coords, Sex = pupfish$Sex,
Pop = pupfish$Pop)
fit <- procD.lm(coords ~ Pop * Sex, data = gdf, print.progress = FALSE)
pc.plot <- plot(fit, type = "PC", pch = 19)
shapeHulls(pc.plot)

pc.plot <- plot(fit, type = "PC", pch = 19)
groups <- interaction(gdf$Pop, gdf$Sex)

shapeHulls(pc.plot, groups = groups, 
group.cols = c("dark red", "dark red", "dark blue", "dark blue"),
group.lwd = rep(2, 4), group.lty = c(2, 1, 2, 1))

legend("topright", levels(groups), 
col = c("dark red", "dark red", "dark blue", "dark blue"),
lwd = rep(2,4), lty = c(2, 1, 2, 1))

pc.plot <- plot(fit, type = "PC", pch = 19)
shapeHulls(pc.plot, groups = gdf$Sex, group.cols = c("black", "black"), 
group.lwd = rep(2, 2), group.lty = c(2, 1))
legend("topright", levels(gdf$Sex), lwd = 2, lty = c(2, 1))

# Via gm.prcomp and plot.gm.prcomp

data(plethspecies) 
Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
pleth.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
summary(pleth.phylo)

pc.plot <- plot(pleth.phylo, phylo = TRUE)
gp <- factor(c(rep(1, 5), rep(2, 4)))
shapeHulls(pc.plot, groups = gp, group.cols = 1:2, 
group.lwd = rep(2, 2), group.lty = c(2, 1))
legend("topright", c("P. cinereus clade", "P. hubrichti clade"), 
col = 1:2, lwd = 2, lty = c(2, 1))



cleanEx()
nameEx("two.b.pls")
### * two.b.pls

flush(stderr()); flush(stdout())

### Name: two.b.pls
### Title: Two-block partial least squares analysis for Procrustes shape
###   variables
### Aliases: two.b.pls
### Keywords: analysis

### ** Examples

data(plethShapeFood) 
Y.gpa<-gpagen(plethShapeFood$land)    #GPA-alignment    

#2B-PLS between head shape and food use data
PLS <-two.b.pls(Y.gpa$coords,plethShapeFood$food,iter=999)
summary(PLS)
P <- plot(PLS)
 
 # Visualize shape at minimum and maximum PLS scores.
 # This is the challenging way
 
 # Block 1 (shape)
 minx <- min(P$plot.args$x)
 maxx <- max(P$plot.args$x)
 preds <- shape.predictor(P$A1, 
 x = P$plot.args$x,
 min = minx, max = maxx)
 plotRefToTarget(mshape(P$A1), preds$min)
 plotRefToTarget(mshape(P$A1), preds$max)
 
 ### Visualize shape variation using picknplot.shape Because picknplot  
 ### requires user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(P)
 




cleanEx()
nameEx("two.d.array")
### * two.d.array

flush(stderr()); flush(stdout())

### Name: two.d.array
### Title: Convert (p x k x n) data array into 2D data matrix
### Aliases: two.d.array
### Keywords: utilities

### ** Examples

data(plethodon) 
plethodon$land    #original data in the form of 3D array

two.d.array(plethodon$land)   # Convert to a 2D data matrix



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
