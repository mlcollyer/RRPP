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

x<-matrix(rnorm(18),nrow=3)  # Random triangles (all coordinates on same row for each triangle)
arrayspecs(x,3,2) 
 
x2<-matrix(rnorm(18),ncol=2) # Random triangles (each landmark on its own row)
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

data(mosquito)
gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, ind=mosquito$ind, side=mosquito$side,
replicate=mosquito$replicate)
mosquito.sym <- bilat.symmetry(A = wingshape, ind = ind, side = side,
replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 499, data = gdf)
summary(mosquito.sym)
plot(mosquito.sym, warpgrids = TRUE)
mosquito.sym$shape.anova # extract just the anova table on shape

# Previous example, performing GPA first
Y.gpa <- gpagen(mosquito$wingshape)
mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = ind, side = side,
replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 499, data = gdf)
summary(mosquito.sym2)
summary(mosquito.sym) # same results

#Example of object symmetry

data(scallops)
gdf <- geomorph.data.frame(shape = scallops$coorddata, ind=scallops$ind)
scallop.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 499)
summary(scallop.sym)

# Previous example, incorporating semilandmarks

scallop.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
curves= scallops$curvslide, surfaces = scallops$surfslide,
land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 499)
summary(scallop.sym)
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

data(larvalMorph) 
head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders)
tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
comb.lm <- combine.subsets(head = head.gpa, tail = tail.gpa, gpa = TRUE)
summary(comb.lm)
# Plot first specimen and color code landmarks 
# (configurations are actual relative size)
comb.lm$coords[,,1]
plot(comb.lm$coords[,,1], pch = 21, bg = c(rep(1,26), rep(2,64)), asp = 1)

# choose to override scaling by relative size 
comb.lm <- combine.subsets(head = head.gpa$coords, 
tail = tail.gpa$coords, gpa = FALSE, CS.sets = NULL)
summary(comb.lm)
plot(comb.lm$coords[,,1], pch = 21, bg = c(rep(1,26), rep(2,64)), asp = 1)



cleanEx()
nameEx("compare.CR")
### * compare.CR

flush(stderr()); flush(stdout())

### Name: compare.CR
### Title: Comparisons of Effect Sizes from Modularity Analyses
### Aliases: compare.CR
### Keywords: analysis

### ** Examples


### NOT RUN
# Example 1: Compare modular signal across datasets

# data(pupfish) 
# Y.gpa<-gpagen(pupfish$coords, print.progress = FALSE)    #GPA-alignment    
## landmarks on the body and operculum
# land.gps<-rep('a',56); land.gps[39:48]<-'b'

# group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
# levels(group)

# coords.gp <- coords.subset(Y.gpa$coords, group)

# modularity.tests <- lapply(1:nlevels(group), function(j) modularity.test(coords.gp[[j]],
#  land.gps, iter = 499, print.progress = FALSE))
## the lapply function performs the modularity test on each 3D array in the lists provided

# modularity.tests[[1]]
# modularity.tests[[2]]
# modularity.tests[[3]]
# modularity.tests[[4]]

# group.Z <- compare.CR(modularity.tests, CR.null = FALSE)
# group.Z  ## NOTE: need a summary function 

# Example 2: Compare alternative modular hypotheses

# land.gps3 <- rep('a',56); land.gps3[39:48]<-'b'; land.gps3[c(6:9,28:38)] <- 'c' 
   #3 module hypothesis (tail now a module)
# land.gps4 <- rep('a',56); land.gps4[39:48]<-'b'; land.gps4[c(6:9,28:38)] <- 'c'; 
 # land.gps4[c(10,49:56)] <- 'd'  #4 module hypothesis (eye now a module)

# m3.test <- modularity.test(coords.gp$Marsh.F,land.gps3, iter = 499, print.progress = FALSE)
# m4.test <- modularity.test(coords.gp$Marsh.F,land.gps4, iter = 499, print.progress = FALSE)

# model.Z <- compare.CR(modularity.tests[[1]],m3.test,m4.test, CR.null = TRUE)
# model.Z 




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

ER<-compare.evol.rates(A=Y.gpa$coords, phy=plethspecies$phy,method="simulation",gp=gp.end,iter=999)
summary(ER)
plot(ER)



cleanEx()
nameEx("compare.multi.evol.rates")
### * compare.multi.evol.rates

flush(stderr()); flush(stdout())

### Name: compare.multi.evol.rates
### Title: Comparing net rates of shape evolution among traits on
###   phylogenies
### Aliases: compare.multi.evol.rates
### Keywords: analysis

### ** Examples


data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
land.gp<-c("A","A","A","A","A","B","B","B","B","B","B")  #mandible and cranium subsets

EMR<-compare.multi.evol.rates(A=Y.gpa$coords,gp=land.gp, 
    Subset=TRUE, phy= plethspecies$phy,iter=999)
summary(EMR)
plot(EMR)



cleanEx()
nameEx("compare.pls")
### * compare.pls

flush(stderr()); flush(stdout())

### Name: compare.pls
### Title: Comparisons of Effect Sizes from Partial Least Squares
### Aliases: compare.pls
### Keywords: analysis

### ** Examples

# Example of comparative morphological integration between pupfish head and body shapes

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

 integ.tests <- Map(function(x,y) integration.test(x, y, iter=499, print.progress = FALSE), 
 head.coords.gp, tail.coords.gp)
# the map function performs the integration test on each 3D array in the lists provided

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
 # define.sliders(scallops$coorddata[,,1], nsliders=11,surfsliders = scallops$surfslide) 
 # here the first specimen is used for plotting purposes only
 
## Examples of AUTO mode 
 ## 1 curve of sliding semilandmark
 # Define sliders for scallopdata
 # sliders = define.sliders(c(5:16,1))

 ## 2 curves of sliding semilandmarks
 # Define sliders for 10 landmarks, where LMs 1, 5, and 10 fixed
 # 2, 3, and 4 are along a curve between 1 and 5
 # and 6, 7, 8, and 9 are along a curve between 5 and 10.
 # sliders = rbind(define.sliders(1:5), define.sliders(5:10)) 



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
#Articulation point is landmark 1, rotate mandibular landmarks (2-5) relative to cranium

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

gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)
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
### Title: Principal components analysis of shape data
### Aliases: gm.prcomp
### Keywords: visualization

### ** Examples

 data(plethspecies) 
 Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
 
 ### PCA 
 pleth.raw <- gm.prcomp(Y.gpa$coords)
 summary(pleth.raw)
 
 ### PCA with phylogeny (result is same as above, but with additional components)
 pleth.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
 summary(pleth.phylo)
 
 #### Plotting
 plot(pleth.raw)
 gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
 par(mar=c(2, 2, 2, 2))
 plot(pleth.raw, pch=22, cex = 1.5, bg = gps) # Modify options as desired
 #  Add things as desired using standard R plotting
 text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", pos = 4, font = 2)
 text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
 legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))
 
 ### Phylomorphospace plot
 plot(pleth.phylo, pch=21, bg=1:nrow(pleth.phylo$x), phylo = TRUE,
      phylo.par = list(edge.color="grey", node.cex=0)); title(main="phylomorphospace")
 text(pleth.phylo$x, labels = labels(pleth.phylo$x)[[1]],
      pos = 2, font = 4) 
 text(pleth.phylo$anc.x, labels = 1:nrow(pleth.phylo$anc.x),
      adj = c(-0.1, -0.1), font = 2) 
 
 ### Visualize shape variation using picknplot.shape Because picknplot requires 
 ### user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(plot(pleth.phylo))
 




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
Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)
summary(Y.gpa)
plot(Y.gpa)

# Example 2: points and semilandmarks on curves
data(hummingbirds)

###Slider matrix
hummingbirds$curvepts

# Using bending energy for sliding
Y.gpa <- gpagen(hummingbirds$land,curves=hummingbirds$curvepts,ProcD=FALSE)   
summary(Y.gpa)
plot(Y.gpa)


# Using Procrustes Distance for sliding
Y.gpa <- gpagen(hummingbirds$land,curves=hummingbirds$curvepts,ProcD=TRUE)   
summary(Y.gpa)
plot(Y.gpa)


# Example 3: points, curves and surfaces
data(scallops)

# Using Procrustes Distance for sliding
Y.gpa <- gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
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
GP1 <- gridPar(pt.bg = "red", pt.size = 1, link.col="blue", link.lwd=2, n.col.cell=50)
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP1, mag=2, 
links=plethodon$links, method="TPS")

# Altering point color
GP2 <- gridPar(pt.bg = "green", pt.size = 1) 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP2, mag=3, method="vector")

# Altering ref and target points
GP3 <- gridPar(pt.bg = "blue", pt.size = 1.5, tar.pt.bg = "orange", tar.pt.size = 1) 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP3, mag=3, method="points")

# Altering outline color
GP4 <- gridPar(tar.out.col = "red", tar.out.cex = 0.3) 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP4, mag=3, 
outline=plethodon$outline, method="TPS")

# Altering text labels
GP5 <- gridPar(txt.pos = 3, txt.col = "red") 
plotRefToTarget(ref,Y.gpa$coords[,,39], gridPars=GP5, mag=3, method="vector", label=TRUE)



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
plot(IT) # PLS plot

 
 ### Visualize shape variation using picknplot.shape Because picknplot requires 
 ### user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(plot(IT))

IT$left.pls.vectors # extracting just the left (first block) singular vectors



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
nameEx("modularity.test")
### * modularity.test

flush(stderr()); flush(stdout())

### Name: modularity.test
### Title: Evaluate the degree of modular signal in shape data
### Aliases: modularity.test
### Keywords: analysis

### ** Examples

data(pupfish) 
Y.gpa<-gpagen(pupfish$coords, print.progress = FALSE)    #GPA-alignment    
 #landmarks on the body and operculum
land.gps<-rep('a',56); land.gps[39:48]<-'b'

MT <- modularity.test(Y.gpa$coords,land.gps,CI=FALSE,iter=99)
summary(MT) # Test summary
plot(MT) # Histogram of CR sampling distribution 
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

data(plethodon)
Y.gpa<-gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment
gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
site = plethodon$site)

# Morphological disparity for entire data set
morphol.disparity(coords ~ 1, groups = NULL, data = gdf, 
iter = 999, print.progress = FALSE)

# Morphological disparity for entire data set, accounting for allometry
morphol.disparity(coords ~ Csize, groups= NULL, data = gdf, 
iter = 999, print.progress = FALSE)

# Morphological disparity without covariates, using overall mean
morphol.disparity(coords ~ 1, groups= ~ species*site, data = gdf, 
iter = 999, print.progress = FALSE)

# Morphological partial disparities for overal mean
morphol.disparity(coords ~ 1, groups= ~ species*site, partial = TRUE, 
data = gdf, iter = 999, print.progress = FALSE)

# Morphological disparity without covariates, using group means
morphol.disparity(coords ~ species*site, groups= ~species*site, 
data = gdf, iter = 999, print.progress = FALSE)

# Morphological disparity of different groups than those 
# described by the linear model
morphol.disparity(coords ~ Csize + species*site, groups= ~ species, 
data = gdf, iter = 999, print.progress = FALSE)

# Extracting components
MD <- morphol.disparity(coords ~ Csize + species*site, groups= ~ species, 
data = gdf, iter = 999, print.progress = FALSE)
MD$Procrustes.var # just the Procrustes variances


### Morphol.disparity can be used with previously-defined 
### procD.lm or lm.rrpp class objects

data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
names(gp.end)<-plethspecies$phy$tip

gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy, gp.end = gp.end)

pleth.ols <- procD.lm(coords ~ Csize + gp.end, 
data = gdf, iter = 999) # ordinary least squares
pleth.pgls <- procD.pgls(coords ~ Csize + gp.end, phy = phy, 
data = gdf, iter = 999) # phylogenetic generalized least squares

summary(pleth.ols)
summary(pleth.pgls)

morphol.disparity(f1 = pleth.ols, groups = ~ gp.end, data = gdf, 
iter = 999, print.progress = FALSE)
morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end, data = gdf, 
iter = 999, print.progress = FALSE)

 ### Morphol.disparity using RRPP
 PW <- pairwise(pleth.ols, groups = gp.end)
 summary(PW, test.type = 'var')




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

mshape(Y.gpa$coords)   #mean (consensus) configuration



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

IT<- phylo.integration(Y.gpa$coords,partition.gp=land.gps,phy=plethspecies$phy,iter=999)
summary(IT) # Test summary
plot(IT) # PLS plot

 ### Visualize shape variation using picknplot.shape Because picknplot requires 
 ### user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(plot(IT))




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

data(plethspecies)
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 

MT <- phylo.modularity(Y.gpa$coords, partition.gp=land.gps, phy=plethspecies$phy, 
CI = FALSE, iter=499)
summary(MT) # Test summary
plot(MT) # Histogram of CR sampling distribution 



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

#Test for phylogenetic signal in size
PS.size <- physignal(A=Y.gpa$Csize,phy=plethspecies$phy,iter=999)
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
# picknplot.shape(plot(pleth.pca), method = "points", mag = 3, links=plethodon$links)

# 2d with phylogeny
# data(plethspecies) 
# Y.gpa <- gpagen(plethspecies$land)
# gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
# pleth.phylo <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
# pleth.phylomorphospace <- plot(pleth.phylo, phylo = TRUE, cex = 2, pch = 22, 
# bg = gps, phylo.par = list(edge.color = "blue", edge.width = 2, edge.lty = 2,
# node.pch = 22, node.bg = "black"))
# links.species <- plethodon$links[-11,]
# links.species[11, 1] <- 11
# picknplot.shape(pleth.phylomorphospace, method = "points", links = links.species)

# 2d allometry 
# gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
# species = plethodon$species) 
# fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, print.progress = FALSE)
# Predline
# PA <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", pch = 19)
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
fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, print.progress = FALSE)

# Predline
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", pch = 19)

# same as
logSize <- log(gdf$Csize)
plot(fit, type = "regression", reg.type = "PredLine", predictor = logSize, pch = 19)

# RegScore
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "RegScore", pch = 19)

# same as
plot(fit, type = "regression", reg.type = "RegScore", predictor = logSize, pch = 19)

# CAC
plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "CAC", pch = 19)

# same (first plot) as
PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
plot(PLS)

# Group Allometries
fit2 <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=0, print.progress = FALSE)

# CAC (should not change from last time; model change has no effect)
plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "CAC", pch = 19)

# Predline
plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))

# RegScore
plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))

# Size-Shape PCA

pc.plot <- plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "size.shape", 
pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
summary(pc.plot$size.shape.PCA)

# Are species' shape differences just a manifestation of shape allometry?

fit3 <- procD.lm(coords ~ species, data=gdf, iter=0, print.progress = FALSE)
plotAllometry(fit3, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
pch = 19, col = as.numeric(gdf$species))

# No evidence this is the case




cleanEx()
nameEx("plotGMPhyloMorphoSpace")
### * plotGMPhyloMorphoSpace

flush(stderr()); flush(stdout())

### Name: plotGMPhyloMorphoSpace
### Title: Plot phylogenetic tree and specimens in tangent space
### Aliases: plotGMPhyloMorphoSpace
### Keywords: visualization

### ** Examples

data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    

plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords)
plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords, 
                 plot.param=list(t.bg="blue",txt.col="red",n.cex=1))
#NOTE: 3D plot also available: plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords, zaxis= "time",
#                 plot.param=list(n.cex=2, n.bg="blue"), shadow=TRUE)



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
out <- plotOutliers(Y$coords) # function returns dimnames and address of all specimens ordered
plotOutliers(Y$coords, inspect.outliers = TRUE) # function also produces plots of 
# identified outlier specimens compared to the mean shape

# example with groups
plotOutliers(Y$coords, groups = plethodon$species, inspect.outliers = TRUE)
 



cleanEx()
nameEx("plotRefToTarget")
### * plotRefToTarget

flush(stderr()); flush(stdout())

### Name: plotRefToTarget
### Title: Plot shape differences between a reference and target specimen
### Aliases: plotRefToTarget
### Keywords: visualization

### ** Examples

# Two dimensional data
data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
ref<-mshape(Y.gpa$coords)
plotRefToTarget(ref,Y.gpa$coords[,,39])
plotRefToTarget(ref,Y.gpa$coords[,,39],mag=2,outline=plethodon$outline)   #magnify by 2X
plotRefToTarget(ref,Y.gpa$coords[,,39],method="vector",mag=3)
plotRefToTarget(ref,Y.gpa$coords[,,39],method="points",outline=plethodon$outline)
plotRefToTarget(ref,Y.gpa$coords[,,39],gridPars=gridPar(pt.bg = "green", pt.size = 1),
method="vector",mag=3)

# Three dimensional data
# data(scallops)
# Y.gpa<-gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide)
# ref<-mshape(Y.gpa$coords)
# plotRefToTarget(ref,Y.gpa$coords[,,1],method="points")
# scallinks <- matrix(c(1,rep(2:16, each=2),1), nrow=16, byrow=TRUE)
# plotRefToTarget(ref,Y.gpa$coords[,,1],gridPars=gridPar(tar.pt.bg = "blue", tar.link.col="blue",
# tar.link.lwd=2), method="points", links = scallinks)




cleanEx()
nameEx("plotTangentSpace")
### * plotTangentSpace

flush(stderr()); flush(stdout())

### Name: plotTangentSpace
### Title: Plot specimens in tangent space
### Aliases: plotTangentSpace
### Keywords: visualization

### ** Examples

data(plethodon) 
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment

gp <- interaction(plethodon$species, plethodon$site) # group must be a factor
plotTangentSpace(Y.gpa$coords, groups = gp) 

## To save and use output
PCA <- plotTangentSpace(Y.gpa$coords, groups = gp, legend=TRUE) 
summary(PCA)
PCA$pc.shapes
PCA$rotation

##To change colors of groups
col.gp <- rainbow(length(levels(gp))) 
   names(col.gp) <- levels(gp)
col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must NOT be a factor
plotTangentSpace(Y.gpa$coords, groups = col.gp)

## To plot residual shapes from an allometry regression (note: must add mean back in!) 
plotTangentSpace(arrayspecs(resid(lm(two.d.array(Y.gpa$coords)~log(Y.gpa$Csize)))+
         predict(lm(two.d.array(Y.gpa$coords)~1)),12,2))



cleanEx()
nameEx("plotspec")
### * plotspec

flush(stderr()); flush(stdout())

### Name: plotspec
### Title: Plot 3D specimen, fixed landmarks and surface semilandmarks
### Aliases: plotspec
### Keywords: visualization

### ** Examples

# data(scallopPLY)
# ply <- scallopPLY$ply
# digitdat <- scallopPLY$coords
# plotspec(spec=ply,digitspec=digitdat,fixed=16, centered =TRUE)



cleanEx()
nameEx("procD.lm")
### * procD.lm

flush(stderr()); flush(stdout())

### Name: procD.lm
### Title: Procrustes ANOVA/regression for Procrustes shape variables
### Aliases: procD.lm
### Keywords: analysis

### ** Examples

### ANOVA example for Goodall's F test (multivariate shape vs. factors)
data(plethodon) 
Y.gpa <- gpagen(plethodon$land)    #GPA-alignment  
gdf <- geomorph.data.frame(Y.gpa, 
site = plethodon$site, 
species = plethodon$species) # geomorph data frame

fit1 <- procD.lm(coords ~ species * site, 
data = gdf, iter = 999, 
RRPP = FALSE, print.progress = FALSE) # randomize raw values
fit2 <- procD.lm(coords ~ species * site, 
data = gdf, iter = 999, 
RRPP = TRUE, print.progress = FALSE) # randomize residuals

summary(fit1)
summary(fit2)

# RRPP example applications

coef(fit2)
coef(fit2, test = TRUE)
anova(fit2) # same as summary
anova(fit2, effect.type = "Rsq")
# if species and site were modeled as random effects ...
anova(fit2, error = c("species:site", "species:site", "Residuals"))  
# not run, because it is a large object to print 
# remove # to run
# predict(fit2) 

# TPS plots for fitted values of some specimens

M <- Y.gpa$consensus
plotRefToTarget(M, fit2$GM$fitted[,,1], mag = 3)
plotRefToTarget(M, fit2$GM$fitted[,,20], mag = 3)

### THE FOLLOWING ILLUSTRATES SIMPLER SOLUTIONS FOR 
### THE NOW DEPRECATED advanced.procD.lm FUNCTION AND
### PERFORM ANALYSES ALSO FOUND VIA THE morphol.disparity FUNCTION
### USING THE pairwise FUNCTION

# Comparison of LS means, with log(Csize) as a covariate

# Model fits
fit.null <- procD.lm(coords ~ log(Csize) + species + site, data = gdf, 
iter = 999, print.progress = FALSE)
fit.full <- procD.lm(coords ~ log(Csize) + species * site, data = gdf, 
iter = 999, print.progress = FALSE)

# ANOVA, using anova.lm.rrpp function from the RRPP package (replaces advanced.procD.lm)
anova(fit.null, fit.full, print.progress = FALSE)

# Pairwise tests, using pairwise function from the RRPP package
gp <-  interaction(gdf$species, gdf$site)

PW <- pairwise(fit.full, groups = gp, covariate = NULL)

# Pairwise distances between means, summarized two ways (replaces advanced.procD.lm):
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)

# Pairwise comaprisons of group variances, two ways (same as morphol.disaprity):
summary(PW, test.type = "var", confidence = 0.95, stat.table = TRUE)
summary(PW, test.type = "var", confidence = 0.95, stat.table = FALSE)
morphol.disparity(fit.full, groups = gp, iter=999)

### Regression example
data(ratland)
rat.gpa<-gpagen(ratland)         #GPA-alignment
gdf <- geomorph.data.frame(rat.gpa) # geomorph data frame is easy without additional input

fit <- procD.lm(coords ~ Csize, data = gdf, iter = 999, 
RRPP = TRUE, print.progress = FALSE) 
summary(fit)

### Extracting objects and plotting options
# diagnostic plots
plot(fit, type = "diagnostics") 
# diagnostic plots, including plotOutliers
plot(fit, type = "diagnostics", outliers = TRUE) 

# PC plot rotated to major axis of fitted values
plot(fit, type = "PC", pch = 19, col = "blue") 
# Regression-type plots

# Use fitted values from the model to make prediction lines
plot(fit, type = "regression", 
predictor = gdf$Csize, reg.type = "RegScore", 
pch = 19, col = "green")

# Uses coefficients from the model to find the projected regression scores
rat.plot <- plot(fit, type = "regression", 
predictor = gdf$Csize, reg.type = "RegScore", 
pch = 21, bg = "yellow") 

# TPS grids for min and max scores in previous plot
preds <- shape.predictor(fit$GM$fitted, x = rat.plot$RegScore, 
                        predmin = min(rat.plot$RegScore), 
                        predmax = max(rat.plot$RegScore))
M <- rat.gpa$consensus
plotRefToTarget(M, preds$predmin, mag=2)
plotRefToTarget(M, preds$predmax, mag=2)
                        
attributes(fit)
fit$fitted[1:3, ] # the fitted values (first three specimens)
fit$GM$fitted[,, 1:3] # the fitted values as Procrustes coordinates (same specimens)

### THE FOLLOWING ILLUSTRATES SIMPLER SOLUTIONS FOR 
### THE NOW DEPRECATED nested.update FUNCTION USING
### THE anova GENERIC FUNCTION

data("larvalMorph")
Y.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders,
ProcD = TRUE, print.progress = FALSE)
gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
family = larvalMorph$family)

fit <- procD.lm(coords ~ treatment/family, data = gdf, 
print.progress = FALSE, iter = 199)
anova(fit) # treatment effect not adjusted
anova(fit, error = c("treatment:family", "Residuals")) # treatment effect updated (adjusted)






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
pleth.pgls <- procD.pgls(coords ~ Csize, phy = phy, data = gdf, iter = 999)
anova(pleth.pgls)
summary(pleth.pgls)  #similar output

### Working with procD.pgls objects
predict(pleth.pgls)
plot(pleth.pgls, type="regression", reg.type="RegScore", predictor = gdf$Csize)
attributes(pleth.pgls) # Note the PGLS object
attributes(pleth.pgls$PGLS) # PGLS details embedded within PGLS object
pleth.pgls$LM$Pcov # the projection matrix derived from the phylogenetic covariance matrix
pleth.pgls$pgls.fitted # the PGLS fitted values 
pleth.pgls$GM$pgls.fitted # The same fitted values, in a 3D array



cleanEx()
nameEx("read.ply")
### * read.ply

flush(stderr()); flush(stdout())

### Name: read.ply
### Title: Read mesh data (vertices and faces) from ply files
### Aliases: read.ply
### Keywords: IO

### ** Examples

# If the file has no mesh color, or color is undesirable, user can assign this as follows:
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

# A true example is not possible, as digitizing experiences are unqiue, but here is a general set-up
# myShapes <- readShapes("myDigitizingFile") # data from readShapes from StereoMorph 
# myGMdata <- readland.shapes(myShapes) # just reading in the fixed landmarks
# myGMdata <- readland.shapes(myShapes, 
#      nCurvePts = c(10, 15, 5), 
#      continuous.curve = 2) # fixed landmarks plus curve points for three curves, one closed
# myGPA <- gpagen(myGMdata, ProcD = FALSE) # GPA perfomed with minimized bending energy



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

data("plethodon")

Y.gpa <- gpagen(plethodon$land)    #GPA-alignment    
plotTangentSpace(Y.gpa$coords)

preds <- shape.predictor(Y.gpa$coords, x= NULL, Intercept = FALSE, 
pred1 = -0.1, pred2 = 0.1) # PC 1 extremes, sort of
M <- mshape(Y.gpa$coords)
plotRefToTarget(M, preds$pred1)
plotRefToTarget(M, preds[[1]]) # same result
plotRefToTarget(M, preds$pred2)

PCA <- plotTangentSpace(Y.gpa$coords)
PC <- PCA$pc.scores[,1]
preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
plotRefToTarget(M, preds$pred1)
plotRefToTarget(M, preds$pred2)

PC <- PCA$pc.scores[,1:2]
# user-picked spots - can be anything, but it in this case, apparent groups
preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
                        pred1 = c(0.045,-0.02), pred2 = c(-0.025,0.06), pred3 = c(-0.06,-0.04)) 
plotRefToTarget(M, preds$pred1)
plotRefToTarget(M, preds$pred2)
plotRefToTarget(M, preds$pred3)

# allometry example - straight-up allometry

preds <- shape.predictor(Y.gpa$coords, x= log(Y.gpa$Csize), Intercept = TRUE, 
                        predmin = min(log(Y.gpa$Csize)), predmax = max(log(Y.gpa$Csize))) 

plotRefToTarget(M, preds$predmin, mag=3)
plotRefToTarget(M, preds$predmax, mag=3)

# allometry example - using RegScore or PredLine via procD.lm

gdf <- geomorph.data.frame(Y.gpa)
plethAllometry <- procD.lm(coords ~ log(Csize), data=gdf)
allom.plot <- plot(plethAllometry, 
type = "regression", 
predictor = log(gdf$Csize),
reg.type ="RegScore") # make sure to have a predictor 

preds <- shape.predictor(plethAllometry$GM$fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                        predmin = min(allom.plot$RegScore), 
                        predmax = max(allom.plot$RegScore)) 
plotRefToTarget(M, preds$predmin, mag=3)
plotRefToTarget(M, preds$predmax, mag=3)

allom.plot <- plot(plethAllometry, 
type = "regression", 
predictor = log(gdf$Csize),
reg.type ="PredLine")
preds <- shape.predictor(plethAllometry$GM$fitted, x= allom.plot$PredLine, Intercept = FALSE, 
                        predmin = min(allom.plot$PredLine), 
                        predmax = max(allom.plot$PredLine)) 
plotRefToTarget(M, preds$predmin, mag=3)
plotRefToTarget(M, preds$predmax, mag=3)

# using factors via PCA

gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)
pleth <- procD.lm(coords ~ species*site, data=gdf)
PCA <- prcomp(pleth$fitted)
plot(PCA$x, asp=1, pch=19)

means <- unique(round(PCA$x,3))
means # note: suggests 3 PCs useful enough

preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x= PCA$x[,1:3],
                        Intercept = FALSE,
                        pred1 = means[1,1:3], 
                        pred2 = means[2,1:3],
                        pred3 = means[3,1:3], 
                        pred4 = means[4,1:3])                   
plotRefToTarget(M, preds$pred1, mag=2)
plotRefToTarget(M, preds$pred2, mag=2)
plotRefToTarget(M, preds$pred3, mag=2)
plotRefToTarget(M, preds$pred4, mag=2)

# Using a design matrix for factors

X <- pleth$X
X # includes intercept; remove for better functioning 
X <- X[,-1]
symJord <- c(0,1,0) # design for P. Jordani in sympatry
alloJord <- c(0,0,0) # design for P. Jordani in allopatry
preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x = X, Intercept = TRUE, 
                        symJord=symJord, alloJord=alloJord)
plotRefToTarget(M, preds$symJord, mag=2)
plotRefToTarget(M, preds$alloJord, mag=2)

# PLS Example

data(plethShapeFood) 
Y.gpa<-gpagen(plethShapeFood$land)    #GPA-alignment    

# 2B-PLS between head shape and food use data
PLS <-two.b.pls(A1 = plethShapeFood$food, A2 = Y.gpa$coords, iter=999) 
summary(PLS)
plot(PLS)

preds <- shape.predictor(Y.gpa$coords, plethShapeFood$food, Intercept = FALSE,
                        method = "PLS",
                        pred1 = 2, pred2 = -4, pred3 = 2.5) # using PLS plot as a guide
M <- mshape(Y.gpa$coords)
plotRefToTarget(M, preds$pred1, mag=2)
plotRefToTarget(M, preds$pred2, mag=2)
plotRefToTarget(M, preds$pred3, mag=2)




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
plot(PLS)
 
 ### Visualize shape variation using picknplot.shape Because picknplot requires 
 ### user decisions, the following example
 ### is not run (but can be with removal of #).
 ### For detailed options, see the picknplot help file
 # picknplot.shape(plot(PLS))
 




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
