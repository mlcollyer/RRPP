#' Landmarks on pupfish
#'
#' @name Pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprinodon pecosensis body shapes, with 
#' indication of Sex and
#' Population from which fish were sampled (Marsh or Sinkhole).
#' @details These data were previously aligned with GPA.  Centroid size (CS) 
#' is also provided.  
#' See the \pkg{geomorph} package for details.
#' 
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for 
#' analysis of phenotypic
#' change for phenotypes described by high-dimensional data. Heredity. 113: 
#' doi:10.1038/hdy.2014.75.
NULL

#' Landmarks on pupfish heads
#'
#' @name PupfishHeads
#' @docType data
#' @author Michael Collyer
#' @description Landmark data from Cyprinodon pecosensis head shapes, with 
#' variables for 
#' sex, month and year sampled, locality, head size, and coordinates of 
#' landmarks for head shape,
#' per specimen.  These data are a subset of a larger data set.
#' @details The variable, "coords", are data that were previously aligned
#' with GPA.  The variable, "headSize", is the Centroid size of each vector 
#' of coordinates.
#' See the \pkg{geomorph} package for details.
#' @references Gilbert, M.C. 2016. Impacts of habitat fragmentation on the 
#' cranial morphology of a 
#' threatened desert fish (Cyprinodon pecosensis). Masters Thesis, 
#' Western Kentucky University.
NULL

#' Experimental metabolic rate data from pupfish
#'
#' @name PupfishMR
#' @docType data
#' @author Michael Collyer
#' @description Simulated metabolic rate data from experimental pupfish
#' @details These data are simulated data meant to mimic a subset of the 
#' data from Collyer & Stockwell (2004).  Variables include experimental treatment 
#' (Control versus artificially infected), mass in grams, standard length in mm,
#' and metabolic rate (oxygen consumption, mg O2 per L per hour).  These data do not have
#' a temporal component like the original data, and do not include aberrant or missing
#' observations.  These data are suitable data for demonstrating linear regression
#' and analysis of covariance 
#' (ANCOVA), used in the book, Permutational Biometry (Collyer & Adams, 2025).
#' 
#' @references Collyer, M. L., & Stockwell, C. A. (2004). Experimental 
#' evidence for costs of parasitism for a threatened species, 
#' White Sands pupfish (Cyprinodon tularosa). Journal of 
#' Animal Ecology, 73(5), 821-830.
#' 
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
#' 
NULL

#' Plethodon comparative morphological data 
#'
#' @name PlethMorph
#' @docType data
#' @author Michael Collyer and Dean Adams
#' @keywords datasets
#' @description Data for 37 species of plethodontid salamanders.  
#' Variables include snout to vent length
#' (SVL) as species size, tail length, head length, snout to eye length, 
#' body width, forelimb length,
#' and hind limb length, all measured in mm.  A grouping variable is also 
#' included for functional guild size.  A variable for species names is also 
#' included.
#' The data set also includes a phylogenetic covariance matrix based on a 
#' Brownian model of evolution, to assist in 
#' generalized least squares (GLS) estimation.
#' @details The covariance matrix was estimated with the vcv.phylo function 
#' of the R package, ape, based on the tree
#' described in Adams and Collyer (2018).
#' @references Adams, D.C and Collyer, M.L. 2018. Multivariate phylogenetic 
#' anova: group-clade aggregation, biological 
#' challenges, and a refined permutation procedure. Evolution, 72: 1204-1215.
NULL


#' Plethodon comparative SVL data 
#'
#' @name PlethCinHoff
#' @docType data
#' @author Dean Adams
#' @keywords datasets
#' @description Snout-to-vent-length (SVL) data for two species of salamanders 
#' @details The data set is a subset of data from two species of salamanders 
#' (Plethodon cinereus and 
#' P. hoffmani) in eastern North America (based on Adams and Rohlf, 2000; 
#' Adams, 2000), collected from the Smithsonian Museum (Washington DC, USA).  Data 
#' include species, locality (whether in sympatry or allopatry), and SVL in mm.
#' These data are suitable data for demonstrating single factor and
#' factorial analysis of variance (ANOVA).
#' (ANCOVA), used in the book, Permutational Biometry (Collyer & Adams, 2025).
#' @references Adams, D. C. (2000). Divergence of trophic morphology and resource use
#' among populations of Plethodon cinereus and P. hoffmani in Pennsylvania.
#' In Bruce, R. C., Jaeger, R. G., and Houck, L. D., editors, The Biology of
#' Plethodontid Salamanders, pages 383–394. Springer US.
#' @references Adams, D. C. and Rohlf, F. J. (2000). Ecological character displacement in
#' Plethodon: Biomechanical differences found from a geometric morphometric
#' study. Proceedings of the National Academy of Sciences, 97:4106–4111.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Reptile Abundance Data in Western Australia
#'
#' @name ReptAbund
#' @docType data
#' @keywords datasets
#' @description Reptile abundance data from several sites in Western Australia,
#' originally collected to examine recovery after a large fire. 
#' @details The data set is from a study examining reptile abundance from
#' several sites in Western Australia (Davis and Doherty, 2015). The original 
#' study examined the recovery of reptile communities in ten remnant
#' localities following a large fire in 2009. Half of localities were in the burn zone 
#' and the remaining half were unaffected by the fire. Data were sampled over five years.
#' These data are suitable data for demonstrating mixed linear models, used in the book, 
#' Permutational Biometry (Collyer & Adams, 2025).
#' @references Davis, R. A. and Doherty, T. S. (2015). Rapid recovery of an urban remnant
#' reptile community following summer wildfire. PLOS ONE, 10(5):e0127925.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Mammal Home Range Data
#'
#' @name mammalHR
#' @docType data
#' @keywords datasets
#' @description Home range and body mass data for 49 species of mammals.
#' @details The data set is from a study of mammal species home range size 
#' (squared kilometers) and its covariation with body mass or whether mammals
#' are carnivores or ungulates (Garland et al., 1993).  the data set also includes a
#' covariance matrix based on phylogenetic relatedness among species, to account for the
#' non-independence of observations.
#' These data are suitable data for demonstrating generalized least squares (GLS)
#' estimation, used in the book, 
#' Permutational Biometry (Collyer & Adams, 2025).
#' @references Garland, T., Dickerman, A. W., Janis, C. M., and Jones, J. A. (1993). 
#' Phylogenetic analysis of covariance by computer simulation. Systematic Biology,
#' 42:265–292.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Phylogenetic Data for the Mammal Home Range Data
#'
#' @name mammalHRphy
#' @docType data
#' @keywords datasets
#' @description Phylogenetic tree data that can be used along with 
#' \code{\link{mammalHR}}.
#' @details A class \code{phylo} object containing data for the phylogenetic tree
#' that pertains to the \code{\link{mammalHR}} data set.  Although these data are
#' not required for the generalized least squares (GLS)
#' estimation example, used in the book, 
#' Permutational Biometry (Collyer & Adams, 2025), they provide additional 
#' graphical support for the example.

#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Simulated motion paths
#'
#' @name motionpaths
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for 
#' the analysis of phenotypic
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @keywords datasets
NULL

#' Simulated Fish Data for Measurement Error Analysis
#'
#' @name fishy
#' @docType data
#' @author Michael Collyer
#' @references Collyer, M.L, and D.C. Adams. 2024. 
#' Interrogating random and systematic measurement 
#' error in morphometric data. Evolutionary Biology. 
#' In press.
#' @details Data as simulated in Collyer and Adams (2024),
#' resembling a fish shape, comprising Procrustes coordinates
#' for 11 anatomical landmarks.  Data represent 120 
#' configurations for 60 subjects, each with two replicates of
#' measurement.  The 60 subjects represent 20 subjects each 
#' from three groups.
#' @keywords datasets
NULL

#' Phytoplankton Biomass and Eutrophication
#'
#' @name phytoplankton
#' @docType data
#' @keywords datasets
#' @description A data set containing variables for phytoplankton biomass,
#' phosphorous levels, trophic status, and locations (coordinates) of 70 temperate lakes.
#' @details The data in this data set are samples from a larger study 
#' that compared phytoplantkton biomass
#' from temperate lakes in Alberta, Canada, with different levels of eutrophication
#' (Loewen et al., 2020).  These 
#' data are suitable for considering the spatial non-independence of observations in 
#' generalized least squares (GLS) estimation, for the example in the book,
#' Permutational Biometry (Collyer & Adams, 2025).
#' @references Loewen, C. J. G., Wyatt, F. R., Mortimer, C. A., Vinebrooke, R. D., and
#' Zurawell, R. W. (2020). Multiscale drivers of phytoplankton communities
#' in north‐temperate lakes. Ecological Applications, 30(5):e02102.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Flower Production and Temperature
#'
#' @name flower
#' @docType data
#' @keywords datasets
#' @description A data set containing variables for flower production and
#' temperature, averaged yearly for 30 years, from a dry tropical forest.
#' @details The data in this data set are samples from a larger study 
#' that compared flower production between forest types over 30 years, 
#' based on climate variables.  The data in this set are from a dry tropical forest,
#' include variables for annual flower production, temperature, and year. These 
#' data are suitable for considering the temporal non-independence of observations in 
#' generalized least squares (GLS) estimation, for the example in the book,
#' Permutational Biometry (Collyer & Adams, 2025).
#' @references Pau, S., Wolkovich, E. M., Cook, B. I., Nytch, C. J., Regetz, J., Zimmerman,
#' J. K., and Joseph Wright, S. (2013). Clouds and temperature drive dynamic
#' changes in tropical flower production. Nature Climate Change, 3:838–842.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Snake Head Size and Regional Variation
#'
#' @name snakeHS
#' @docType data
#' @keywords datasets
#' @description A data set containing variables for snake head size, 
#' snout-to-vent-length, sex, and geographic information.
#' 
#' @details These data were collected from prairie rattlesnakes (Crotalus viridis viridis),
#' from several geographic locations in the upper Midwest of the United States (Smith and 
#' Collyer, 2008).  Variables include a spatial measure of head size (HS), 
#' snout-to-vent-length (SVL, cm), sex, and the site and region where 107
#' snakes were captured.  These 
#' data are suitable for demonstrating model selection  in the book,
#' Permutational Biometry (Collyer & Adams, 2025).
#' @references Smith, M. T. and Collyer, M. L. (2008). The biology of rattlesnakes, chapter
#' Regional variation and sexual dimorphism in head form of the prairie
#' rattlesnale (Crotalus viridid viridis): Comparisons using new analytical techniques
#' and collection methods, pages 79–90. Loma Linda University Press.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL

#' Podarcis Wall Lizards from the Mediterranean Basin 
#'
#' @name podarcis
#' @docType data
#' @keywords datasets
#' @description A data set of body size for males and females form different lizard 
#' species in different geographic locations.
#' 
#' @details These data are a subset of the data collected originally by 
#' Kaliontzopoulou et al.(2015).  Variables include snout-to-vent-length (SVL),
#' sex, species, and habitat preference.  These 
#' data are suitable for comparing effect sizes in the book,
#' Permutational Biometry (Collyer & Adams, 2025).
#' @references Kaliontzopoulou, A., Carretero, M. A., and Adams, D. C. (2015). Ecomorphological
#' variation in male and female wall lizards and the macroevolution
#' of sexual dimorphism in relation to habitat use. Journal of Evolutionary
#' Biology, 28:80–94.
#' @references Collyer, M. L., & Adams D. C. (2025). Permutational Biometry: Volume
#' 1, Univariate Data. Iowa State University Digital Press.  NOT YET PUBLISHED.
NULL