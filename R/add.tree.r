#' Plot tool to add phylogenetic trees to ordination plots
#'
#' Function adds a tree based on a description of edges from a class phylo object 
#' to an existing plot made from an ordinate object.
#' 
#' With some \code{\link{ordinate}} plots, it might be desirable to add a tree 
#' connecting points in a prescribed way, which would be tedious using 
#' \code{\link{points}} or \code{\link{lines}}.  This function will project a 
#' tree from an object of class phylo into a plot with class, 
#' \code{\link{plot.ordinate}}.  Using an edges matrix from a phylo object, 
#' this function will systematically connect plot points with lines that pass 
#' through estimated ancestral character points in the same plot space.  
#' Ancestral states are estimated assuming a Brownian motion model 
#' of evolutionary divergence.

#' @param OP An object with class \code{\link{plot.ordinate}}.
#' @param tree An object of class phylo.
#' @param edge.col A single value or vector equal to the number of edges for edge colors.
#' @param edge.lty A single value or vector equal to the number of edges for edge line type
#' @param edge.lwd A single value or vector equal to the number of edges for edge line weight.
#' @param anc.pts A logical value for whether to add points for ancestral values.
#' @param return.ancs A logical value for whether ancestral values should be printed.
#' @param ... Arguments passed onto \code{\link{points}}, used only for ancestral points.
#' @keywords graphics
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{lines}} and \code{\link{points}}
#' @examples
#' 
#' # Examples use residuals from a regression of salamander morphological 
#' # traits against body size (snout to vent length, SVL).
#' # Observations are species means and a phylogenetic covariance matrix
#' # describes the relatedness among observations.
#'
#' data("PlethMorph")
#' Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
#' "Snout.eye", "BodyWidth", 
#' "Forelimb", "Hindlimb")])
#' Y <- as.matrix(Y)
#' R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
#' iter = 0, print.progress = FALSE)$LM$residuals
#' 
#' PCA <- ordinate(R, scale. = TRUE)
#' pc.plot <- plot(PCA, pch = 19, col = "blue")
#' 
#' add.tree(pc.plot, tree = PlethMorph$tree, anc.pts = TRUE, 
#' pch = 19, cex = 0.5, col = "red")
#' 
add.tree <- function(OP, tree, 
                     edge.col = 1, edge.lty = 1, edge.lwd = 1,
                     anc.pts = FALSE, return.ancs = FALSE, ...) {
  
  if(!inherits(OP, c("plot.ordinate", "plot.kcomp")))
    stop("\nThe plot must be a plot object obtained from ordinate or kcomp.\n",
         call. = FALSE)
  
  if(!inherits(tree, c("phylo")))
    stop("\ntree is not a class phylo object.\n",
         call. = FALSE)

  pts <- as.matrix(OP$points)
  ancs <- anc.BM(tree, pts)
  
  if(is.null(rownames(pts)))
    stop("Plot points do not have associated taxa names.\n", 
         call. = FALSE)
  
  ind <- match(tree$tip.label, rownames(pts))
  
  if(any(is.na(ind)))
    stop("\nRow names of data and tree tip names do not match.\n",
         call. = FALSE)
  
  pts <- pts[ind, ]
  z <- rbind(pts, ancs)
  
  edges <- as.matrix(tree$edge)
  
  for (i in 1:NROW(edges)) {
    pts <- z[edges[i,], ]
    points(pts, type = "l", col = edge.col, 
            lwd = edge.lwd, lty = edge.lty)
  }
  
  if(anc.pts) points(ancs, ...)
  
  do.call(points, OP$plot_args)
  
  if(return.ancs) return(ancs)
  
}
  