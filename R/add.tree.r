#' Plot tool to add trees to ordination plots
#'
#' Function adds a tree based on a description of edges to an existing
#' plot made from an ordinate object.
#' 
#' With some \code{\link{ordinate}} plots, it might be desirable to add a tree connecting points
#' in a prescribed way, which would be tedious using \code{\link{points}} or \code{\link{lines}}.
#' This function will project a tree into a plot with class, \code{\link{plot.ordinate}}.  New data 
#' can olso be projected into the plot (but it is assumed that such data will not require transformation).
#' Using an edges matrix, this function will systematically connect plot points with lines,
#' based on a numeric identification of points.  This funcrion might only work after determining the order
#' of points from the \code{\link{plot.ordinate}} object.
#' 
#' Generally, if new data are used, this function will only be employed after the points
#' in a \code{\link{plot.ordinate}} object have been summarized, somehow.  For example, if one wants
#' to project a phylogenetic tree into a principal component (PC) plot, the ancestral states (nodes of tree) would 
#' have to be first estimated with respect to the PC scores that were generated.  Thsi function
#' will not convert ancestral states estimated from original data into PC scores.
#' 
#' @param OP An object with class \code{\link{plot.ordinate}}.
#' @param newdata An optional matrix or data frame with two columns only, corresponding to the x
#' and y values in the plot.  It is assumed and data transformations have been performed, already.
#' @param edges A matrix or data frame with two columns, with the first value in each row indicating
#' the point to start and the second value the point to end an edge connection.
#' @param edge.col A single value or vector equal to the number of edges for edge colors.
#' @param edge.lty A single value or vector equal to the number of edges for edge line type
#' @param edge.lwd A single value or vector equal to the number of edges for edge line weight.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{lines}} and \code{points} 
#' @examples 
#' data("PlethMorph")
#' Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
#' "Snout.eye", "BodyWidth", "Forelimb", "Hindlimb")])
#' Y <- as.matrix(Y)
#' R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
#' iter = 0, print.progress = FALSE)$LM$residuals
#' 
#' PCA <- ordinate(R, scale. = TRUE)
#' pc.plot <- plot(PCA, pch = 19, col = "blue")
#' 
#' # To run the following, do so without the # tags
#' # (this example uses a function that RRPP does not import)
#' # library(phytools)
#' # ancs <- sapply(1:2, function(j) 
#' # fastAnc(PlethMorph$tree, pc.plot$points[,j]))
#' # add.tree(pc.plot, newdata = ancs, edges = PlethMorph$tree$edge)

add.tree <- function(OP, newdata = NULL, edges, 
                     edge.col = 1, edge.lty = 1, edge.lwd = 1) {
  
  if(!inherits(OP, c("plot.ordinate")))
    stop("\nThe plot must be an object with class ordinate.\n",
         call. = FALSE)
  
  if(!is.null(newdata)) {
    
    if(!inherits(newdata, c("matrix", "data.frame")))
      stop("\nnewdata is not an object with useable data, such as class matrix or data.frame.\n",
           call. = FALSE)
    
    if(NCOL(newdata) != 2)
      stop("\nYour newnewdata can only have two variables and match the variables for the points in the plot.\n",
           call. = FALSE)
    
    z <- rbind(OP$points, newdata)
    
  } else z <- OP$points
  
  if(!inherits(edges, c("matrix", "data.frame")))
    stop("\nedges is not an object with useable data, such as class matrix or data.frame.\n",
         call. = FALSE)
  
  if(NCOL(edges) != 2)
    stop("\nThe edges must be a matrix with two columns.\n",
         call. = FALSE)
  
  edges <- as.matrix(edges)
  
  if(max(edges) > NROW(z))
    stop("\nThe edges refer to more points than are possible in the plot.\n",
         call. = FALSE)
  
  for (i in 1:NROW(edges)) {
    pts <- z[edges[i,], ]
    points(pts, type = "l", col = edge.col, 
            lwd = edge.lwd, lty = edge.lty)
  }
  
  do.call(points, OP$plot.args)
  
}
  