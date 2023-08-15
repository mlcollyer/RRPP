#' Convert RRPP plots to ggplot objects
#'
#' Function attempts to coerce plot information from an RRPP plot object to an 
#' amenable ggplot object.  
#'
#' This function will attempt to use the plot arguments from an RRPP plot object 
#' to make a ggplot that can be additionally updated, as desired.  Not all plot 
#' characteristics can be converted.  For example, text arguments are not currently
#' passed to \code{\link{ggplot}}, as the \code{\link{text}} function and geom_text
#' arguments do not easily align.  However, one can use text arguments produced by
#' a RRPP plot object and geom_text to augment a ggplot object the way they like.
#' 
#' This function assumes no responsibility for arguments made by \code{\link{ggplot}}.
#' It merely produces a ggplot object that should resemble an RRPP plot default.  Any 
#' augmentation of ggplot objects can be done either by direct intervention of the ggplot 
#' produced or reformatting the initial RRPP plot produced.  One should not expect direct
#' correspondence between R base plot parameters and ggplot parameters.  For example,
#' error bars will generally appear as different widths, without an easy way to control them,
#' changing from one format to the other.
#' 
#' @param object A plot object produced from \code{\link{plot.lm.rrpp}} and type 
#' equals either "PC" or "regression, \code{\link{plot.predict.lm.rrpp}}, or 
#' \code{\link{plot.ordinate}}.  Essentially, any RRPP plot except a series of diagnostic
#' plots should work.
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @examples
#' 
#' ### Linear Model Example
#' 
#' data(Pupfish)
#' fit <- lm.rrpp(coords ~ log(CS) + Sex*Pop, SS.type = "I", 
#' data = Pupfish, print.progress = FALSE) 
#' 
#' # Predictions (holding alternative effects constant)
#' 
#' shapeDF <- expand.grid(Sex = levels(Pupfish$Sex), 
#' Pop = levels(Pupfish$Pop))
#' rownames(shapeDF) <- paste(shapeDF$Sex, shapeDF$Pop, sep = ".")
#' 
#' shapePreds <- predict(fit, shapeDF)
#' summary(shapePreds, PC = TRUE)
#' 
#' # Plot prediction
#' 
#' P <- plot(shapePreds, PC = TRUE, ellipse = TRUE)
#' convert2ggplot(P)
#' 
#' ### Ordination Example
#' 
#' data("PlethMorph")
#' 
#' Y <- as.data.frame(PlethMorph[c("TailLength", "HeadLength", 
#'                                "Snout.eye", "BodyWidth", 
#'                                "Forelimb", "Hindlimb")])
#'Y <- as.matrix(Y)
#'R <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
#'              print.progress = FALSE)$LM$residuals
#'
#' # PCA (on correlation matrix)
#'
#'
#' PCA.ols <- ordinate(R, scale. = TRUE)
#' PCA.ols$rot
#' prcomp(R, scale. = TRUE)$rotation # should be the same
#'
#'
#' PCA.gls <- ordinate(R, scale. = TRUE, 
#'                    transform. = FALSE, 
#'                    Cov = PlethMorph$PhyCov)
#'                
#' P <- plot(PCA.gls)
#' convert2ggplot(P)                   
#' 

convert2ggplot <- function(object){
  x <- y <- NULL
  pa <- object$plot_args
  if(!is.null(object$arrow_args)) {
    arrowed <- TRUE
    ellipses <- FALSE
    ea <- object$arrow_args
  } else if(!is.null(object$ellipse.points)) {
    arrowed <- FALSE
    ellipses <- TRUE
    ea <- object$ellipse.points
  } else {
    arrowed <- FALSE
    ellipses <- FALSE
  }
  
  if(is.null(pa)) 
    stop("Plot object does not appear to be a type that can be converted to ggplot.\n",
         call. = FALSE)
  
  if(arrowed){
    
    xbar <- ybar <- FALSE
    if(!identical(ea$x0, ea$x1)) xbar <- TRUE
    if(!identical(ea$y0, ea$y1)) ybar <- TRUE
    
    xb <- if(xbar) abs(ea$x1 - ea$x0) else NULL
    yb <- if(ybar) abs(ea$y1 - ea$y0) else NULL
    
    df <- data.frame(x = ea$x0, y = ea$y0)
    if(xbar) df$xb <- xb
    if(ybar) df$yb <- yb
    
  }  else if(ellipses) {
    df <- as.data.frame(ea$means)
    colnames(df) <- c("x", "y")
  } else  df <- data.frame(x = pa$x, y = pa$y)
                                
  if(is.null(pa$cex)) pa$cex <- 1
  if(is.null(pa$pch)) pa$pch <- 19
  if(is.null(pa$bg)) pa$bg <- NA
  if(is.null(pa$col)) pa$col <- 1
  if(is.null(pa$lty)) pa$lty <- 1
  if(is.null(pa$lwd)) pa$lwd <- 1
  
  if(arrowed || ellipses) {
    if(is.null(ea$cex)) ea$cex <- 1
    if(is.null(ea$pch)) ea$pch <- 19
    if(is.null(ea$bg)) ea$bg <- NA
    if(is.null(ea$col)) ea$col <- 1
    if(is.null(ea$lty)) ea$lty <- 1
    if(is.null(ea$lwd)) ea$lwd <- 1
    if(!is.null(ea$length) && ea$length < 0.5) ea$length <- 0.9
    if(is.null(ea$length)) ea$length <- 0
    
    if(is.null(pa$xlim)) xLim <- c(min(df$x), max(df$x))
    if(is.null(pa$ylim)) yLim <- c(min(df$y), max(df$y))
  }

  if(ellipses) {
    xLim <- ea$ellipse.points$pc1lim
    yLim <- ea$ellipse.points$pc2lim
  }
  if(!is.null(pa$xlim)) xLim <- pa$xlim
  if(!is.null(pa$ylim)) yLim <- pa$ylim
  
  # 5% span boost for graphics
  xboost <- 0.05 *(xLim[2] - xLim[1])
  yboost <- 0.05 *(yLim[2] - yLim[1])
  
  xLim[1] <- xLim[1] - xboost * abs(xLim[1])
  xLim[2] <- xLim[2] + xboost * abs(xLim[2])
  yLim[1] <- yLim[1] - yboost * abs(yLim[1])
  yLim[2] <- yLim[2] + yboost * abs(yLim[2])

  if(!is.null(pa$asp)){
    g <- ggplot(df, aes(x, y)) + 
      labs(title = pa$main, x = pa$xlab, y = pa$ylab) +
      xlim(xLim) + ylim(yLim) + coord_fixed(ratio = pa$asp)
  } else {
    g <- ggplot(df, aes(x, y)) + 
      labs(title = pa$main, x = pa$xlab, y = pa$ylab) +
      xlim(xLim) + ylim(yLim) 
    
  }

  if(ellipses) {
    ep <- ea$ellP
    k <- dim(ep)[3]
    for(i in 1:k){
      pts <- as.data.frame(ep[,,i])
      colnames(pts) <- c("x", "y")
      g <- g + geom_polygon(data = pts, aes(x, y), linetype = ea$lty, size = ea$lwd,
                         color = ea$col, fill=NA)
    }
  }
  
  if(arrowed){
    check <- match(names(df), c("x", "y", "xb", "yb"))
    if(length(check) == 4) {
      g <- g + geom_errorbar(data = df, aes(xmin = x - xb, 
                                            xmax = x + xb), 
                             size = ea$lwd,
                             linetype = ea$lty,
                             color = ea$col,
                             width = ea$length
      )
      g <- g + geom_errorbar(data = df, aes(ymin = y - yb,
                                            ymax = y + yb), 
                             size = ea$lwd,
                             linetype = ea$lty,
                             color = ea$col,
                             width = ea$length
      )
    }

    
    if(length(check) == 3) {
      if(check[3] == 3)
        g <- g + geom_errorbar(data = df, aes(xmin = x - xb, 
                                              xmax = x + xb),
                               size = ea$lwd,
                               linetype = ea$lty,
                               color = ea$col,
                               width = ea$length
        )
      if(check[3] == 4)
        g <- g + geom_errorbar(data = df, aes(ymin = y - yb, 
                                              ymax = y + yb),
                               size = ea$lwd,
                               linetype = ea$lty,
                               color = ea$col,
                               width = ea$length
        )
    }
  }
  
  pts <- data.frame(x = pa$x, y = pa$y)
  
  g <- g + geom_point(data = pts, aes(x, y), 
                      shape = pa$pch, size = pa$cex * 5, col = pa$col, fill = pa$bg) 
  g
}