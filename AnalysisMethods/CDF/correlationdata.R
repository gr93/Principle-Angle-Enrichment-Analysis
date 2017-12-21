correlationPlot<-function(
      true, ##matrix 
      est,
      plot.genes=sample(1:dim(true)[1], 18), ## index of which genes to plot
      boxes=TRUE, ## draw boxes FALSE/TRUE
      title, 
      line=-1)
{
  
  if(is.matrix(true)==FALSE) stop("true needs to be a matrix object.")
  if(is.matrix(est)==FALSE) stop("est needs to be a matrix object.")
  
  if(length(plot.genes)!=18 & boxes==TRUE){
    warning("No boxes can be drawn when not exactly 18 genes are specified.")
    boxes=FALSE
  }
  
  est.true<-mashUp(true, est, plot.genes)
  ## lower triangle estimated values, upper triangle true values
  
  corrplot(est.true, method = c("color"), order = "original", tl.col = "black")
  ## draw object
  mtext(paste(title), side=3, line=line, cex = 1)
  if(boxes){
    rect(0.5, 12.5, 6.5, 18.5, border = "black", lwd = 2)
    rect(6.5, 6.5, 12.5, 12.5, border="black", lwd=2)
    rect(12.5, 0.5, 18.5, 6.5, border="black", lwd=2)
  }
}