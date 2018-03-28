##Test
ECDFPlot<-function(
      X, ## matrix estimated correlations or list of matrices
      Y, ## matrix reference correlations
      index="all", ##vector
      col.X="red",
      col.Y="black",
      title,
      legend
){ 
  
  if(class(X)=="matrix"){
    if(index[1]!="all"){
      X<-X[index, index]
      Y<-Y[index, index]
    }
    
    X<-X[lower.tri(X)]
    Y<-Y[lower.tri(Y)]
    
    ECDF.X<-ecdf(abs(X))
    ECDF.Y<-ecdf(abs(Y))
    
    plot(ECDF.Y, xlim=c(0, 1), ylim=c(0, 1), col=col.Y, 
    xlab="|Correlation|", ylab="Proportion of |Correlation|", lwd=2,
    main=paste(title), verticals = TRUE, do.p=FALSE, bty="l")
    lines(ECDF.X, col = col.X, lwd=2, verticals=TRUE, do.p=FALSE)
    legend("bottomright", legend=paste(legend), col = c(col.X, col.Y), 
    bty="n", lty=1, lwd=3, cex=0.95)
  }
  
  if(class(X)=="list"){
    if(index[1]!="all"){
      X<-lapply(X, function(x) x[index, index])
      Y<-Y[index, index]
    }
    
    X<-lapply(X, function(x) x[lower.tri(x)])
    Y<-Y[lower.tri(Y)]
    
    ECDF.X<-lapply(X, function(x) ecdf(abs(x)))	
    ECDF.Y<-ecdf(abs(Y))
    
    if(length(col.X)!=length(X)){
      warning("Specified colors are no longer valid.")
      
      col.X<-hcl(h = seq(0, 360, round(360/length(X), 2)), c=45, l=70)[1:length(X)]
    }
    plot(ECDF.Y, xlim=c(0, 1), ylim=c(0, 1), col=col.Y, 
    xlab="|Correlation|", ylab="Proportion of |Correlation|", lwd=2,
    main=paste(title), verticals=TRUE, do.p=FALSE, bty="l")
    for (i in 1:length(X)){
      lines(ECDF.X[[i]], col=col.X[i], lwd=2, verticals=TRUE, do.p=FALSE)
    }
    legend("bottomright", legend=paste(legend), col=c(col.X, col.Y), 
    bty="n", lty=1, lwd=3, cex=0.95)
    
  }
  
}