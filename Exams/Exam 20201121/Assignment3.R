kernelMaker = function(sigmaF=1,l=0.3){
  
  expKernel <- function(x1,x2){
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaF^2*exp(- 0.5*( (x1-x2[i])/l )^2 )
    }
    return(K)
  }
  
  return(expKernel)
}





plotGP = function(fmean,
                  fvar = NULL,
                  Xstar,
                  X,
                  y,
                  sigmaN = 0,
                  xlab = "time",
                  ylab = "temp") {
  subheading = sprintf("%i samples  %i observed",
                       length(Xstar),
                       length(X))
  legendPosition = "topright"
  # IF VARIANCE GIVEN
  if (!is.null(fvar)) {
    upperConf = (fmean + sqrt(diag(fvar)) * 1.96)[, 1]
    lowerConf = (fmean - sqrt(diag(fvar)) * 1.96)[, 1]
    
    #Conf for y
    legendText = c("data", "post mean", "95% intervals for f")
    legendCols = c("green", "blue",  rgb(0, 0, 0, 0.3))
    if (!sigmaN == 0) {
      legendText = c("data",
                     "post mean",
                     "95% intervals for f",
                     "95% intervals for y")
      legendCols = c("green", "blue",  rgb(0, 0, 0, 0.3), rgb(0, 0, 0.8, 0.3))
      Y_upperConf = (fmean + sqrt(diag(fvar) + sigmaN ^ 2) * 1.96)[, 1]
      Y_lowerConf = (fmean - sqrt(diag(fvar) + sigmaN ^ 2) * 1.96)[, 1]
    }
    plot(
      x = Xstar,
      y = fmean,
      col = "blue",
      type = "l",
      ylim = c(min(lowerConf), max(upperConf)),
      lwd = 2,
      sub = subheading,
      main = "GP",
      ylab = ylab,
      xlab = xlab
    ) +
      points(x = X, y = y, col = "green")
    legend(
      legendPosition,
      inset = 0.02,
      legend = legendText,
      col = legendCols,
      pch = c('o', NA, NA, NA),
      lty = c(NA, 1, 1, 1),
      lwd = 2,
      cex = 0.55
    )
    
    polygon(
      x = c(rev(Xstar), Xstar),
      y = c(rev(upperConf), lowerConf),
      col = rgb(0, 0, 0, 0.3)
    )
    if (!sigmaN == 0) {
      polygon(
        x = c(rev(Xstar), Xstar),
        y = c(rev(Y_upperConf), Y_lowerConf),
        col = rgb(0, 0, 0.8, 0.3)
      )
    }
  }
  #IF NO VARIANCE GIVEN
  else {
    plot(
      x = Xstar,
      y = fmean,
      col = "blue",
      type = "l",
      lwd = 2,
      sub = subheading,
      main = "GP",
      ylab = ylab,
      xlab = xlab
    ) +
      points(x = X, y = y, col = "green")
    legend(
      legendPosition,
      inset = 0.02,
      legend = c("data", "post mean"),
      col = c("green", "blue"),
      pch = c('o', NA, NA, NA),
      lty = c(NA, 1, 1, 1),
      lwd = 2,
      cex = 0.55
    )
    
  }
  
}
posteriorGP = function (X, y, Xstar, sigmaNoise, k){
  kstar = k(X,Xstar)
  L = chol(k(X,X)+sigmaNoise^2*diag(length(X)))
  L = t(L)
  alpha = solve(t(L), solve(L,y))
  fMean = t(kstar) %*% alpha
  v = solve(L,kstar)
  fVar = k(Xstar, Xstar)-t(v)%*%v
  logMarginalLike = -0.5*(t(y)%*%alpha)-sum(diag(L))-length(y)/2*log(2*pi)
  return (list(fmean=fMean, fvar=fVar, logLike=logMarginalLike))
}

graphics.off()
par(mfrow=c(1,1))
X = c(-1.0, -0.6, -0.2, 0.4, 0.8)
y = c(0.768, -0.044, -0.940, 0.719, -0.664)

sigmaNoise = 0
Xstar = seq(-1,1,by=0.01)

gpParam = posteriorGP(X, y, Xstar, sigmaNoise, kernelMaker(sigmaF=1, l=0.3))
plotGP(gpParam$fmean, gpParam$fvar, Xstar, X, y)


# (1)

## The covariance from x=0 to all other xs is found in the variance for the corresponding column of
# x=0 
plot(x=Xstar, y=gpParam$fvar[which(Xstar==0),], type="l", xlab="x", ylab="covariance", main="Covariance from x'=0 to x", col="blue", lwd=2)
# as we can see, maximum covariance is found around 0, but slightly shifted.
# this is because the points in the data influences the covariance by making x values with similiar y values more correlated
# for instance, the training point at 0.8 is ~ 0.3 away from the training point at -0.2 in y value. This means that the
# corresponding x values for theese points are slightly correlated, explainging the curvature upwards towards x=1 in the covariance plot


# (2)
plot(x=Xstar, y=gpParam$fvar[which(Xstar==0),], ylim=c(min(y),max(y)),
     type="l", xlab="x", ylab="covariance", main="Posterior Covariance from x'=0 to x", col="blue", lwd=2)
points(x=X, y=y, col="green", lwd=3)

# At all the plotted training points, the covariance is 0. This because the the estimated f value at 0
# has no correlation with the x values that has a training point, as this covariance only depends on the
# training point at that x value

# (3)
# because as we mentioned earlier, the covariance will be influenced by the training points, and thus training
# points close to 0 will have correlation with other points where the training points for 0 and that x value
# are similiar


# (4) the prior for the covariance of x and x' is given by the exponential kernel
varPrior = kernelMaker(sigmaF=1, l=0.3)(Xstar, 0)
plot(x=Xstar, y=varPrior, xlab="x", ylab="covariance", type="l", main="Prior Covariance from x'=0 to x", col="blue", lwd=2)
