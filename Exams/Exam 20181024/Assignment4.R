# PLOT GP plots gp given function mean, variance, and x values. Plots variance for y if sigmaN is provided
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


library(kernlab)

kernelWrapper = function(sigmaF, l) {
  kernel <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA, n1, n2)
    for (i in 1:n2) {
      r = abs(x1-x2[i])
      K[, i] <-
        sigmaF^2*(1+sqrt(3)*r/l)*exp(-sqrt(3)*r/l)
    }
    return(K)
  }
  class(kernel) <- 'kernel'
  
  return(kernel)
}

# Posterior of gaussian process
posteriorGP = function (X, y, Xstar, sigmaNoise, k){
  kstar = k(X,Xstar)
  L = chol(k(X,X)+sigmaNoise^2*diag(length(X)))
  L = t(L)
  alpha = solve(t(L), solve(L,y))
  fMean = t(kstar) %*% alpha
  v = solve(L,kstar)
  fVar = k(Xstar, Xstar)-t(v)%*%v
  logMarginalLike = -0.5*t(y)%*%alpha-sum(sum(log(L)))-length(y)/2*log(2*pi)
  return (list(fmean=fMean, fvar=fVar, logLike=logMarginalLike))
}

sigmaF_list = list(1, sqrt(5))
graphics.off()
par(mfrow=c(1,2))
for(sigmaF in sigmaF_list){
  l = 0.5
  
  kernel = kernelWrapper(sigmaF, l)
  zGrid = seq(0.01, 1, by=0.01)
  
  k0z = kernel(0, zGrid)
  plot(x=zGrid, y=k0z, type="l", lwd=2, col="blue", main = "K(0,z)")
}

### (2)
lidarData <- read.table('https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/LidarData', 
                        header = T)
sigmaN = 0.05

GPfit <- gausspr(
  lidarData$Distance,
  lidarData$LogRatio,
  kpar = list(sigmaF = 1, l = 1),
  kernel = kernelWrapper,
  var = sigmaN ^ 2
)

# a()
fstar = predict(GPfit, lidarData$Distance)
plotGP(fstar, NULL, lidarData$Distance, lidarData$Distance, lidarData$LogRatio)

# b() c()
# BETTER ALTERNATIVE
gpParam = posteriorGP(scale(lidarData$Distance), scale(lidarData$LogRatio),
                      scale(lidarData$Distance), sigmaN/sd(lidarData$LogRatio)^2, kernelWrapper(sigmaF = 1, l = 1))


graphics.off()
par(mfrow=c(2,1))
for (length in list(1,5)){
  gpParam = posteriorGP(scale(lidarData$Distance), scale(lidarData$LogRatio),
                        scale(lidarData$Distance), sigmaN/sd(lidarData$LogRatio)^2, kernelWrapper(sigmaF = 1, l = length))
  fmean_unscaled = mean(lidarData$LogRatio) + gpParam$fmean*sd(lidarData$LogRatio)
  fvar_unscaled = gpParam$fvar*sd(lidarData$LogRatio)^2
  plotGP(fmean_unscaled, fvar_unscaled, lidarData$Distance, lidarData$Distance, lidarData$LogRatio, sigmaN)
}

# Second alternative, variance from own function with unscaled data and
# fstar from gausspr

sigmaN = 0.05

gpParam = posteriorGP(lidarData$Distance, lidarData$LogRatio,
                      lidarData$Distance, sigmaN, kernelWrapper(sigmaF = 1, l = 1))

plotGP(fstar, gpParam$fvar, lidarData$Distance, lidarData$Distance, lidarData$LogRatio, sigmaN)
