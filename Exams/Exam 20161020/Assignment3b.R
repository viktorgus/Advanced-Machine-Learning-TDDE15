# PLOT GP plots gp given function mean, variance, and x values. Plots variance for y if sigmaN is provided
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
  legendPosition = "bottomright"
  # IF VARIANCE GIVEN
  if (!is.null(fvar)) {
    upperConf = (fmean + sqrt((fvar)) * 1.96)
    lowerConf = (fmean - sqrt((fvar)) * 1.96)
    
    #Conf for y
    legendText = c("data", "post mean", "95% intervals for f")
    legendCols = c("green", "blue",  rgb(0, 0, 0.6, 0.05))
    if (!sigmaN == 0) {
      legendText = c("data",
                     "post mean",
                     "95% intervals for f",
                     "95% intervals for y")
      legendCols = c("green", "blue",  rgb(0, 0, 0.6, 0.05), rgb(0, 0, 0.1, 0.05))
      Y_upperConf = (fmean + sqrt((fvar) + sigmaN ^ 2) * 1.96)
      Y_lowerConf = (fmean - sqrt((fvar) + sigmaN ^ 2) * 1.96)
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
      points(x = X, y = y, col = "green", lwd=2)
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
      col = rgb(0, 0, 0.6, 0.05)
    )
    if (!sigmaN == 0) {
      polygon(
        x = c(rev(Xstar), Xstar),
        y = c(rev(Y_upperConf), Y_lowerConf),
        col = rgb(0, 0, 0.1, 0.05)
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



K1 = function(sigmaF, l) {
  k1 = function(x,y) {
    r = abs(x-y)
    sigmaF^2 * exp(-0.5 * (r / l)^2)
  }
  class(k1) = "kernel"
  return(k1)
}

K3 = function(sigmaF, l) {
  k3 = function(x, y) {
    r = abs(x-y)
    sigmaF^2 * (1 + sqrt(3) * r / l) * exp(-sqrt(3) * r / l)
  }
  class(k3) = "kernel"
  return(k3)
}


getVar = function (k, x, xstar, sigmaN){
  n = length(x)
  I = diag(n)
  
  Kss = kernelMatrix(kernel = k, x = xstar, y = xstar)
  Ksx = kernelMatrix(kernel = k, x = xstar, y = x)
  Kxs = t(Ksx)
  Kxx = kernelMatrix(kernel = k, x = x, y = x)
  V = Kss - Ksx %*% solve( Kxx + sigmaN^2 * I ) %*% Kxs
  
  return(diag(V))
}

load("GPdata.RData")

sigmaN = 0.5
sigmaF = 1
l = 1
xStar = seq(0,1,length=200)

graphics.off()
par(mfrow=c(1,2))
for (kern in list(K1, K3)){
  GPfit <- gausspr(
    x,
    y,
    kpar = list(sigmaF = sigmaF, l = l),
    kernel = kern,
    var = sigmaN ^ 2
  )
  fvar = getVar(K1(sigmaF, l), x, xStar, sigmaN)
  fmean = predict(GPfit, xStar)[,1]
  
  plotGP(fmean, fvar, xStar, x, y, sigmaN, "x", "y")
}
