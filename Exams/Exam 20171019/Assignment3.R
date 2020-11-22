library("mvtnorm")

kernelMaker = function(sigmaF = 1, l = 0.3) {
  expKernel <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA, n1, n2)
    for (i in 1:n2) {
      K[, i] <- sigmaF ^ 2 * exp(-0.5 * ((x1 - x2[i]) / l) ^ 2)
    }
    return(K)
  }
  # if used with gp library
  class(expKernel) <- 'kernel'
  
  return(expKernel)
}

# samples from the GP prior over a grid of X given a kernelwrapper
sample_and_plot_F = function(xgrid, kernel, sigmaF, l, number_of_F) {
  cov_matrix = kernel(sigmaF, l)(xgrid, xgrid)
  f_samples = rmvnorm(number_of_F, sigma = cov_matrix)
  
  plot.new()
  box()
  plot.window(xlim = c(min(xgrid), max(xgrid)), ylim = c(min(f_samples), max(f_samples)))
  
  axis(1)
  axis(2)
  title(
    main = sprintf("%i samples of F for sigmaF=%g and l=%g", number_of_F, sigmaF, l),
    xlab = "x",
    ylab = "f"
  )
  colors = rainbow(number_of_F)
  
  for (i in 1:number_of_F) {
    lines(
      x = xgrid,
      y = f_samples[i, ],
      col = colors[i],
      lwd = 2
    )
  }
}

# Assignment a)

xgrid = seq(-1, 1, by = 0.1)


# Sample 5 different functions from each covariance matrix for l=0.2 and l=1
number_of_F = 5
graphics.off()
par(mfrow = c(2, 1))
sample_and_plot_F(xgrid, kernelMaker, 1, 0.2, number_of_F)
sample_and_plot_F(xgrid, kernelMaker, 1, 1, number_of_F)

for (length in list(0.2, 1)) {
  corr1 = kernelMaker(1, length)(0, 0.1)
  corr2 = kernelMaker(1, length)(0, 0.5)
  print(sprintf(
    "L=%g, corr(f(0),f(0.1))=%g  corr(f(0),f(0.5))=%g",
    length,
    corr1,
    corr2
  ))
}


# Assignment b)
posteriorGP = function (X, y, Xstar, sigmaNoise, k) {
  kstar = k(X, Xstar)
  L = chol(k(X, X) + sigmaNoise ^ 2 * diag(length(X)))
  L = t(L)
  alpha = solve(t(L), solve(L, y))
  fMean = t(kstar) %*% alpha
  v = solve(L, kstar)
  fVar = k(Xstar, Xstar) - t(v) %*% v
  logMarginalLike = -0.5 * t(y) %*% alpha - sum(sum(log(L))) - length(y) /
    2 * log(2 * pi)
  return (list(
    fmean = fMean,
    fvar = fVar,
    logLike = logMarginalLike
  ))
}

# PLOT GP plots gp given function mean, variance, and x values. Plots variance for y if sigmaN is provided
plotGP = function(fmean,
                  fvar = NULL,
                  Xstar,
                  X,
                  y,
                  sigmaN = 0,
                  xlab = "time",
                  ylab = "temp") {
  subheading = sprintf("sigmaN = %g  %i samples  %i observed",
                       sigmaNoise,
                       length(Xstar),
                       length(X))
  legendPosition = "bottomright"
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
load("GPdata.RData")

graphics.off()
par(mfrow = c(2, 1))
for (length in list(0.2, 1)) {
  sigmaNoise = 0.2
  Xstar = seq(min(x), max(x), by = 0.001)
  var = posteriorGP(x, y, Xstar, sigmaNoise, kernelMaker(sigmaF = 1, l =
                                                           0.2))$fvar
  GPfit <- gausspr(
    x,
    y,
    kpar = list(sigmaF = 20, l = length),
    kernel = kernelMaker,
    var = sigmaNoise ^ 2
  )
  fposterior = predict(GPfit, Xstar)
  plotGP(fposterior, var, Xstar, x, y, sigmaNoise)
}
