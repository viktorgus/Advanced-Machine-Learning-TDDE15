library(kernlab)

par(mfrow = c(1, 1))

kernelWrapPeriod = function(sigmaF, l1, l2, d) {
  expKernel <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA, n1, n2)
    for (i in 1:n2) {
      K[, i] <-
        sigmaF ^ 2 * exp(-2 * (sin(pi * abs(x1 - x2[i]) / d) / l1) ^ 2) *
        exp(-0.5 * (abs(x1 - x2[i]) / l2) ^ 2)
    }
    return(K)
  }
  class(expKernel) <- 'kernel'
  
  return(expKernel)
}

kernelWrap = function(sigmaF, l) {
  expKernel <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA, n1, n2)
    for (i in 1:n2) {
      K[, i] <- sigmaF ^ 2 * exp(-0.5 * ((x1 - x2[i]) / l) ^ 2)
    }
    return(K)
  }
  class(expKernel) <- 'kernel'
  
  return(expKernel)
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


getGPvar = function (X, y, Xstar, sigmaNoise, k) {
  y = scale(y)
  X = scale(X)
  Xstar = scale(Xstar)
  kstar = k(X, Xstar)
  L = chol(k(X, X) + sigmaNoise ^ 2 * diag(length(X)))
  L = t(L)
  v = solve(L, kstar)
  fVar = k(Xstar, Xstar) - t(v) %*% v
  return (fVar)
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
  logMarginalLike = -0.5*(t(y)%*%alpha)-sum(diag(L))-length(y)/2*log(2*pi)
  return (list(fmean=fMean, fvar=fVar, logLike=logMarginalLike))
}

# Posterior of gaussian process
loglike = function (params, X, y, Xstar, sigmaNoise){
  sigmaF = params[1]
  l = params[2]
  k = kernelWrap(sigmaF, l)
  kstar = k(X,Xstar)
  L = chol(k(X,X)+sigmaNoise^2*diag(length(X)))
  L = t(L)
  alpha = solve(t(L), solve(L,y))
  logMarginalLike = -0.5*(t(y)%*%alpha)-sum(diag(L))-length(y)/2*log(2*pi)
  return (logMarginalLike)
}


# preprocessing

data = read.csv(
  "https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/
Code/TempTullinge.csv",
  header = TRUE,
  sep = ";"
)

dates = as.Date(data$date, "%d/%m/%y")
rel_days = sapply(dates, function(x) {
  as.numeric(x - dates[1]) %% 365 + 1
})
abs_days = sapply(dates, function(x) {
  as.numeric(x - dates[1]) + 1
})
data$day = rel_days
data$time = abs_days

small.data = c()
small.data$temp = data$temp[seq(1, length(data$temp), by = 5)]
small.data$day = rel_days[seq(1, length(rel_days), by = 5)]
small.data$time = abs_days[seq(1, length(abs_days), by = 5)]


polyFit <-
  lm(scale(small.data$temp) ~  scale(small.data$time) + I(scale(small.data$time) ^ 2))

unscaled_polyFit <-
  lm(small.data$temp ~  small.data$time + I(small.data$time ^ 2))

sigmaNoise = sd(polyFit$residuals)
sigmaNoise_unscaled = sd(unscaled_polyFit$residuals)



optimal = optim(par = c(1,0.1),
      fn = loglike, X=scale(small.data$time),y=scale(small.data$temp),Xstar=scale(small.data$time),sigmaNoise=sigmaNoise,
      method="L-BFGS-B",
      lower = c(.Machine$double.eps, .Machine$double.eps),
      control=list(fnscale=-1))

gpParam = posteriorGP(scale(small.data$time), scale(small.data$temp), scale(small.data$time), sigmaNoise, kernelWrap(sigmaF=optimal$par[1],l=optimal$par[2]))
fmean_unscaled = mean(small.data$temp) + gpParam$fmean*sd(small.data$temp)
fvar_unscaled = gpParam$fvar*sd(small.data$temp)^2
plotGP(fmean_unscaled, fvar_unscaled, small.data$time, small.data$time, small.data$temp, sigmaNoise_unscaled)

