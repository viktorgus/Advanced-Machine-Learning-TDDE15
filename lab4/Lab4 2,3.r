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

plotGP = function(fmean, fvar, Xstar, X, y) {
  subheading = sprintf("sigmaN = %g  %i samples  %i observed",
                       sigmaNoise,
                       length(Xstar),
                       length(X))
  
  if (!is.null(fvar)) {
    upperConf = (fmean + diag(fvar) * 1.96)[, 1]
    lowerConf = (fmean - diag(fvar) * 1.96)[, 1]
    plot(
      x = Xstar,
      y = fmean,
      col = "blue",
      type = "l",
      ylim = c(min(lowerConf), max(upperConf)),
      lwd = 2,
      sub = subheading,
      main = "GP",
      ylab = "temp",
      xlab = "time"
    ) +
      lines(x = X, y = y, col = "black")
    
    polygon(
      x = c(rev(Xstar), Xstar),
      y = c(rev(upperConf), lowerConf),
      col = rgb(0, 0, 0, 0.3)
    )
  } else {
    plot(
      x = Xstar,
      y = fmean,
      col = "blue",
      type = "l",
      lwd = 2,
      sub = subheading,
      main = "GP",
      ylab = "temp",
      xlab = "time"
    ) +
      lines(x = X, y = y, col = "black")
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


## Assignment 1

X = c(1, 3, 4)
Xstar = c(2, 3, 4)
kernelMatrix(kernel = kernelWrap(sigmaF = 20, l = 0.2),
             x = X,
             y = Xstar)

## Assignment 2

polyFit <-
  lm(small.data$temp ~  small.data$time + I(small.data$time ^ 2))
sigmaNoise = sd(polyFit$residuals)

GPfit <- gausspr(
  small.data$time,
  small.data$temp,
  kpar = list(sigmaF = 20, l = 0.2),
  kernel = kernelWrap,
  var = sigmaNoise ^ 2
)

fstar = predict(GPfit, small.data$time)
plotGP(fstar, NULL, small.data$time, small.data$time, small.data$temp)

## Assignment 3
fvar = getGPvar(
  small.data$time,
  small.data$temp,
  small.data$time,
  sigmaNoise,
  kernelWrap(sigmaF = 20, l = 0.2)
)
plotGP(fstar, fvar, small.data$time, small.data$time, small.data$temp)
  
## Assignment 4


day.polyFit <-
  lm(small.data$temp ~  small.data$day + I(small.data$day ^ 2))
day.sigmaNoise = sd(day.polyFit$residuals)

day.GPfit <- gausspr(
  small.data$day,
  small.data$temp,
  kpar = list(sigmaF = 20, l = 0.2),
  kernel = kernelWrap,
  var = day.sigmaNoise ^ 2
)

day.fstar = predict(day.GPfit, small.data$day)
lines(x = small.data$time,
      day.fstar,
      col = "red",
      lwd = 2)

## Assignment 5

period.GPfit <- gausspr(
  small.data$time,
  small.data$temp,
  kpar = list(
    sigmaF = 20,
    l1 = 1,
    l2 = 10,
    d = 365 / sd(small.data$time)
  ),
  kernel = kernelWrapPeriod,
  var = sigmaNoise ^ 2
)

period.fstar = predict(period.GPfit, small.data$time)
period.fvar = getGPvar(
  small.data$time,
  small.data$temp,
  small.data$time,
  sigmaNoise,
  kernelWrapPeriod(
    sigmaF = 20,
    l1 = 1,
    l2 = 10,
    d = 365 / sd(small.data$time)
  )
)

plotGP(period.fstar,
       period.fvar,
       small.data$time,
       small.data$time,
       small.data$temp)
