library(kernlab)

par(mfrow = c(1, 1))

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
  subheading = sprintf("%i samples  %i observed",
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

getMargLike = function (X, y, Xstar, sigmaNoise, k){
  post = posteriorGP(X, y, Xstar, sigmaNoise, k)
  return ((post$logLike))
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



k = kernelWrap(sigmaF = 20, l = 0.2)

## Grid search over sigma param
best_sigma = 0
best_marg = -Inf


sigmas_grid = seq(0.1, 10, by=0.01)
for (sigma_test in sigmas_grid){
  currMarg = getMargLike(small.data$time, small.data$temp, small.data$time, sigma_test, k)

  if(currMarg>best_marg){
    best_sigma = sigma_test
    best_marg = currMarg
  }
}

posterior = posteriorGP(small.data$time, small.data$temp,small.data$time, best_sigma, k)
plotGP(posterior$fmean, posterior$fvar, small.data$time, small.data$time, small.data$temp)
