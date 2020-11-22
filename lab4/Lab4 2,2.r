
## Assignment 1

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

plotGP = function(fmean, fvar, Xstar, X, y){
  upperConf = (fmean+diag(fvar)*1.96)[,1]
  lowerConf = (fmean-diag(fvar)*1.96)[,1]
  subheading = sprintf("sigmaN = %g  %i samples  %i observed", 
                       sigmaNoise, length(Xstar),length(X))
  plot(x=Xstar, y=fmean, col="blue", type="l", lwd=2,
       ylim=c(min(lowerConf),max(upperConf)), 
       sub=subheading, main="GP", ylab="y", xlab="x") +
      points(x=X, y=y, col="green")
  
  polygon(x = c(rev(Xstar), Xstar),
          y = c(rev(upperConf), lowerConf),
          col = rgb(0, 0, 0, 0.3))
  
}

par(mfrow=c(1,1))
#par(mfrow = c(2,2))

## Assignment 2
X = c(0.4)
y = c(0.719)
sigmaNoise = 0.1
Xstar = seq(0,1,length=10)

gpParam = posteriorGP(X, y, Xstar, sigmaNoise, kernelMaker())
plotGP(gpParam$fmean, gpParam$fvar, Xstar, X, y)


## Assignment 3
X = c(0.4, -0.6)
y = c(0.719, -0.044)
sigmaNoise = 0.1
Xstar = seq(-1,1,length=10)

gpParam = posteriorGP(X, y, Xstar, sigmaNoise, kernelMaker())
plotGP(gpParam$fmean, gpParam$fvar, Xstar, X, y)

## Assignment 4

X = c(-1.0, -0.6, -0.2, 0.4, 0.8)
y = c(0.768, -0.044, -0.940, 0.719, -0.664)
sigmaNoise = 0.1
Xstar = seq(-1,1,length=10)
gpParam = posteriorGP(X, y, Xstar, sigmaNoise, kernelMaker())
plotGP(gpParam$fmean, gpParam$fvar, Xstar, X, y)

## Assignment 5
gpParam = posteriorGP(X, y, Xstar, sigmaNoise, kernelMaker(sigmaF=1,l=1))
plotGP(gpParam$fmean, gpParam$fvar, Xstar, X, y)