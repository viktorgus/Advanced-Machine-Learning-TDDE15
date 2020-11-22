

K1 = function(sigmaF, l) {
  k1 = function(x,y) {
    r = abs(x-y)
    sigmaF^2 * exp(-0.5 * (r / l)^2)
  }
  class(k1) = "kernel"
  return(k1)
}
K2 = function(sigmaF, l, alpha) {
  k2 = function(x,y) {
    r = abs(x-y)
    sigmaF^2 * (1 + r^2 / (2 * alpha * l^2))^-alpha
  }
  class(k2) = "kernel"
  return(k2)
}
K3 = function(sigmaF, l) {
  k3 = function(x, y) {
    r = abs(x-y)
    sigmaF^2 * (1 + sqrt(3) * r / l) * exp(-sqrt(3) * r / l)
  }
  class(k3) = "kernel"
  return(k3)
}

library("kernlab")
r = seq(0,4,by=0.01)
sigmaF = 1
l = 1

cov1 = kernelMatrix(kernel=K1(sigmaF,1),x=r, y=0)
cov3 = kernelMatrix(kernel=K3(sigmaF,1),x=r, y=0)

graphics.off()
par(mfrow=c(2,2))

alphas = c(1/2, 2, 20)
for(alpha in alphas){
cov2 = kernelMatrix(kernel=K2(sigmaF,1, alpha),x=r, y=0)
plot(x=r, y=cov1, type="l", lwd=2, col="blue", main ="Kernels", sub=sprintf("Blue: k1, red: k2 with alpha=%g, green: k3",alpha))+
  lines(x=r,y=cov2, type="l", lwd=2, col="red")+
  lines(x=r,y=cov3, type="l", lwd=2, col="green")}