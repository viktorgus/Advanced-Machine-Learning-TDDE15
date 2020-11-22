library(bnlearn)
library(gRain)

# C: Car door
# D: your door choice
# M: monthy hall door choice

structure = model2network("[A][B][C|A:B]")
plot(structure)

par.A = matrix(c(1, 1)/2, ncol = 2, dimnames = list(NULL, c(0,1)))
par.B = matrix(c(1, 1)/2, ncol = 2, dimnames = list(NULL, c(0,1)))

par.C = c(1, 0, 0, 1 ,0, 1, 1, 0)

dim(par.C) = c(2, 2, 2)
dimnames(par.C) = list("C" = c(0,1), "A" = c(0,1), "B" = c(0,1))
par.C

fitted_bn = custom.fit(structure, dist = list(A = par.A, B = par.B, C = par.C))

samples = rbn(fitted_bn, 1000)
## Find model graph
graphics.off()
par(mfrow=c(2,2))
for (i in 1:4){
  random_indexes = sample(1:1000, 50)
  network = hc(samples[random_indexes,], score ="aic")
  plot(network)
}
samples
