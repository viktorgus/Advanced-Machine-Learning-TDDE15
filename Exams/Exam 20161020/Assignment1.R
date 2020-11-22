library("bnlearn")
library(gRain)
if (!require(gRain)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RBGL")
}


data("asia")

getFittedBN = function(trainData){
  graphStructure = pc.stable(trainData)
  fittedBN = bn.fit(cextend(graphStructure), trainData)
  return(fittedBN)
}

### Exact p(A|X=TRUE, B=TRUE)
set.seed(12345)
fitted_BN = getFittedBN(asia)
grainBN = as.grain(fitted_BN)
grainTree = compile(grainBN)


evidence <- setEvidence(object = grainTree,
                        nodes = c("X", "B"),
                        states = c("yes", "yes"))


exact_probs = querygrain(object = evidence,
                   nodes = "A")$A
print("Exact:")
print(exact_probs)
## Approx
approx = cpquery(fitted_BN, (A=="yes"), ((X=="yes") & (B=="yes")))
print(sprintf("Approximate: no: %g  yes: %g", 1-approx, approx))

## Approx 
Approximate = prop.table(table(cpdist(fitted_BN, "S", (X == "yes" & B == "yes"))))
## Approx 2
sims = rbn(fitted_BN, 4000)
approx2 = which(sims$A=="yes" && sims$X=="yes" && sims$B=="yes")/nrow(sims)
