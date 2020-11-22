library(bnlearn)
library(gRain)

# C: Car door
# D: your door choice
# M: monthy hall door choice

structure = model2network("[C][D][M|D:C]")
plot(structure)


cStates = c("CAR_BEHIND_DOOR1", "CAR_BEHIND_DOOR2", "CAR_BEHIND_DOOR3")
dStates = c("CHOOSE_DOOR1", "CHOOSE_DOOR2", "CHOOSE_DOOR3")
mStates = c("OPEN_DOOR1", "OPEN_DOOR2", "OPEN_DOOR3")
par.C = matrix(c(1, 1, 1)/3, ncol = 3, dimnames = list(NULL, cStates))
par.D = matrix(c(1, 1, 1)/3, ncol = 3, dimnames = list(NULL, dStates))

par.M = c(0, 0.5, 0.5, 0, 0, 1, 0, 1, 0, 
          0, 0, 1, 0.5, 0, 0.5, 1, 0, 0, 
          0, 1, 0, 1, 0, 0, 0.5, 0.5, 0)

dim(par.M) = c(3, 3, 3)
dimnames(par.M) = list("M" = mStates, "C" =  cStates, "D" = dStates)

fitted_bn = custom.fit(structure, dist = list(M = par.M, D = par.D, C = par.C))


grainBN = as.grain(fitted_bn)
grainTree <- compile(grainBN)

get_probs = function(grainTree, nodes, states, goalNode){
  
  # Probability for first door
  evidence <- setEvidence(object = grainTree,
                          nodes = nodes,
                          states = states)
  
  probs = querygrain(object = evidence,
                     nodes = goalNode)[goalNode]
  return(probs)
}


print("probs C: ")
print(par.C)

probCar1 = get_probs(grainTree,c("D", "M"), c("CHOOSE_DOOR1", "OPEN_DOOR2"), "C")
probCar2 = get_probs(grainTree,c("D", "M"), c("CHOOSE_DOOR1", "OPEN_DOOR3"), "C")

