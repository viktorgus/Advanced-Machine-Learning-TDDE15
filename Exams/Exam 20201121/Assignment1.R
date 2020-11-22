#### Assignment a)

getRandom = function(){
  structure = model2network("[C][A|C][Y|A:C][D|C]")
  
  
  # C
  cparams = runif(1)
  par.C = matrix(c(cparams, 1-cparams), ncol = 1, nrow = 2, dimnames=list(C=c(0,1)))
  
  # y
  yparams = runif(4)
  par.Y = c(yparams[1], 1-yparams[1],
            yparams[2], 1-yparams[2],
            yparams[3], 1-yparams[3],
            yparams[4], 1-yparams[4])
  dim(par.Y) = c(2,2,2)
  dimnames(par.Y) = list("Y"=c(0,1), "C"=c(0,1), "A" =c(0,1))
  
  # A
  aparams = runif(2)
  par.A = matrix(c(aparams[0], 1-aparams[0],
                   aparams[1], 1-aparams[1]), ncol = 2, nrow = 2, dimnames=list(A=c(0,1), C=c(0,1)))
  
  # D
  dparams = runif(2)
  par.D = matrix(c(dparams[1], 1-dparams[1],
                   dparams[2], 1-dparams[2]), ncol = 2, nrow = 2)
  dimnames(par.D) = list("D"=list(0,1), "C"=list(0,1))
  
  # Y
  yparams = runif(4)
  par.Y = c(yparams[1], 1-yparams[1],
            yparams[2], 1-yparams[2],
            yparams[3], 1-yparams[3],
            yparams[4], 1-yparams[4])
  dim(par.Y) = c(2,2,2)
  dimnames(par.Y) = list("Y"=c(0,1), "C"=c(0,1), "A" =c(0,1))
  fitted_bn = custom.fit(structure, dist = list("C" = par.C, "Y" = par.Y, "A" = par.A, "D"=par.D))
  return(fitted_bn)
}

get_probs = function(grainTree, nodes, states, goalNode){
  
  # Probability for first door
  evidence <- setEvidence(object = grainTree,
                          nodes = nodes,
                          states = states)
  
  probs = querygrain(object = evidence,
                     nodes = goalNode)[goalNode]
  return(probs)
}



check_monotone =function(grainTree, over){
  # p(Y=1|A=1, over=1)
  prob1 = get_probs(grainTree, c("A",over), c("1","1"), "Y")$Y[2]
  
  # p(Y=1|A=1, over=0)
  prob2 = get_probs(grainTree, c("A",over), c("1","0"), "Y")$Y[2]
  
  # p(Y=1|A=0, over=0)
  prob3 = get_probs(grainTree, c("A",over), c("1","0"), "Y")$Y[2]
  
  # p(Y=1|A=0, over=1)
  prob4 = get_probs(grainTree, c("A",over), c("1","0"), "Y")$Y[2]
  non_decreasing = prob1>=prob2 && prob2>=prob3
  non_increasing = prob1<=prob2 && prob2<=prob3
  monotone = non_decreasing | non_increasing
  return(monotone)
}

library(bnlearn)
library(gRain)

ss = 1000

# counts of graphs that are monotone in C but not D
count1 = 0

# counts of graphs that are monotone in D but not C
count2 = 0

for(i in 1:ss){
  fitted_bn = getRandom()
  
  grainBN = as.grain(fitted_bn)
  grainTree <- compile(grainBN)
  
  # check i)
  monoC = check_monotone(grainTree, "C")
  monoD = check_monotone(grainTree, "D")
  
  is_1 = monoC && !monoD
  is_2 = monoD && !monoC
  count1 = count1 + as.numeric(is_1)
  count2 = count2 + as.numeric(is_2)
}

print(sprintf("%i of case (i)   %i of case (ii)", count1, count2))