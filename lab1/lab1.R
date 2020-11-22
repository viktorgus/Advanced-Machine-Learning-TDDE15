library("bnlearn")

data("asia")

set.seed(12345)
network1 = cpdag(hc(asia, score ="aic"))
network2 = cpdag(hc(asia, score ="bde"))

plot(network1, main="AIC Network")
plot(network2, main="BDeu Network")


print(all.equal(network1,network2))

## Question 2
if (!require(gRain)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RBGL")
}

library(gRain)


getFittedBN = function(trainData){
  graphStructure = pc.stable(trainData)
  fittedBN = bn.fit(cextend(graphStructure), trainData)
  return(fittedBN)
}

getTrueBN = function(trainData){
  graphStructure = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
  fittedBN = bn.fit(cextend(graphStructure), trainData)
  return(fittedBN)
}

predictBN = function(fittedBN, testData, missingNode){
  grainBN = as.grain(fittedBN)
  grainTree <- compile(grainBN)
  preds = rep(NA,nrow(testData))
  state = NULL
  for(i in 1:nrow(testData)){
    for (j in 1:ncol(testData)){
      state[j] = ifelse(testData[i,j]=="yes","yes","no")
    }
    evidence <- setEvidence(object = grainTree,
                            nodes = colnames(testData),
                            states = state)
    

    probs = querygrain(object = evidence,
               nodes = missingNode)$S

    preds[i] = ifelse(probs["yes"]>0.5,"yes","no")
  }
  return(preds)
}

index = sample(seq(0,nrow(asia)), round(nrow(asia)*0.8), replace = F)
train = asia[index,]
test = asia[-index,]
testWithMissing = subset(test, select=-S)


fittedBN = getFittedBN(train)
fittedBN.true = getTrueBN(train)

confMatrix = table(test$S,predictBN(fittedBN, testWithMissing, "S"))
confMatrix.true = table(test$S,predictBN(fittedBN.true, testWithMissing, "S"))

print("Network structure learned using P.C")
print(confMatrix)
print(paste("Accuracy: ", sum(diag(confMatrix))/sum(confMatrix)))

print("True network structure")
print(confMatrix.true)
print(paste("Accuracy: ", sum(diag(confMatrix.true))/sum(confMatrix.true)))

