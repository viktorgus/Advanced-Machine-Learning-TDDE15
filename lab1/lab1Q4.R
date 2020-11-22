library("bnlearn")

data("asia")

set.seed(12345)

if (!require(gRain)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RBGL")
}


library(gRain)

getFittedBN = function(trainData){
  graphStructure = model2network("[S][A|S][T|S][L|S][B|S][D|S][E|S][X|S]")
  plot(graphStructure, main="Bayes Classifier Network")
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

confMatrix = table(test$S,predictBN(fittedBN, testWithMissing, "S"))

print("Bayes Classifier")
print(confMatrix)
print(paste("Accuracy: ", sum(diag(confMatrix))/sum(confMatrix)))