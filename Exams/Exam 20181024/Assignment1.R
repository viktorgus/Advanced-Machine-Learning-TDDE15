library(bnlearn)
library(gRain)

set.seed(567)
data("asia")
ind <- sample(1:5000, 4000)
tr <- asia[ind,]
te <- asia[-ind,]

train_lengths = c(10,20,50,100,1000,2000)

make_naives_structure = function(target_node, other_nodes){
  sub_strings = c(sprintf("[%s]", target_node))
  for(dependant in other_nodes){
    if (!dependant==target_node){
    node_dependency = sprintf("[%s|%s]", dependant,target_node)
    sub_strings = c(sub_strings, node_dependency)}
    
  }

  structure_string = paste(sub_strings, collapse='')
  print(structure_string)
  return(model2network(structure_string))
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

naive_accs = c()
inv_naive_accs = c()
naive_structure = make_naives_structure("S", colnames(tr))
inv_naive_structure = model2network("[T][X][A][B][D][E][L][S|T:X:A:B:D:E:L]")

plot(naive_structure, main="Bayes Classifier")
plot(inv_naive_structure, main="Bayes Classifier")

for(trainLen in train_lengths){
  fittedBN = bn.fit(naive_structure, tr[1:trainLen,], method="bayes")
  predictions = predictBN(fittedBN, subset(te, select=-S), c("S"))
  conf_matrix = table(te$S, predictions)
  accuracy = sum(diag(conf_matrix))/sum(conf_matrix)  
  naive_accs = c(naive_accs, accuracy)
  
  invfittedBN = bn.fit(inv_naive_structure, tr[1:trainLen,], method="bayes")
  invpredictions = predictBN(invfittedBN, subset(te, select=-S), c("S"))
  inv_conf_matrix = table(te$S, invpredictions)
  inv_accuracy = sum(diag(inv_conf_matrix))/sum(inv_conf_matrix)  
  inv_naive_accs = c(inv_naive_accs, inv_accuracy)
}

plot(x=train_lengths, y=naive_accs, main="Accuracies", type="p", col="Green", ylim=c(0,1), lwd=2, sub="Green: Naive Bayes, Red: Inverse Naive Bayes ")
points(x=train_lengths, y=inv_naive_accs, col="Red", lwd=2)
