library(AtmRay)
library(kernlab)

#Prints accuracy and confusion tables given predictions and targets
printAccAndConf = function(label, preds, real) {
  conf = table(preds, real)
  acc = sum(diag(conf)) / sum(conf)
  print(sprintf("%s Confusion table", label))
  print(conf)
  print(sprintf("%s Accuracy: %g", label, acc))
}

#Prints accuracy and confusion tables given predictions and targets
getAcc = function(preds, real) {
  conf = table(preds, real)
  acc = sum(diag(conf)) / sum(conf)
  return(acc)
}

testAcc = function(param=c(0.1), train, validation){
  GPfit <- gausspr(fraud ~  ., data = train, kernel = "rbfdot", kpar=list(sigma=param[1]))
  val_preds = predict(GPfit, validation)
  return(getAcc(val_preds, validation$fraud))
}

data <-
  read.csv(
    "https://github.com/STIMALiU/AdvMLCourse/raw/master/
GaussianProcess/Code/banknoteFraud.csv",
    header = FALSE,
    sep = ","
  )

names(data) <-
  c("varWave", "skewWave", "kurtWave", "entropyWave", "fraud")
data[, 5] <- as.factor(data[, 5])

set.seed(111)
trainIndex <- sample(1:dim(data)[1], size = 1000,
                     replace = FALSE)
train = data[trainIndex, ]
validation = data[-trainIndex, ]


## Grid search over sigma param
best_par = 0.1
best_acc = 0
pars = seq(0, 10, by=0.5)
for (par in pars){
  currAcc = testAcc(par, train, validation)
  if(currAcc>best_acc){
    best_par = par
    best_acc = currAcc
  }
}
print(sprintf("Best sigma: %g  Accuracy: %g", best_par, best_acc))

GPfitAll <- gausspr(fraud ~  ., data = train, kernel="rbfdot", kpar=list(sigma=best_par))
all_test_predict = predict(GPfitAll, validation)
printAccAndConf("Test", all_test_predict, validation$fraud)

#optimal = optim(par = c(0.1),
#                fn = testAcc, train = train, validation = validation, 
#                method="L-BFGS-B",
#                lower = c(.Machine$double.eps),
#                control=list(fnscale=-1))



