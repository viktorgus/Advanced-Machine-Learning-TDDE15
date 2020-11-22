library(AtmRay)
library(kernlab)

# plots cov1, cov2 against target with contours given a trained GP modell GPfit
makeGridPlot = function(GPfit,
                        targetName,
                        cov1Name,
                        cov2Name,
                        data) {
  cov1Col = which(names(data) == cov1Name)
  cov2Col = which(names(data) == cov2Name)
  targetCol = which(names(data) == targetName)
  
  x1 <- seq(min(data$varWave), max(data$varWave), length = 100)
  x2 <- seq(min(data$skewWave), max(data$skewWave), length = 100)
  gridPoints <- meshgrid(x1, x2)
  gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
  gridPoints <- data.frame(gridPoints)
  names(gridPoints) <- c(names(data)[1], names(data)[2])
  
  probPreds <- predict(GPfit, gridPoints, type = "probabilities")
  contour(
    x1,
    x2,
    matrix(probPreds[, 2], 100, byrow = TRUE),
    20,
    xlab = cov1Name,
    ylab = cov2Name,
    main = sprintf('Prob(%s) - 1:blue  0:red', targetName)
  )
  
  points(data[data[, targetCol] == "0", cov1Col], data[data[, targetCol] == "0", cov2Col], col =
           "red")
  points(data[data[, targetCol] == "1", cov1Col], data[data[, targetCol] == "1", cov2Col], col =
           "blue")
}

#Prints accuracy and confusion tables given predictions and targets
printAccAndConf = function(label, preds, real) {
  conf = table(preds, real)
  acc = sum(diag(conf)) / sum(conf)
  print(sprintf("%s Confusion table", label))
  print(conf)
  print(sprintf("%s Accuracy: %g", label, acc))
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
test = data[-trainIndex, ]

## Assignment 1
GPfit <- gausspr(fraud ~  varWave + skewWave, data = train)
makeGridPlot(GPfit, "fraud", "varWave", "skewWave", train)

# conf matrix
train_predict = predict(GPfit, train)
printAccAndConf("Train", train_predict, train$fraud)

## Assignment 2
test_predict = predict(GPfit, test)
printAccAndConf("Test", test_predict, test$fraud)

# Assignment 3
GPfitAll <- gausspr(fraud ~  ., data = train)
all_test_predict = predict(GPfitAll, test)
printAccAndConf("Test", all_test_predict, test$fraud)
makeGridPlot(GPfitAll, "fraud", "varWave", "skewWave", test)

