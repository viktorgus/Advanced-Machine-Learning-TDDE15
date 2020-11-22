
library("HMM")


transProbs=t(matrix(c(0.5,0.5,0,0,0,0,0,0,0,0,
                    0,0.5,0.5,0,0,0,0,0,0,0,
                    0,0,0.5,0.5,0,0,0,0,0,0,
                    0,0,0,0.5,0.5,0,0,0,0,0,
                    0,0,0,0,0.5,0.5,0,0,0,0,
                    0,0,0,0,0,0.5,0.5,0,0,0,
                    0,0,0,0,0,0,0.5,0.5,0,0,
                    0,0,0,0,0,0,0,0.5,0.5,0,
                    0,0,0,0,0,0,0,0,0.5,0.5,
                    0.5,0,0,0,0,0,0,0,0,0.5),10))

emissionProbs = t(matrix(c(0.2,0.2,0.2,0,0,0,0,0,0.2,0.2,
                         0.2,0.2,0.2,0.2,0,0,0,0,0,0.2,
                         0.2,0.2,0.2,0.2,0.2,0,0,0,0,0,
                         0,0.2,0.2,0.2,0.2,0.2,0,0,0,0,
                         0,0,0.2,0.2,0.2,0.2,0.2,0,0,0,
                         0,0,0,0.2,0.2,0.2,0.2,0.2,0,0,
                         0,0,0,0,0.2,0.2,0.2,0.2,0.2,0,
                         0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,
                         0.2,0,0,0,0,0,0.2,0.2,0.2,0.2,
                         0.2,0.2,0,0,0,0,0,0.2,0.2,0.2),10))

ringStates = seq(1,10,1)
sensorStates = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
model = initHMM(States=ringStates, Symbols = sensorStates, transProbs= transProbs, emissionProbs = emissionProbs)

getAccuracies = function(nSims, model){

  
  sims = simHMM(model, nSims)
  obs = sims$observation
  
  alphas = exp(forward(model,obs))
  betas = exp(backward(model,obs))
  
  probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/') 
  probs.smoothed = sweep(alphas*betas, 2, apply(alphas*betas,2,sum), '/') 
  
  predStates.filtered = apply(probs.filtered,2,which.max)
  predStates.smoothed = apply(probs.smoothed,2,which.max)
  maxProbPath = viterbi(model, obs)
  
  conf.maxProbPath = table(maxProbPath, sims$states)
  conf.filtered = table(predStates.filtered, sims$states)
  conf.smoothed = table(predStates.smoothed, sims$states)
  
  
  acc.filtered = sum(diag(conf.filtered))/sum(conf.filtered)
  acc.smoothed = sum(diag(conf.smoothed))/sum(conf.smoothed)
  acc.maxProbPath = sum(diag(conf.maxProbPath))/sum(conf.maxProbPath)
  
  return(t(c(acc.filtered, acc.smoothed, acc.maxProbPath)))
  }

sims = seq(10,400,20)
maxSims = sims[length(sims)]
minSims = sims[1]

graphics.off()
par(mfrow=c(3,3))
accuracies = matrix(0,nrow=3, ncol=(length(sims)))
labels = c("Filtered", "Smoothed", "Max probable path")

for (i in 1:length(sims)){
  accuracies[,i] = getAccuracies(sims[i],model)
  print(sims[i]/maxSims)
  flush.console()                          
}

graphics.off()
plot(0,0, xlab="n Simulations", ylab="Accuracy", main="Accuracies for different measures",xlim=c(sims[1],sims[length(sims)]),ylim=c(0,1))
for(i in 1:3 ){
  lines(x=sims, y=accuracies[i,], type="b", col = rainbow(3)[i])
  text(maxSims-maxSims/4*i,0.9,labels[i],col= rainbow(3)[i])
  
}
