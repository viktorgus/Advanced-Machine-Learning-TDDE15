
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

nSims = 100
sims = simHMM(model, nSims)
obs = sims$observation

alphas = exp(forward(model,obs))
probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/')

probs101 = t(transProbs) %*% probs.filtered[,100]
probs101 = probs101[,1]/sum(probs101[,1])

barplot(probs101,names.arg = ringStates, ylab="probs", ylim=c(0,1),xlab="State 101")
print("### Probabilities ###")
for(i in 1:length(probs101)){
  print(sprintf("State %i: %g",i,probs101[i]))
}
