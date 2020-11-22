
library("HMM")
library("circlize")


transProbs=matrix(c(0.5,0.5,0,0,0,0,0,0,0,0,
                    0,0.5,0.5,0,0,0,0,0,0,0,
                    0,0,0.5,0.5,0,0,0,0,0,0,
                    0,0,0,0.5,0.5,0,0,0,0,0,
                    0,0,0,0,0.5,0.5,0,0,0,0,
                    0,0,0,0,0,0.5,0.5,0,0,0,
                    0,0,0,0,0,0,0.5,0.5,0,0,
                    0,0,0,0,0,0,0,0.5,0.5,0,
                    0,0,0,0,0,0,0,0,0.5,0.5,
                    0.5,0,0,0,0,0,0,0,0,0.5),10)

emissionProbs = matrix(c(0.2,0.2,0.2,0,0,0,0,0,0.2,0.2,
                         0.2,0.2,0.2,0.2,0,0,0,0,0,0.2,
                         0.2,0.2,0.2,0.2,0.2,0,0,0,0,0,
                         0,0.2,0.2,0.2,0.2,0.2,0,0,0,0,
                         0,0,0.2,0.2,0.2,0.2,0.2,0,0,0,
                         0,0,0,0.2,0.2,0.2,0.2,0.2,0,0,
                         0,0,0,0,0.2,0.2,0.2,0.2,0.2,0,
                         0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,
                         0.2,0,0,0,0,0,0.2,0.2,0.2,0.2,
                         0.2,0.2,0,0,0,0,0,0.2,0.2,0.2),10)

ringStates = seq(1,10,1)
sensorStates = seq(1,10,1)
model = initHMM(States=ringStates, Symbols = sensorStates, transProbs= transProbs, emissionProbs = emissionProbs)
nSims = 20

sims = simHMM(model, nSims)
obs = sims$observation

alphas = exp(forward(model,obs))
betas = exp(backward(model,obs))

probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/') 
predStates.filtered = apply(probs.filtered,2,which.max)



entropies = rep(0,nSims)
for (i in 1:nSims){
  entropies[i] = entropy.empirical(probs.filtered[,i])
}
plot(x=1:nSims, y=entropies, type="l", col="red", xlab="t", ylab="entropy", main="Entropy")
