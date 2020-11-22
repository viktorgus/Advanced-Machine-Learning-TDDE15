## Libraries
library("HMM")


## Transition probability matrix
## Rows corresponding to starting state, column to next state
trans_probs=matrix(c(0,1,0,0,0,0,0,0,0,0,
                     0,0.5,0.5,0,0,0,0,0,0,0,
                     0,0,0,1,0,0,0,0,0,0,
                     0,0,0,0.5,0.5,0,0,0,0,0,
                     0,0,0,0,0,1,0,0,0,0,
                     0,0,0,0,0,0.5,0.5,0,0,0,
                     0,0,0,0,0,0,0,1,0,0,
                     0,0,0,0,0,0,0,0.5,0.5,0,
                     0,0,0,0,0,0,0,0,0,1,
                     0.5,0,0,0,0,0,0,0,0,0.5),nrow=10,ncol=10, byrow=TRUE)


emissionProbs = matrix(c(0.2,0.2,0.2,0,0,0,0,0,0.2,0.2,
                         0.2,0.2,0.2,0.2,0,0,0,0,0,0.2,
                         0.2,0.2,0.2,0.2,0.2,0,0,0,0,0,
                         0,0.2,0.2,0.2,0.2,0.2,0,0,0,0,
                         0,0,0.2,0.2,0.2,0.2,0.2,0,0,0,
                         0,0,0,0.2,0.2,0.2,0.2,0.2,0,0,
                         0,0,0,0,0.2,0.2,0.2,0.2,0.2,0,
                         0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,
                         0.2,0,0,0,0,0,0.2,0.2,0.2,0.2,
                         0.2,0.2,0,0,0,0,0,0.2,0.2,0.2),nrow=10, ncol=10, byrow=TRUE)



## Initiate HMM Model with hidden states and emission states
## Transition probabilities and emission probabilities
ringStates = seq(1,10,1)
sensorStates = seq(1,10,1)
model = initHMM(States=ringStates, Symbols = sensorStates, 
                transProbs= trans_probs, emissionProbs = emissionProbs)

## Simulate from HMM
nSims = 100
sims = simHMM(model, nSims)
obs = sims$observation
