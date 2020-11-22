
## Libraries
library("HMM")

#states: rain,rain,   rain,sun  sun,rain    sun,sun

## Transition probability matrix
## Rows corresponding to starting state, column to next state


trans_probs = matrix(  c(0.75, 0.25, 0, 0,
                         0, 0, 0.5, 0.5,
                         0.5, 0.5, 0, 0,
                         0, 0, 0.25, 0.75),nrow=4, ncol=4, byrow=TRUE)
colnames(trans_probs) = c("RAIN, RAIN", "RAIN, SUN", "SUN, RAIN", "SUN, SUN")

# col1: rain col2: sun
emission_probs=matrix(c(0.9, 0.1,
                        0.1, 0.9,
                        0.9, 0.1,
                        0.1, 0.9),nrow=4,ncol=2, byrow=TRUE)
colnames(emission_probs) = c("RAIN", "SUN")

## Initiate HMM Model with hidden states and emission states
## Transition probabilities and emission probabilities
model = initHMM(States=colnames(trans_probs), Symbols = colnames(emission_probs), 
                transProbs= trans_probs, emissionProbs = emission_probs)

## Simulate from HMM
nSims = 10
sims = simHMM(model, nSims)
obs = sims$observation
barplot(table(obs))
