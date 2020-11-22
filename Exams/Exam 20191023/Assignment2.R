library(HMM)

## Transition probability matrix
## Rows corresponding to starting state, column to next state
transProbs=matrix(c(0.5,0.5,0,0,0,0,0,0,0,0,
                      0,0.5,0.5,0,0,0,0,0,0,0,
                      0,0,0.5,0.5,0,0,0,0,0,0,
                      0,0,0,0.5,0.5,0,0,0,0,0,
                      0,0,0,0,0.5,0.5,0,0,0,0,
                      0,0,0,0,0,0.5,0.5,0,0,0,
                      0,0,0,0,0,0,0.5,0.5,0,0,
                      0,0,0,0,0,0,0,0.5,0.5,0,
                      0,0,0,0,0,0,0,0,0.5,0.5,
                      0.5,0,0,0,0,0,0,0,0,0.5),nrow=10,ncol=10, byrow=TRUE)

emissionProbs = matrix(c(0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,0.5,
                           0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,0.5,
                           0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,0.5,
                           0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,0.5,
                           0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,0.5,
                           0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,0.5,
                           0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,0.5,
                           0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.5,
                           0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,0.5,
                           0.1,0.1,0,0,0,0,0,0.1,0.1,0.1,0.5),nrow=10, ncol=11, byrow=TRUE)


## Initiate HMM Model with hidden states and emission states
## Transition probabilities and emission probabilities
ringStates = seq(1,10,1)
sensorStates = seq(1,11,1)
model = initHMM(States=ringStates, Symbols = sensorStates, 
                transProbs= transProbs, emissionProbs = emissionProbs)

obs = c(1,11,11,11)
vit.path = viterbi(model, obs)


alphas = exp(forward(model,obs))
betas = exp(backward(model,obs))

graphics.off()
par(mfrow=c(2,1))
probs.smoothed = sweep(alphas*betas, 2, apply(alphas*betas,2,sum), '/') 
smooth.path = apply(probs.smoothed, 2, which.max)
plot(x=obs, y = vit.path, type="l", col="blue", main="Most probable paths", sub="Blue: Viterbi  Red: Smooth")
plot(x=obs, y= smooth.path, type="l", col="red")
