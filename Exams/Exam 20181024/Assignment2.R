
transProbs=t(matrix(c(0.5,0.5,0,0,0,0,0,0,0,0,
                      0,0.5,0.5,0,0,0,0,0,0,0,
                      0,0,0.5,0.5,0,0,0,0,0,0,
                      0,0,0,0.5,0.5,0,0,0,0,0,
                      0,0,0,0,0.5,0.5,0,0,0,0,
                      0,0,0,0,0,0.5,0.5,0,0,0,
                      0,0,0,0,0,0,0.5,0.5,0,0,
                      0,0,0,0,0,0,0,0.5,0.5,0,
                      0,0,0,0,0,0,0,0,0.5,0.5,
                      0.5,0,0,0,0,0,0,0,0,0.5),nrow=10,ncol=10))

emissionProbs = t(matrix(c(0.2,0.2,0.2,0,0,0,0,0,0.2,0.2,
                           0.2,0.2,0.2,0.2,0,0,0,0,0,0.2,
                           0.2,0.2,0.2,0.2,0.2,0,0,0,0,0,
                           0,0.2,0.2,0.2,0.2,0.2,0,0,0,0,
                           0,0,0.2,0.2,0.2,0.2,0.2,0,0,0,
                           0,0,0,0.2,0.2,0.2,0.2,0.2,0,0,
                           0,0,0,0,0.2,0.2,0.2,0.2,0.2,0,
                           0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,
                           0.2,0,0,0,0,0,0.2,0.2,0.2,0.2,
                           0.2,0.2,0,0,0,0,0,0.2,0.2,0.2),nrow=10, ncol=10))

# Get observations
ringStates = seq(1,10,1)
sensorStates = seq(1,10,1)
model = initHMM(States=ringStates, Symbols = sensorStates, transProbs= transProbs, emissionProbs = emissionProbs)
nSims = 20

sims = simHMM(model, nSims)
obs = sims$observation
real_alphas = exp(forward(model,obs))


#a(zt) = p(x0:t, zt) = p(xt|zt)Sum(a(zt-1)*p(zt|zt-1))
# a(z0) = p(x0, z0) = p(x0|z0)*p(z0)
# x: emissions z:hidden state

# Uniform distrubution over first state
# make matrix with rows: states, columns: time
alphas = matrix(0, nrow(transProbs), length(obs))

initial_state_prob = 1/nrow(transProbs)
alpha_zero = initial_state_prob*emissionProbs[,obs[1]]
alphas[,1] = alpha_zero
for (i in 2:length(obs)){
  # p(xt|zt)
  p1 = emissionProbs[,obs[i]]

  summedPrevTransitions = 0
  for (j in 1:nrow(alphas)){
    summedPrevTransitions = summedPrevTransitions + alphas[j,i-1]*transProbs[j,]
  }
  alphas[,i] = p1*summedPrevTransitions
}

# Check if equal with filtered distrubutions p(zt|x0:t)
probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/') 
real_probs.filtered = sweep(real_alphas, 2, apply(real_alphas,2,sum), '/') 

probs.filtered
real_probs.filtered
probs.filtered- real_probs.filtered
