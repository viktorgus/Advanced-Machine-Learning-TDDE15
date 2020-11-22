## GENERAL R ---------------------------------------------------------

# remove column, select
subset(test, select=-S)

# Confusion matrix
table(true_vals, predicted)

# Accuracy and miss classification rate (MCR)
conf_matrix = table(true_vals, predicted)
accuracy = sum(diag(conf_matrix))/sum(conf_matrix)
mcr = 1-accuracy

# Prints
print(paste("accuracy: ", acc, "mcr: ", mcr))
print(sprintf("some integer %i: %g",i,some_float))

# Plots

# type = p:points, l:lines, b:both, h:histogram
plot(0,0, xlab="n Simulations", ylab="Accuracy", 
     main="Accuracies for different measures",
     xlim=c(sims[1],sims[length(sims)]),ylim=c(0,1))

barplot(c(acc.filtered, acc.smoothed, acc.maxProbPath), 
        names.arg=c("Filtered", "Smoothed", "Most Probable Path"), 
        ylab="Accuracy", ylim=c(0,1))
#lwd = linewidth
lines()
points()
text()


## BAYESIAN NETOWORKS / DAGS / GRAPHS ---------------------------------

## Libraries
library(gRain)
library(bnlearn)


## Learn the structure using HC (Hill Climbing) with a dataset
## hill-climbing algorithm is a heuristic that performs local search optimization
## Can use aic, bic, bde (BDeu) and loglik (equal to entropy measure)
network = hc(data, score ="aic")

## Get equivalent graph (equivalence class) construct its moral graph
## complete acyclic partially directed graph
## The CPDAGs have the same skeleton as the original networks,
## but with undirected arcs for every reversible arc in the equivalence class
equivalent_network = cpdag(network)

## Check if graphs are equal
is_equal = all.equal(network1, network2)

## Learn Graph Structure I.E equivalence class of a DAG with PC
eq_graph_structure = pc.stable(data)

## Set known Graph Sctructure
graph_structure = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")


## Find the equivalence class and v-structure of a bayesian network
## Constructs moral graph
graph_structure = cextend(eq_graph_structure)

## Find parameters of graph structure given data
fitted_graph = bn.fit(graph_structure, data)

## Query Graph, get probabilities of nodes given observations
## set observation to state
grainBN = as.grain(fitted_graph)
grainTree <- compile(grainBN)

node_names = colnames(testData)
evidence <- setEvidence(object = grainTree,
                        nodes = node_names,
                        states = state)

probs = querygrain(object = evidence,
                   nodes = missingNode)$S

## Random graph
vars = c("A","B","D","E","L","S","T","X") 

set.seed(1234567890) 
random_start = random.graph(vars)
random.graph(c("A","B","C","D","E"),num=ss,method="melancon",every=50,burn.in=30000)

## Empty graph
vars = c("A","B","D","E","L","S","T","X") 
empty_start = empty.graph(vars) 

# Make naive bayes classifier graph structure from nodes
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
naive_structure = make_naives_structure("S", colnames(data))


## PREDICT from bn
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

Approximate = prop.table(table(cpdist(fitted_BN, "S", (X == "yes" & B == "yes"))))

## HIDDEN MARKOV MODELS ------------------------------------------------

## Libraries
library("HMM")


## Transition probability matrix
## Rows corresponding to starting state, column to next state
trans_probs=matrix(c(0.5,0.5,0,0,0,0,0,0,0,0,
                    0,0.5,0.5,0,0,0,0,0,0,0,
                    0,0,0.5,0.5,0,0,0,0,0,0,
                    0,0,0,0.5,0.5,0,0,0,0,0,
                    0,0,0,0,0.5,0.5,0,0,0,0,
                    0,0,0,0,0,0.5,0.5,0,0,0,
                    0,0,0,0,0,0,0.5,0.5,0,0,
                    0,0,0,0,0,0,0,0.5,0.5,0,
                    0,0,0,0,0,0,0,0,0.5,0.5,
                    0.5,0,0,0,0,0,0,0,0,0.5),nrow=10,ncol=10, byrow=TRUE)

emissionProbs = matrix(c(0.1,0.1,0.1,0,0,0,0,0,0.1,0.1,
                         0.1,0.1,0.1,0.1,0,0,0,0,0,0.1,
                         0.1,0.1,0.1,0.1,0.1,0,0,0,0,0,
                         0,0.1,0.1,0.1,0.1,0.1,0,0,0,0,
                         0,0,0.1,0.1,0.1,0.1,0.1,0,0,0,
                         0,0,0,0.1,0.1,0.1,0.1,0.1,0,0,
                         0,0,0,0,0.1,0.1,0.1,0.1,0.1,0,
                         0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,
                         0.1,0,0,0,0,0,0.1,0.1,0.1,0.1,
                         0.1,0.1,0,0,0,0,0,0.1,0.1,0.1),nrow=10, ncol=10, byrow=TRUE)


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

barplot(table(obs))
# Calculate filtered p(zt|x0:t) = a(zt)/sum(a(zt))
#           and smoothed p(zt|x0:T) = a(zt)*b(zt)/sum(a(zt)*b(zt))
#
# a(zt) = p(x0:t, zt) = p(xt|zt)Sum(a(zt-1)*p(zt|zt-1))   /sum over zt-1
# b(zt) = p(xt+1:T|zt) = sum(b(zt+1)*p(xt+1|zt+1)*p(zt+1|z))
posterior(hmm,obs)

alphas = exp(forward(model,obs))
betas = exp(backward(model,obs))

probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/') 
probs.smoothed =sweep(alphas*betas, 2, apply(alphas*betas,2,sum), '/') 

## Maximum probable path (viterbi)
maxProbPath = viterbi(model, obs)

## Plot filtered probability for each state 
graphics.off()
par(mfrow=c(3,4))
for(i in ringStates ){
  plot(probs.filtered[i,], type="b", col = rainbow(10)[i], xlab="t", ylab=sprintf("p(zt=%i, x1:t)",i), main=ifelse(i==1,"Filtered Distrubutions",""))
  plot(probs.smoothed[i,], type="b", col = rainbow(10)[i], xlab="t", ylab=sprintf("p(zt=%i| x1:T)",i), main=ifelse(i==1,"Smoothed Distrubutions",""))
}

## Plot most probable path
plot(maxProbPath,type="b", xlab="t", ylab="state", main="most probable path", col="blue")

## Find most probable state for filtered and smoothed distrubution
predStates.filtered = apply(probs.filtered,2,which.max)
predStates.smoothed = apply(probs.smoothed,2,which.max)

## Get entropy for probabilities
entropy.empirical(probs.filtered[,i])


## REINFORCMENT LEARNING ----------------------------------------------
GreedyPolicy <- function(x, y){
  action = which(q_table[x,y,]==max(q_table[x,y,]))
  if(length(action)>1){
    action = sample(action, 1)
  }
  return(action)
}

EpsilonGreedyPolicy <- function(x, y, epsilon){
  if (runif(1)>epsilon){
    return(GreedyPolicy(x,y))
  }else{
    return(sample(1:4))
  }
}

q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                       beta = 0){

  current_state = start_state
  repeat{
    # Follow policy, execute action, get reward.
    
    # Q-table update.
    episode_correction = 0
    action = EpsilonGreedyPolicy(current_state[1],current_state[2], epsilon)
    new_state = transition_model(current_state[1],current_state[2], action=action, beta=beta)
    reward = reward_map[new_state[1],new_state[2]]
    
    # Calculate correction  
    correction = alpha*(gamma*max(q_table[new_state[1],new_state[2],]) +
                          reward-q_table[current_state[1],current_state[2],action])
    
    # Update current Q-table
    q_table[current_state[1],current_state[2],action] <<- 
      q_table[current_state[1],current_state[2],action] + correction
    
    #Sum all corrections
    episode_correction = episode_correction + correction
    current_state = new_state
    
    if(reward!=0)
      # End episode.
      return (c(reward,episode_correction))
    
  }
}


## GAUSSIAN PROCESSES ------------------------------------------------
# library
library(kernlab)

# Kernel with wrapper, calcs entire matrix for vectors x1, x2
# Squared exponential kernel
kernelMaker = function(sigmaF=1,l=0.3){
  
  expKernel <- function(x1,x2){
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaF^2*exp(- 0.5*( (x1-x2[i])/l )^2 )
    }
    return(K)
  }
  # if used with
  #class(expKernel) <- 'kernel'
  
  return(expKernel)
}

# Periodic kernel l1: length of periodic part l2: length of distance
kernelWrapPeriod = function(sigmaF, l1, l2, d) {
  expKernel <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA, n1, n2)
    for (i in 1:n2) {
      K[, i] <-
        sigmaF ^ 2 * exp(-2 * (sin(pi * abs(x1 - x2[i]) / d) / l1) ^ 2) *
        exp(-0.5 * (abs(x1 - x2[i]) / l2) ^ 2)
    }
    return(K)
  }
  class(expKernel) <- 'kernel'
  
  return(expKernel)
}

# Posterior of gaussian process
posteriorGP = function (X, y, Xstar, sigmaNoise, k){
  kstar = k(X,Xstar)
  L = chol(k(X,X)+sigmaNoise^2*diag(length(X)))
  L = t(L)
  alpha = solve(t(L), solve(L,y))
  fMean = t(kstar) %*% alpha
  v = solve(L,kstar)
  fVar = k(Xstar, Xstar)-t(v)%*%v
  logMarginalLike = -0.5*(t(y)%*%alpha)-sum(diag(L))-length(y)/2*log(2*pi)
  return (list(fmean=fMean, fvar=fVar, logLike=logMarginalLike))
}

gpParam = posteriorGP(X, y, Xstar, sigmaNoise, kernelMaker(sigmaF=1,l=1))

## Get Fmean, fvar, with scaling, unscale and plot
polyFit <-
  lm(scale(y) ~  scale(x) + scale(I(x ^ 2)))
sigmaN_scaled = sd(polyFit$residuals)
sigmaN = sigmaN_scaled*sd(y)^2 # Eller kör polyfit fast på data utan scale

sigmaN = 0.05
sigmaN_scaled = sigmaN/sd(y)^2

gpParam = posteriorGP(scale(x), scale(y),
                      scale(xStar), sigmaN_scaled, kernelWrapper(sigmaF = 1, l = 1))
fmean_unscaled = mean(y) + gpParam$fmean*sd(y)
fvar_unscaled = gpParam$fvar*sd(y)^2
plotGP(fmean_unscaled, fvar_unscaled, xStar, x, y, sigmaN)




# PLOT GP plots gp given function mean, variance, and x values. Plots variance for y if sigmaN is provided
plotGP = function(fmean,
                  fvar = NULL,
                  Xstar,
                  X,
                  y,
                  sigmaN = 0,
                  xlab = "time",
                  ylab = "temp") {
  subheading = sprintf("%i samples  %i observed",
                       length(Xstar),
                       length(X))
  legendPosition = "topright"
  # IF VARIANCE GIVEN
  if (!is.null(fvar)) {
    upperConf = (fmean + sqrt(diag(fvar)) * 1.96)[, 1]
    lowerConf = (fmean - sqrt(diag(fvar)) * 1.96)[, 1]
    
    #Conf for y
    legendText = c("data", "post mean", "95% intervals for f")
    legendCols = c("green", "blue",  rgb(0, 0, 0, 0.3))
    if (!sigmaN == 0) {
      legendText = c("data",
                     "post mean",
                     "95% intervals for f",
                     "95% intervals for y")
      legendCols = c("green", "blue",  rgb(0, 0, 0, 0.3), rgb(0, 0, 0.8, 0.3))
      Y_upperConf = (fmean + sqrt(diag(fvar) + sigmaN ^ 2) * 1.96)[, 1]
      Y_lowerConf = (fmean - sqrt(diag(fvar) + sigmaN ^ 2) * 1.96)[, 1]
    }
    plot(
      x = Xstar,
      y = fmean,
      col = "blue",
      type = "l",
      ylim = c(min(lowerConf), max(upperConf)),
      lwd = 2,
      sub = subheading,
      main = "GP",
      ylab = ylab,
      xlab = xlab
    ) +
      points(x = X, y = y, col = "green")
    legend(
      legendPosition,
      inset = 0.02,
      legend = legendText,
      col = legendCols,
      pch = c('o', NA, NA, NA),
      lty = c(NA, 1, 1, 1),
      lwd = 2,
      cex = 0.55
    )
    
    polygon(
      x = c(rev(Xstar), Xstar),
      y = c(rev(upperConf), lowerConf),
      col = rgb(0, 0, 0, 0.3)
    )
    if (!sigmaN == 0) {
      polygon(
        x = c(rev(Xstar), Xstar),
        y = c(rev(Y_upperConf), Y_lowerConf),
        col = rgb(0, 0, 0.8, 0.3)
      )
    }
  }
  #IF NO VARIANCE GIVEN
  else {
    plot(
      x = Xstar,
      y = fmean,
      col = "blue",
      type = "l",
      lwd = 2,
      sub = subheading,
      main = "GP",
      ylab = ylab,
      xlab = xlab
    ) +
      points(x = X, y = y, col = "green")
    legend(
      legendPosition,
      inset = 0.02,
      legend = c("data", "post mean"),
      col = c("green", "blue"),
      pch = c('o', NA, NA, NA),
      lty = c(NA, 1, 1, 1),
      lwd = 2,
      cex = 0.55
    )
    
  }
  
}


## Get noise in data
day.polyFit <-
  lm(small.data$temp ~  small.data$day + I(small.data$day ^ 2))
day.sigmaNoise = sd(day.polyFit$residuals)

## kernlab prediction and gp fit
GPfit <- gausspr(
  small.data$time,
  small.data$temp,
  kpar = list(sigmaF = 20, l = 0.2),
  kernel = kernelWrap,
  var = sigmaNoise ^ 2
)

fstar = predict(GPfit, small.data$time)


# Get variance over xstar
getVar = function (k, x, xstar, sigmaN){
  n = length(x)
  I = diag(n)
  
  Kss = kernelMatrix(kernel = k, x = xstar, y = xstar)
  Ksx = kernelMatrix(kernel = k, x = xstar, y = x)
  Kxs = t(Ksx)
  Kxx = kernelMatrix(kernel = k, x = x, y = x)
  V = Kss - Ksx %*% solve( Kxx + sigmaN^2 * I ) %*% Kxs
  
  return(diag(V))
}