
# Assignment 1
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

ringStates = seq(1,10,1)
sensorStates = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
model = initHMM(States=ringStates, Symbols = sensorStates, transProbs= transProbs, emissionProbs = emissionProbs)

## Assignment 2
nSims = 100
sims = simHMM(model, nSims)
obs = sims$observation

# Assignment 3
# Calculate filtered p(zt|x0:t) = a(zt)/sum(a(zt))
#           and smoothed p(zt|x0:T) = a(zt)*b(zt)/sum(a(zt)*b(zt))
#
# a(zt) = p(x0:t, zt) = p(xt|zt)Sum(a(zt-1)*p(zt|zt-1))   /sum over zt-1
# b(zt) = p(xt+1:T|zt) = sum(b(zt+1)*p(xt+1|zt+1)*p(zt+1|z))

alphas = exp(forward(model,obs))
betas = exp(backward(model,obs))

probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/') 
probs.smoothed =sweep(alphas*betas, 2, apply(alphas*betas,2,sum), '/') 
maxProbPath = viterbi(model, obs)


#plots FILTERED
graphics.off()
par(mfrow=c(3,4))
for(i in ringStates ){
  plot(probs.filtered[i,], type="b", col = rainbow(10)[i], xlab="t", ylab=sprintf("p(zt=%i, x1:t)",i), main=ifelse(i==1,"Filtered Distrubutions",""))
}

#plots SMOOTHED
graphics.off()
par(mfrow=c(3,4))
for(i in ringStates ){
  plot(probs.smoothed[i,], type="b", col = rainbow(10)[i], xlab="t", ylab=sprintf("p(zt=%i| x1:T)",i), main=ifelse(i==1,"Smoothed Distrubutions",""))
}

#plots MOST PROBABLE PATH
graphics.off()
par(mfrow=c(1,1))
plot(maxProbPath,type="b", xlab="t", ylab="state", main="most probable path", col="blue")



# Assignment 4 

predStates.filtered = apply(probs.filtered,2,which.max)
predStates.smoothed = apply(probs.smoothed,2,which.max)

conf.maxProbPath = table(maxProbPath, sims$states)
conf.filtered = table(predStates.filtered,sims$states)
conf.smoothed = table(predStates.smoothed, sims$states)


acc.filtered = sum(diag(conf.filtered))/sum(conf.filtered)
acc.smoothed = sum(diag(conf.smoothed))/sum(conf.smoothed)
acc.maxProbPath = sum(diag(conf.maxProbPath))/sum(conf.maxProbPath)

print(paste("Accuracy filtered: ", acc.filtered, "Accuracy smoothed: ", acc.smoothed, "Accuracy max probable path: ", acc.maxProbPath))
barplot(c(acc.filtered, acc.smoothed, acc.maxProbPath), names.arg=c("Filtered", "Smoothed", "Most Probable Path"), ylab="Accuracy", ylim=c(0,1))
