library("HMM")


## Transition probability matrix
## Rows corresponding to starting state, column to next state
trans_matrix = matrix(0, 100, 100)
for (i in 1:nrow(trans_matrix)){
  trans_matrix[i,i]=0.1
  if(i==nrow(trans_matrix)){
    trans_matrix[i,1]=0.9
  }else{
    trans_matrix[i,i+1]=0.9
  }
}

# Emission for door, output of 1 equals door output of
emission_matrix = matrix(NA, 100, 2)
emission_matrix[,1] = 0.9
emission_matrix[,2] = 0.1
for(i in list(10,11,12,20,21,22,30,31,32)){
  emission_matrix[i,1]=0.1
  emission_matrix[i,2]=0.9
}


## Initiate HMM Model with hidden states and emission states
## Transition probabilities and emission probabilities
segments = seq(1,100,1)
door_states = c("NO_DOOR","DOOR")
model = initHMM(States=segments, Symbols = door_states, 
                transProbs= trans_matrix, emissionProbs = emission_matrix)



## 2)
no_doors = rep("NO_DOOR", 30)
doors = rep("DOOR",3)
obs = c(no_doors, doors)
alphas = exp(forward(model,obs))
probs.filtered = sweep(alphas, 2, apply(alphas,2,sum), '/') 
plot(x=1:100, y=probs.filtered[,length(obs)], type="l", lwd="2", 
     main="Probability of Segment",ylab="p(z|x1:t)", xlab="segment")
likelySegment = which.max(probs.filtered[,length(obs)])

