
library(ggplot2)

arrows <- c("^", ">", "v", "<")
action_deltas <- list(c(1,0), # up
                      c(0,1), # right
                      c(-1,0), # down
                      c(0,-1)) # left

MovingAverage <- function(x, n){
  
  cx <- c(0,cumsum(x))
  rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  
  return (rsum)
}


GreedyPolicy <- function(x, y){

  action = which(q_table[x,y,]==max(q_table[x,y,]))
  if(length(action)>1){
    action = sample(action, 1)
  }
  return(action)
}

alt_GreedyPolicy <- function(x, y){
  
  action = which(alt_q_table[x,y,]==max(alt_q_table[x,y,]))
  if(length(action)>1){
    action = sample(action, 1)
  }
  return(action)
}

alt_EpsilonGreedyPolicy <- function(x, y, epsilon){
  
  # Get an epsilon-greedy action for state (x,y) from q_table.
  #
  # Args:
  #   x, y: state coordinates.
  #   epsilon: probability of acting randomly.
  # 
  # Returns:
  #   An action, i.e. integer in {1,2,3,4}.
  
  # Your code here.
  if (runif(1)>epsilon){
    return(alt_GreedyPolicy(x,y))
  }else{
    return(sample(1:4,1))
  }
}

EpsilonGreedyPolicy <- function(x, y, epsilon){
  
  # Get an epsilon-greedy action for state (x,y) from q_table.
  #
  # Args:
  #   x, y: state coordinates.
  #   epsilon: probability of acting randomly.
  # 
  # Returns:
  #   An action, i.e. integer in {1,2,3,4}.
  
  # Your code here.
  if (runif(1)>epsilon){
    return(GreedyPolicy(x,y))
  }else{
    return(sample(1:4,1))
  }
}

transition_model <- function(x, y, action, beta){

  
  delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
  final_action <- ((action + delta + 3) %% 4) + 1
  foo <- c(x,y) + unlist(action_deltas[final_action])
  foo <- pmax(c(1,1),pmin(foo,c(H,W)))
  
  return (foo)
}

q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                       beta = 0){
  
  summed_rewards = 0
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
    summed_rewards = summed_rewards+ reward
    if(reward>=0)
      # End episode.
      return (c(summed_rewards,episode_correction))
    
  }
  
}

alt_q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                       beta = 0){
  
  current_state = start_state
  summed_rewards = 0
  repeat{
    # Follow policy, execute action, get reward.
    
    # Q-table update.
    episode_correction = 0
    
    ## Get action and reward for current state
    action = alt_EpsilonGreedyPolicy(current_state[1],current_state[2], epsilon)
    new_state = transition_model(current_state[1],current_state[2], action=action, beta=beta)
    reward = reward_map[new_state[1],new_state[2]]
    
    # get a' from new_states current action and s' from new state to get
    # g(s', a')
    next_state_action = alt_EpsilonGreedyPolicy(new_state[1],new_state[2], epsilon)
    qnew = alt_q_table[new_state[1], new_state[2], next_state_action]
    
    # Calculate correction  
    # correction is altered to the alternative algorithm where the max action value
    # of the next state in the q table is replaced by the epsilonGreedy policy action
    # of the new state (qnew)
    correction = alpha*(gamma*qnew +
                          reward-alt_q_table[current_state[1],current_state[2],action])
    
    # Update current Q-table
    alt_q_table[current_state[1],current_state[2],action] <<- 
      alt_q_table[current_state[1],current_state[2],action] + correction
    
    #Sum all corrections
    episode_correction = episode_correction + correction
    current_state = new_state
    
    summed_rewards = summed_rewards + reward
    if(reward>=0)
      # End episode.
      return (c(summed_rewards,episode_correction))
    
  }
}
  
alt_q_test <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                           beta = 0){
  
  summed_rewards = 0
  current_state = start_state
  repeat{
    # Follow policy, execute action, get reward.
    
    
    ## Get action and reward for current state
    action = alt_GreedyPolicy(current_state[1],current_state[2])
    new_state = transition_model(current_state[1],current_state[2], action=action, beta=beta)
    reward = reward_map[new_state[1],new_state[2]]
    
    current_state = new_state
    
    summed_rewards = summed_rewards +  reward
    if(reward>=0)
      # End episode.
      return (summed_rewards)
    
  }
}

q_test <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                       beta = 0){
  
  summed_rewards = 0
  current_state = start_state
  repeat{
    # Follow policy, execute action, get reward.
    
    
    ## Get action and reward for current state
    action = GreedyPolicy(current_state[1],current_state[2])
    new_state = transition_model(current_state[1],current_state[2], action=action, beta=beta)
    reward = reward_map[new_state[1],new_state[2]]
    
    current_state = new_state
    
    summed_rewards = summed_rewards +  reward
    if(reward>=0)
      # End episode.
      return (summed_rewards)
    
  }
}


# Environment C (the effect of beta).

H <- 3
W <- 6

epsilon = 0.5
gamma = 1
beta = 0
alpha = 0.1

# -1 for all positions, except 1,2:5 = -10 and 1,6 = 10
reward_map <- matrix(-1, nrow = H, ncol = W)
reward_map[1,2:5] <- -10
reward_map[1,6] <- 10

q_table <- array(0,dim = c(H,W,4))
alt_q_table <- array(0,dim = c(H,W,4))

rewards = c()
alt_rewards = c()


# Train (2)
for(i in 1:5000){
  rew <- q_learning(epsilon=epsilon, alpha=alpha, gamma = gamma,  beta = beta, start_state = c(1,1))
  alt_rew <- alt_q_learning(epsilon=epsilon, alpha=alpha, gamma = gamma,  beta = beta, start_state = c(1,1))
  rewards = c(rewards, rew[1])
  alt_rewards = c(alt_rewards, alt_rew[1])
}
    
graphics.off()
par(mfrow=c(2,1))
plot(MovingAverage(rewards,100),type = "l", xlab="Episode", ylab="Reward", main= "Regular Q-Learning")
plot(MovingAverage(alt_rewards,100),type = "l", xlab="Episode", ylab="Reward", main="Alternative Q-Learning")

print("q_table:")
print(q_table)
print("alt_q_table")
print(alt_q_table)    
mean(rewards)
mean(alt_rewards)
# Test (3)

test.rewards = c()
test.alt_rewards = c()
for(i in 1:5000){
  test.rew <- q_test(epsilon=epsilon, alpha=alpha, gamma = gamma,  beta = beta, start_state = c(runif(1,1,H),runif(1,1,W)))
  test.alt_rew <- alt_q_test(epsilon=epsilon, alpha=alpha, gamma = gamma,  beta = beta, start_state =  c(runif(1,1,H),runif(1,1,W)))
  test.rewards = c(test.rewards, test.rew)
  test.alt_rewards = c(test.alt_rewards, test.alt_rew)
}

graphics.off()
par(mfrow=c(2,1))
plot(MovingAverage(test.rewards,100),type = "l", xlab="Episode", ylab="Reward", main= "Regular Q-Learning")
plot(MovingAverage(test.alt_rewards,100),type = "l", xlab="Episode", ylab="Reward", main="Alternative Q-Learning")

