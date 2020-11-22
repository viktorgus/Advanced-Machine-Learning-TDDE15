
library(gRain)
library(bnlearn)

data = asia
train_split = 0.8
train_index = sample(1:nrow(data),train_split*nrow(data))
train = data[train_index,]
test = data[-train_index,]

# Gets network structure with pc algorithm
# Lern parameters
get_fitted_bn = function(data){
  eq_graph_structure = pc.stable(data)
  graph_structure = cextend(eq_graph_structure)
  fitted_bn = bn.fit(graph_structure, data)
  return(list(bn=fitted_bn, structure=graph_structure))
}

# Check that node is independent on independent_node with given_node
# I.E confirm D-seperation
# check that p(node|given_node, independent_node) = p(node|given_node)
check_independence = function(tree, node, given_node, independent_node){
  grainBN = as.grain(tree)
  grainTree <- compile(grainBN)

  all_states = c("yes", "no")
  for(i in 1:length(all_states)){
    for(j in 1:length(all_states)){
      #p(node|given_node, independent_node)
      evidence_joint <- setEvidence(object = grainTree,
                              nodes = c(given_node,independent_node),
                              states = c(all_states[i],all_states[j]))
      
      
      probs_joint = querygrain(object = evidence_joint,
                         nodes = node)
      #p(node|given_node)
      evidence_indep <- setEvidence(object = grainTree,
                                    nodes = c(given_node),
                                    states = c(all_states[i]))
      
      
      probs_indep = querygrain(object = evidence_indep,
                               nodes = node)
      # p(node|given_node, independent_node) - p(node|given_node)
      # if zero for all evidence sets then independent
      
      probability_diff = get(node, probs_joint)["yes"]-get(node, probs_indep)["yes"]
      if (probability_diff){
        print(sprintf("p(%s|%s,%s) != p(%s|%s)",node, given_node, independent_node, node, given_node))
        return(FALSE)
      }
    }
  }
  print(sprintf("p(%s|%s,%s) = p(%s|%s)",node, given_node, independent_node, node, given_node))
  print(sprintf("%s independent on %s given %s", node, independent_node, given_node))
  return(TRUE)
  }
fitted = get_fitted_bn(data)
check_independence(fitted$bn, "D", "B", "S")
plot(cpdag(fitted$bn))

set.seed(123)
ss<-5000
x<-random.graph(c("A","B","C","D","E"),num=ss,method="melancon",every=50,burn.in=30000)

y<-unique(x)
z<-lapply(y,cpdag)

r=0
for(i in 1:length(y)) {
  if(all.equal(y[[i]],z[[i]])==TRUE)
    r<-r+1
}
r/length(y)

