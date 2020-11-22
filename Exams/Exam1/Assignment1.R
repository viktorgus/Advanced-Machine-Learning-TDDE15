library(bnlearn)

# Sample graphs
ss = 1000
graphs = random.graph(c("A","B","C","D","E"),num=ss,method="melancon",every=50,burn.in=30000)
graphs<-unique(graphs)

# A DAG has a UAG with the same independencies if no unshielded colliders are present i the DAG
# The moralization of a DAG does not any edges if no unshielded colliders are present
hasMarkov<-sapply(graphs,function(graph){all.equal(moral(graph),skeleton(graph))==TRUE})

num_graphs_has_markov = length(which(hasMarkov==TRUE))

print(sprintf("Fraction of %i graphs that have a markov network with the same independencies: %g",
              ss, num_graphs_has_markov/ss))

# "Fraction of 10000 graphs that have a markov network with the same independencies: 0.2516"