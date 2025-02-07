using Graphs

# need to start by understanding how graphs work in julia 

# ALGORITHM (very heavy) pseudocode 

# start by randomly choosing a vertex, since all vertices will have entropy 0. Put it in set 1. 

# then we need to propogate, that is, make the obvious choices: 

# what constitutes an "obvious choice"? think about this more

mutable struct G # probably going to need a partitioned set?
    graph::Graph
    cut::int
    # something here for the two partitioned subsets
    
end
function observe(graph)
    # input: a graph and its edge weights
    # output: 
end

function collapse(graph)
    # calculate the cut value when v is placed into either set
    # use probability to determine whether you're going to use the correct or incorrect choice
    #test
    c = 4
end  

function propogate(graph) 
    # make "obvious choices"
    # update the entropy of the vertices
end