import numpy as np
import networkx as nx
from networkx.algorithms import community


adj_matrix = np.array([[0, 1, 1, 1, 1, 0, 0, 0, 0],
                       [1, 0, 1, 0, 1, 0, 0, 0, 0],
                       [1, 1, 0, 1, 1, 0, 0, 0, 0],
                       [1, 0, 1, 0, 1, 1, 1, 0, 0],
                       [1, 1, 1, 1, 0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0, 0, 1, 1, 1],
                       [0, 0, 0, 1, 1, 1, 0, 1, 1],
                       [0, 0, 0, 0, 0, 1, 1, 0, 1],
                       [0, 0, 0, 0, 0, 1, 1, 1, 0]])


G = nx.from_numpy_array(adj_matrix)

groups = community.label_propagation_communities(G) # there also exists a fast label propagation
print(groups)

subgraphs = []

for _, value_set in enumerate(groups): 
    subgraphs.append(G.subgraph(value_set)) # create an induced subgraph for each community 



