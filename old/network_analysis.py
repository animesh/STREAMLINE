import networkx as nx
import numpy as np
#import scipy
#import powerlaw

# load network
A = np.loadtxt('C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Python/scale_free_1500nodes_30000edges_7.5outexp_2.0inexp.txt', usecols=range(1500))
graph=nx.from_numpy_matrix(A, create_using=nx.DiGraph)

### Most important

# Average clustering coefficient
nx.average_clustering(graph)

# Average shortest path length
nx.average_shortest_path_length(graph)

# Global efficiency
nx.global_efficiency(nx.DiGraph.to_undirected(graph))


### Addtitional

# Number of triangles
sum(nx.triangles(nx.DiGraph.to_undirected(graph)))/3

# Average degree connectivity
nx.average_degree_connectivity(graph)


# Small-world coefficients
## nx.sigma(nx.DiGraph.to_undirected(graph))
## nx.omega(nx.DiGraph.to_undirected(graph))

# Power Law fit
#### Use R