import networkx as nx
import numpy as np
import scipy
import powerlaw

# load network
A = np.loadtxt('C:/Users/Niclas Popp/Documents/Niclas/Studium/8. SS 21/Helmholtz/Python/scale_free_1500nodes_30000edges_7.5outexp_2.0inexp.txt', usecols=range(1500))
graph=nx.from_numpy_matrix(A, create_using=nx.DiGraph)

# Average clustering coefficient
nx.average_clustering(graph)

# Average shortest path length
nx.average_shortest_path_length(graph)

# Number of triangles
sum(nx.triangles(nx.DiGraph.to_undirected(graph)))/3


# Small-world coefficients
## nx.sigma(nx.DiGraph.to_undirected(graph))
## nx.omega(nx.DiGraph.to_undirected(graph))

# Power Law fit
degrees = {}
for node in graph.nodes_iter():
    key = len(graph.neighbors(node))
    degrees[key] = degrees.get(key, 0) + 1

max_degree = max(degrees.keys(), key=int)
num_nodes = []
for i in range(1, max_degree + 1):
    num_nodes.append(degrees.get(i, 0))

fit = powerlaw.Fit(nodes)
print(fit.power_law.alpha)