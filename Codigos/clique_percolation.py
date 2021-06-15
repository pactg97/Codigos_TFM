import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations

## Obtengo Las comunidades a traves de k-cliques vecinos, es decir, que comparte k-1 nodos.
## Los nodos que no están contenidos en ningún k-clique aparecen como comunidades individuales aisladas.
def get_percolated_cliques(G, k):
    cliques = list(frozenset(c) for c in nx.find_cliques(G) if len(c) >= k)
    perc_graph = nx.Graph()
    perc_graph.add_nodes_from(cliques)
    for c1, c2 in combinations(cliques, 2):
        if len(c1.intersection(c2)) >= (k - 1):
            perc_graph.add_edge(c1, c2)

    for component in nx.connected_components(perc_graph):
        yield(frozenset.union(*component))
        
def clustering_k_clique(G,k):
    nodos_aislados=set(range(G.order()))
    comunidades=[set(i) for i in get_percolated_cliques(G,k)]
    for i in comunidades:
        nodos_aislados=nodos_aislados - i
    for i in nodos_aislados:
        comunidades.append({i})
    return comunidades

## Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G2=nx.Graph()
G2.add_nodes_from(range(34))
G2.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])

plt.figure(1)
nx.draw_shell(G2,with_labels=True)
print("Ejemplo")
print(clustering_k_clique(G2,3))