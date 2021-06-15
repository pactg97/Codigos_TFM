# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:07:58 2021

@author: 34625
"""


import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from itertools import chain, combinations

def normalized_cut(G,c):
    A=nx.adjacency_matrix(G)
    n=G.order()
    corte=sum([sum([A[i,j] for j in range(n) if(j in c[1])]) for i in range(n) if(i in c[0])])
    assoc1=sum([sum([A[i,j] for j in range(n)]) for i in range(n) if(i in c[0])])
    assoc2=sum([sum([A[i,j] for j in range(n)]) for i in range(n) if(i in c[1])])
    return corte/assoc1 + corte/assoc2

## A partir de nuestro grafo calculamos la matriz L y el autovector asociado al segundo autovalor más pequeño
def normalized_min_cut(graph):
    """Clusters graph nodes according to normalized minimum cut algorithm.
    All nodes must have at least 1 edge. Uses zero as decision boundary. 
    
    Parameters
    -----------
        graph: a networkx graph to cluster
        
    Returns
    -----------
        vector containing -1 or 1 for every node
    References
    ----------
        J. Shi and J. Malik, *Normalized Cuts and Image Segmentation*, 
        IEEE Transactions on Pattern Analysis and Machine Learning, vol. 22, pp. 888-905
    """
    nodos=list(graph.nodes())
    if(nx.is_connected(graph)):
        m_adjacency = np.array(nx.to_numpy_matrix(graph))
    
        D = np.diag(np.sum(m_adjacency, 0))
        D_half_inv = np.diag(1.0 / np.sqrt(np.sum(m_adjacency, 0)))
        M = np.dot(D_half_inv, np.dot((D - m_adjacency), D_half_inv))
    
        (w, v) = np.linalg.eig(M)
        #find index of second smallest eigenvalue
        index = np.argsort(w)[1]
    
        v_partition = v[:, index]
        v_partition = np.sign(v_partition)
        A=set()
        B=set()
        colors=[]
        for i in range(len(v_partition)):
            if(v_partition[i]>0): 
                A.add(nodos[i])
                colors.append("red")
            else: 
                B.add(nodos[i])
                colors.append("green")
        return [B,A],colors
    else:
        comunidades=list(nx.connected_components(graph))
        colores=["red" for i in range(graph.order())]
        return comunidades,colores

def comunidades_por_etapa(graph):
    n=graph.order()
    comunidades_inicial,_=normalized_min_cut(graph)
    conjunto_particiones=[tuple(comunidades_inicial)]
    while(len(conjunto_particiones[-1])<= n-1):
        ultima_particion=conjunto_particiones[-1]
        nueva_particion=[]
        for comunidad in list(ultima_particion):
            if(len(comunidad)==1):
                nueva_particion.append(comunidad)
            else:
                subgrafo=graph.subgraph(list(comunidad))
                biparticion,_=normalized_min_cut(subgrafo)
                nueva_particion.append(biparticion[0])
                nueva_particion.append(biparticion[1])
        conjunto_particiones.append(tuple(nueva_particion))
    return conjunto_particiones

## Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"

## Represento gráficamente el grafo original
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
G7=nx.relabel_nodes(G6,mapping)
plt.figure(figsize=(8, 6))
biparticion,colores=normalized_min_cut(G7)
nx.draw_shell(G7,with_labels=True,node_color=colores,node_size=700,font_size=18)
#plt.savefig("zachary_norm_min_cut.jpg")

communities=comunidades_por_etapa(G7)
# building initial dict of node_id to each possible subset:
node_id = 0
init_node2community_dict = {node_id: communities[0][0].union(communities[0][1])}
for comm in communities:
    for subset in list(comm):
        if subset not in init_node2community_dict.values():
            node_id += 1
            init_node2community_dict[node_id] = subset

# turning this dictionary to the desired format in @mdml's answer
node_id_to_children = {e: [] for e in init_node2community_dict.keys()}
for node_id1, node_id2 in combinations(init_node2community_dict.keys(), 2):
    for node_id_parent, group in init_node2community_dict.items():
        if len(init_node2community_dict[node_id1].intersection(init_node2community_dict[node_id2])) == 0 and group == init_node2community_dict[node_id1].union(init_node2community_dict[node_id2]):
            node_id_to_children[node_id_parent].append(node_id1)
            node_id_to_children[node_id_parent].append(node_id2)

# also recording node_labels dict for the correct label for dendrogram leaves
node_labels = dict()
for node_id, group in init_node2community_dict.items():
    if len(group) == 1:
        node_labels[node_id] = list(group)[0]
    else:
        node_labels[node_id] = ''

# also needing a subset to rank dict to later know within all k-length merges which came first
subset_rank_dict = dict()
rank = 0
for e in communities[::-1]:
    for p in list(e):
        if tuple(p) not in subset_rank_dict:
            subset_rank_dict[tuple(sorted(p))] = rank
            rank += 1
subset_rank_dict[tuple(sorted(chain.from_iterable(communities[-1])))] = rank

# my function to get a merge height so that it is unique (probably not that efficient)
def get_merge_height(sub):
    sub_tuple = tuple(sorted([node_labels[i] for i in sub]))
    n = len(sub_tuple)
    other_same_len_merges = {k: v for k, v in subset_rank_dict.items() if len(k) == n}
    min_rank, max_rank = min(other_same_len_merges.values()), max(other_same_len_merges.values())
    range = (max_rank-min_rank) if max_rank > min_rank else 1
    return float(len(sub)) + 0.8 * (subset_rank_dict[sub_tuple] - min_rank) / range

# finally using @mdml's magic, slightly modified:
G           = nx.DiGraph(node_id_to_children)
nodes       = G.nodes()
leaves      = set( n for n in nodes if G.out_degree(n) == 0 )
inner_nodes = [ n for n in nodes if G.out_degree(n) > 0 ]

# Compute the size of each subtree
subtree = dict( (n, [n]) for n in leaves )
for u in inner_nodes:
    children = set()
    node_list = list(node_id_to_children[u])
    while len(node_list) > 0:
        v = node_list.pop(0)
        children.add( v )
        node_list += node_id_to_children[v]
    subtree[u] = sorted(children & leaves)

inner_nodes.sort(key=lambda n: len(subtree[n])) # <-- order inner nodes ascending by subtree size, root is last

# Construct the linkage matrix
leaves = sorted(leaves)
index  = dict( (tuple([n]), i) for i, n in enumerate(leaves) )
Z = []
k = len(leaves)
for i, n in enumerate(inner_nodes):
    children = node_id_to_children[n]
    x = children[0]
    for y in children[1:]:
        z = tuple(sorted(subtree[x] + subtree[y]))
        i, j = index[tuple(sorted(subtree[x]))], index[tuple(sorted(subtree[y]))]
        Z.append([i, j, get_merge_height(subtree[n]), len(z)]) # <-- float is required by the dendrogram function
        index[z] = k
        subtree[z] = list(z)
        x = z
        k += 1

# dendrogram
plt.figure(figsize=(14, 7))
dendrogram(Z,leaf_font_size=13, labels=[node_labels[node_id] for node_id in leaves])
#plt.savefig('dendrogram_norm_min_cut.jpg')