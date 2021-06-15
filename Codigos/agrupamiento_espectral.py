# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:42:11 2020

@author: 34625
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from itertools import chain, combinations

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## Transformo el problema de partición en uno de bipartición donde mis variables de decisión serán:
## s_i=1 si i se clasifica en la comunidad 1 y
## s_i=-1 si i se clasifica en la comunidad 2
## De esta forma nuestra función objetivo a maximizar sería s_traspuesta*B*s
## donde B es la matriz utilizada para la modularidad, B_{ij}= A_{ij}-k(i)k(j)/2m.
## Este algoritmo aproxima la solución con el autovector asociado al autovalor máximo.
## Continuo haciendo biparticiones de las comunidades obtenidas con el mismo método mientras mejore 
## la modularidad.

def vector_colores(c,n):
    colores=["lightgreen","green","lightcoral","red","black","white","purple","orange","brown"]
    col=["red" for i in range(n)]
    n_comp=len(c)
    for i in range(n_comp):
        comp=list(c[i])
        for j in comp:
            col[j]=colores[i]
    return col

## Calcula el autovector asociado al autovalor máximo de B
def autovector_max(B,n):
    (autovalores,autovectores)=np.linalg.eig(B)
    maxi=max(autovalores)
    posicion_max=np.where(autovalores==maxi)
    return list(autovectores[:,posicion_max])

## Calcula la matriz B de la modularidad
def calculo_B(A,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-k[i]*k[j]/(2*m)
    return B

## Calcula la modularidad para una matriz B y una partición c
def modularity(B,c,m,n):
    mod=0
    for comp in c:
        for i in comp:
            for j in comp:
                mod=mod+B[i,j]
    mod=mod/(2*m)
    return mod

def mod_local(B,c):
    mod=0
    for comp in c:
        for i in comp:
            for j in comp:
                mod=mod+B[i,j]
    return mod

def modularity2(G,c):
    n=G.order()
    k=[G.degree(i) for i in range(n)]
    aristas=list(G.edges())
    n_aristas=G.number_of_edges()
    modularidad=0
    for comunidad in c:
        n_aristas_dentro=len([(i,j) for (i,j) in aristas if(i in comunidad and j in comunidad)])
        ds=sum([k[i] for i in range(n) if(i in comunidad)])
        modularidad=modularidad+n_aristas_dentro/n_aristas-((ds/(2*n_aristas))*(ds/(2*n_aristas)))
    return modularidad

## Realiza una bipartición sobre cada componente de una partición c de forma que 
## si no mejora la modularidad no realiza la bipartición
def clustering_cada_componente_en_dos(B,m,n,c):
    nueva_c=[]
    for comp in c:
        if(len(comp)!=1):
            ci=list(comp)
            n1=len(ci)
            B1=np.zeros([n1,n1])
            for i in range(n1):
                for j in range(n1):
                    if(i==j): B1[i,j]=B[ci[i],ci[j]] - sum([B[ci[i],k] for k in ci])
                    else: B1[i,j]=B[ci[i],ci[j]]
            s=autovector_max(B1,n1)
            A1=set()
            A2=set()
            for i in range(n1):
                if(s[i]>0): A1.add(ci[i])
                else: A2.add(ci[i])
            if(A1==set()):
                A12={i for i in A2}
                nueva_c=nueva_c+[A12]
            else:
                if(A2==set()):
                    A12={i for i in A1}
                    nueva_c=nueva_c+[A12]
                else:
                    nueva_c=nueva_c+[A1,A2]
        else:
            nueva_c=nueva_c+[comp]
    return nueva_c

## Aplicada el algoritmo haciendo biparticiones sobre el conjunto de nodos hasta que la modularidad
## no aumente
def clustering(h):
    g=preproceso(h)
    A=nx.adjacency_matrix(g)
    n=g.order()
    m=g.number_of_edges()
    c=[set(range(n))]
    k=[g.degree(i) for i in range(n)]
    B=calculo_B(A,n,m,k)
    particiones=[]
    while(len(c)<=4):
        c=clustering_cada_componente_en_dos(B,m,n,c)
        particiones.append(c)
    c_final=[{list(h.nodes())[i] for i in list(com)} for com in c]
    return c_final,modularity(B,c,m,n),particiones
## Ejemplo 1: Cuatro K5 unidos por una arista dos a dos
G=nx.Graph()
for i in range(20):
    G.add_node(i)
for i in range(5):
    for j in range(5):
        if(i!=j):
            G.add_edge(i,j)
for i in range(5,10):
    for j in range(5,10):
        if(i!=j):
            G.add_edge(i,j)
for i in range(10,15):
    for j in range(10,15):
        if(i!=j):
            G.add_edge(i,j)
for i in range(15,20):
    for j in range(15,20):
        if(i!=j):
            G.add_edge(i,j)
G.add_edge(0,5)
G.add_edge(5,10)
G.add_edge(10,15)
G.add_edge(15,0)

## Representamos gráficamente el grafo original y el grafo coroleado en función de su
## partición final
plt.figure(1)
nx.draw_shell(G,with_labels=True)
plt.figure(2)
print("Ejemplo1")


## Segundo Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G1=nx.Graph()
G1.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])

G2=nx.Graph()
G2.add_nodes_from(range(34))
G2.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
    
## Representamos gráficamente el grafo original y el grafo coroleado en función de su
## partición final
plt.figure(3)
nx.draw_shell(G1,with_labels=True)
plt.figure(4)
H=preproceso(G1)
nx.draw_shell(preproceso(G1),with_labels=True)
A=nx.adjacency_matrix(G1)
B=nx.adjacency_matrix(G2)
plt.figure(5)
c1=clustering(G2)
plt.figure(6)
c2=clustering(G1)
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
G7=nx.relabel_nodes(G6,mapping)
plt.figure(7)
plt.figure(figsize=(8, 6))
_,_,particiones=clustering(G6)

communities=[tuple(i) for i in particiones]
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