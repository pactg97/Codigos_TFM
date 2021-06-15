# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 17:14:48 2020

@author: 34625
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## Este algoritmo trata de encontrar la bipartición del grafo que minimice la función de corte normalizada:
## NCut(A,B)=Cut(A,B)/assoc(A,V) + Cut(B,A)/assoc(B,V)
## El algoritmo aproxima el mínimo de corte normalizado mediante el segundo autovalor más pequeño de la matriz:
## L=D^{-1/2}*(D-A)*D^{-1/2}, donde A es la matriz de adyacencia y D la diagonal de los grados de los nodos
## De esta forma la bipartición vendrá dada por los signos de las coordenadas del autovector asociado a dicho
## autovalor

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
    return [A,B],normalized_cut(graph,[A,B]),colors

## Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G1=nx.Graph()
G1.add_nodes_from(range(34))
G1.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
    
G2=nx.Graph()
G2.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
H=preproceso(G2)
## Represento gráficamente el grafo original
plt.figure(1)
nx.draw_shell(G1,with_labels=True)    
## Represento gráficamente la bipartición de los nodos 
## Imprimo por pantalla las comunidades y el autovalor asociado al autovector
print("Ejemplo")
plt.figure(2)
print(normalized_min_cut(G1))
plt.figure(3)
print(normalized_min_cut(H))
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
G7=nx.relabel_nodes(G6,mapping)
plt.figure(4)
plt.figure(figsize=(8, 6))
biparticion,valor_objetivo,colores=normalized_min_cut(G7)
print(biparticion)
nx.draw_shell(G7,with_labels=True,node_color=colores,node_size=700,font_size=18)
plt.savefig("zachary_norm_min_cut.jpg")
